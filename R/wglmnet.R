#' Weighted glmnet prediction models using replicate weights (modified version of wlasso::wlasso())
#'
#'@description This function allows as to fit glmnet prediction (linear or 
#' logistic lasso, ridge, or elastic net) models to complex survey data, 
#' considering sampling weights in the estimation process and selects the 
#' lambda that minimizes the error based on different replicating weights 
#' methods.
#'
#' @param data A data frame with information about the outcome variable and predictors and sampling weights.
#' @param outcome column of the outcome variable. It could be indicated as the column number or the column name.
#' @param predictors vector of columns indicating the covariates. It could be defined by numbers or names of the columns.
#' @param cluster name of the column indicating clusters
#' @param strata name of the column indicating strata
#' @param weights name of the column indicating sampling weights
#' @param design an object of class \code{survey.design}. It could be \code{NULL} if information about \code{cluster}, \code{strata}, \code{weights} and \code{data} are given.
#' @param family family to fit glmnet models, choose between \code{gaussian} or \code{binomial} (default).
#' @param lambda_grid a grid for penalization parameters. If it is not defined by the user, the function will define one.
#' @param alpha the elastic net mixing parameter, with \code{alpha=0} equivalent to ridge regression and \code{alpha=1} equivalent to lasso. Default is \code{alpha=1}.
#' @param method method to be applied to define replicate weights, to choose between one of these: \code{JKn}, \code{dCV}, \code{bootstrap} (default), \code{subbootstrap}, \code{BRR}, \code{split}, \code{extrapolation}
#' @param k number of folds to be defined (only for \code{cv}). Default is \code{k=10}.
#' @param R number of times the sample is partioned (needed and used only for \code{cv}, \code{split} or \code{extrapolation} methods). Default R=1.
#' @param B number of bootstrap resamples (only for \code{bootstrap} and \code{subbootstrap} methods). Default \code{B=200}.
#' @param maxit maximum number of iterations allowed. Default \code{maxit=1E5}.
#' @param thresh convergence threshold for coordinate descent. Default \code{thresh=1E-7}.
#' @param cv_error_ind method for estimating the error for \code{cv} method. FALSE (default) estimates the error for each test set and defines the cross-validated error as the average of all those errors. Option TRUE estimates the cross-validated error as the weighted average of the loss for each unit
#' @param train_prob probability for defining training sets (only for \code{split} and \code{extrapolation} methods)
#' @param method_split \code{cv} or \code{bootstrap} (only for \code{split} method)
#' @param seed define a seed
#' 
#' @importFrom glmnet glmnet
#' @importFrom wlasso replicate.weights error.f
#' @importFrom stats predict
#'
#' @return A glmnet object
#' @export
wglmnet <- function(
  data         = NULL,
  outcome      = NULL,
  predictors   = NULL,
  cluster      = NULL,
  strata       = NULL,
  weights      = NULL,
  design       = NULL,
  family       = "binomial",
  lambda_grid  = NULL,
  alpha        = 1,
  method       = "bootstrap",
  k            = 10,
  R            = 1,
  B            = 200,
  maxit        = 1E5,
  thresh       = 1E-7,
  cv_error_ind = FALSE,
  train_prob   = 0.7,
  method_split = c("dCV", "bootstrap", "subbootstrap"),
  seed         = NULL
) {
  # check if data is a data.frame and if not, convert it
  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }

  # Step 1: Generate replicate weights based on the method
  new_data <- replicate_weights(
    data         = data,
    method       = method,
    cluster      = cluster,
    strata       = strata,
    weights      = weights,
    k            = k,
    R            = R,
    B            = B,
    train_prob   = train_prob,
    method_split = method_split,
    seed         = seed
  )

  # Step 2: if is.null(lambda_grid), then initialize it
  if (is.null(lambda_grid)) {
    model_orig <- glmnet::glmnet(
      y       = as.numeric(new_data[, outcome]),
      x       = as.matrix(new_data[, predictors]),
      weights = as.numeric(new_data[, weights]),
      family  = family,
      alpha   = alpha,
      maxit   = maxit,
      thresh  = thresh
    )
    lambda_grid <- model_orig$lambda
  } else {
    model_orig <- NULL
  }

  # Step 3: Fit the training models and estimate yhat for units in the sample
  replicate_weights_training_columns <- grep("_train", colnames(new_data))
  l_yhat <- list()

  for(weight_col in replicate_weights_training_columns) {

    model <- glmnet::glmnet(
      y       = as.numeric(new_data[, outcome]),
      x       = as.matrix(new_data[, predictors]),
      weights = as.numeric(new_data[, weight_col]),
      lambda  = lambda_grid,
      alpha   = alpha,
      family  = family,
      maxit   = maxit,
      thresh  = thresh
    )

    # Sample yhat
    yhat <- stats::predict(
      model,
      newx = as.matrix(new_data[, predictors]),
      type = "response"
    )
    l_yhat[[length(l_yhat) + 1]] <- yhat
    names(l_yhat)[[length(l_yhat)]] <- paste0("yhat_", colnames(new_data)[weight_col])
  }

  # Step 4: estimate the error in the test sets
  error <- error_f(
    data         = new_data,
    l_yhat       = l_yhat,
    method       = method,
    cv_error_ind = cv_error_ind,
    R            = R,
    k            = k,
    B            = B,
    outcome      = outcome,
    family       = family,
    weights      = weights
  )

  mean_error <- apply(error, 2, mean)
  lambda_min <- lambda_grid[which.min(mean_error)]

  model <- glmnet::glmnet(
    y       = data[, outcome],
    x       = as.matrix(data[, predictors]),
    weights = data[, weights],
    lambda  = lambda_min,
    alpha   = alpha,
    family  = family,
    maxit   = maxit,
    thresh  = thresh
  )

  result <- list(
    lambda_grid   = lambda_grid,
    lambda_min    = lambda_min,
    error_min     = min(mean_error),
    average_error = mean_error,
    all_error     = error,
    model_min     = model
  )
  if (!is.null(model_orig)) {
    result$model_grid <- model_orig
  }

  return(result)
}
