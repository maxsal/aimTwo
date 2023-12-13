#' Weighted glmnet prediction models using replicate weights for selection of alpha (modified version of wlasso::wlasso())
#'
#' @description This function allows as to fit glmnet prediction (linear or
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
#' @param alpha_grid a grid for the elastic net mixing parameter alpha. Default is \code{alpha_grid=seq(0,1,by=0.05)}.
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
#' @param workers number of workers to be used in parallelization
#' @param plan_strategy strategy to be used for parallelization
#'
#' @importFrom glmnet glmnet
#' @importFrom wlasso replicate.weights error.f
#' @importFrom stats predict
#' @importFrom future plan multicore
#' @importFrom furrr future_map
#' @importFrom utils object.size
#'
#' @return A glmnet object
#' @export
wglmnet_alpha <- function(
  data          = NULL,
  outcome       = NULL,
  predictors    = NULL,
  cluster       = NULL,
  strata        = NULL,
  weights       = NULL,
  design        = NULL,
  family        = "binomial",
  lambda_grid   = NULL,
  alpha_grid    = seq(0, 1, by = 0.05),
  method        = "dCV",
  k             = 10,
  R             = 1,
  B             = 200,
  maxit         = 1E5,
  thresh        = 1E-4,
  cv_error_ind  = FALSE,
  train_prob    = 0.7,
  method_split  = c("dCV", "bootstrap", "subbootstrap"),
  seed          = NULL,
  workers       = 1,
  plan_strategy = future::multicore
) {
  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }
  if (is.null(options("future.globals.maxSize")[[1]])) {
    if (utils::object.size(data) > (450*1024^2)) {
        options(future.globals.maxSize = utils::object.size(data) + 50 * 1024^2)
    }
  } else if (options("future.globals.maxSize")[[1]] < utils::object.size(data) + 50 * 1024^2) {
    options(future.globals.maxSize = utils::object.size(data) + 50 * 1024^2)
  }

  future::plan(plan_strategy, workers = workers)
  alpha_grid |>
    furrr::future_map(
      \(x) {
        wglmnet(
          data         = data,
          outcome      = outcome,
          predictors   = predictors,
          cluster      = cluster,
          strata       = strata,
          weights      = weights,
          design       = design,
          family       = family,
          lambda_grid  = lambda_grid,
          alpha        = x,
          method       = method,
          k            = k,
          R            = R,
          B            = B,
          maxit        = maxit,
          thresh       = thresh,
          cv_error_ind = cv_error_ind,
          train_prob   = train_prob,
          method_split = method_split,
          seed         = seed
        )
      },
      .progress = TRUE
    )
}
