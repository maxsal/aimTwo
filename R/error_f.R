#' Error function
#' @param data Data frame with the data
#' @param l_yhat List of matrices with the predicted values
#' @param method Method used to generate the replicate weights
#' @param cv_error_ind Logical indicating if the cross-validation error should be calculated
#' @param R Number of replicate weights
#' @param k Number of folds
#' @param B Number of bootstrap samples
#' @param outcome Name of the outcome variable
#' @param family Family of the outcome variable
#' @param weights Name of the weights variable
#' @return Matrix with the error for each lambda
#' @export
error_f <- function(
  data,
  l_yhat,
  method,
  cv_error_ind = c(TRUE, FALSE),
  R = NULL,
  k = NULL,
  B = NULL,
  outcome,
  family,
  weights = NULL
) {

  if (method == "JKn") {
    R            <- 1
    k            <- length(grep("rw", colnames(data)))
    initial_name <- "sw"
    cv_error_ind <- TRUE
  }

  if (method %in% c("bootstrap", "subbootstrap")) {
    R            <- 1
    k            <- B
    weights_name <- weights
    cv_error_ind <- FALSE
  }

  if (method == "dCV") {
    if (cv_error_ind) {
      initial_name <- "sw"
    } else {
      initial_name <- "rw"
    }
  }

  if (method == "BRR") {
    R            <- 1
    k            <- length(l_yhat)
    initial_name <- "rw"
    cv_error_ind <- FALSE
  }

  if (method %in% c("split", "extrapolation")) {
    k            <- 1
    initial_name <- "rw"
    cv_error_ind <- FALSE
  }

  if (!cv_error_ind) {
    error_lambda_r           <- matrix(NA, nrow = R * k, ncol = dim(l_yhat[[1]])[2])
    rownames(error_lambda_r) <- 1:(R * k)
  }

  l_loss <- l_loss_w <- list()

  for (r in 1:R) {
    for (kk in 1:k) {
      # Define the names of the weights' columns
      if (method %in% c("JKn", "dCV", "BRR", "bootstrap", "subbootstrap")) {
        yhat_name <- paste0("yhat_rw_r_", r, "_train_", kk)
      }

      if (method %in% c("JKn", "dCV", "BRR")) {
        weights_name <- paste0(initial_name, "_r_", r, "_test_", kk)
      }

      if (method %in% c("split", "extrapolation")) {
        yhat_name    <- paste0("yhat_rw_r_", r, "_train")
        weights_name <- paste0(initial_name, "_r_", r, "_test")
      }

      # Calculate the loss and the weighted loss
      l_loss[[yhat_name]]   <- apply(l_yhat[[yhat_name]], 2, loss_f, y = as.numeric(data[,outcome]), family = family)
      l_loss_w[[yhat_name]] <- apply(l_loss[[yhat_name]], 2, function(x) {x * data[, weights_name] })

      # Calculate the error
      if(!cv_error_ind){
        sumwi_k                                    <- sum(data[, weights_name])
        error_lambda_r[(r-1)*k + kk,]              <- apply(l_loss_w[[yhat_name]], 2, sum) / sumwi_k
        rownames(error_lambda_r)[(r - 1) * k + kk] <- paste0(method, "_r_", r, "_k_", kk)
      }
    }
  }

  if (cv_error_ind) {
    error_lambda_r <- matrix(NA, nrow = R, ncol = dim(l_yhat[[1]])[2])

    for (r in 1:R) {
      routcomehat     <- grep(paste0("r_", r, "_"), names(l_loss_w))
      rcol_sw_newdata <- grep(paste0(initial_name,"_r_",r,"_test"), names(data))

      wi_li_r <- matrix(NA, nrow = length(routcomehat), ncol = dim(l_yhat[[1]])[2])
      for (ind in 1:length(routcomehat)) {
        wi_li_r[ind,] <- apply(l_loss_w[[routcomehat[ind]]], 2, sum)
      }
      sum_wi_li_r <- apply(wi_li_r, 2, sum)

      sum_wi_r            <- sum(data[, rcol_sw_newdata])
      error_lambda_r[r, ] <- sum_wi_li_r / sum_wi_r
    }
  }

  return(error_lambda_r)
}

#' Loss function
#' @param y_est Estimated values
#' @param y Observed values
#' @param family Family of the outcome variable
#' @return Vector with the loss for each observation
#' @export
loss_f <- function(
  y_est,
  y,
  family = "binomial"
) {

  if (family == "gaussian") {
    l <- (y_est - y)^2
  }

  if (family == "binomial") {
    l <- rep(NA, length(y))
    l[which(y==1)] <- -log(y_est[which(y==1)])
    l[which(y==0)] <- -log(1-y_est[which(y==0)])
  }

  return(l)
}
