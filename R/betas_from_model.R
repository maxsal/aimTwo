#' Extract model information from glm models
#' @param glm_model A fitted glm model object
#' @param intercept Whether to include the intercept in the output
#' @importFrom data.table data.table
#' @importFrom stats coef qnorm
#' @return A data.table with the following columns:
#' \itemize{
#'  \item \code{predictor}: The name of the predictor
#' \item \code{beta}: The beta coefficient
#' \item \code{se_beta}: The standard error of the beta coefficient
#' \item \code{p_value}: The p-value of the beta coefficient
#' \item \code{or}: The odds ratio
#' \item \code{or_lower}: The lower bound of the odds ratio
#' \item \code{or_upper}: The upper bound of the odds ratio
#' }
#' @export
betas_from_glm <- function(glm_model, intercept = FALSE) {
  glm_summary <- stats::coef(summary(glm_model))

  predictors <- rownames(glm_summary)
  betas <- glm_summary[, 1]
  se_betas <- glm_summary[, 2]
  p_values <- glm_summary[, 4]
  ors <- exp(betas)
  or_lowers <- exp(betas - qnorm(0.975) * se_betas)
  or_uppers <- exp(betas + qnorm(0.975) * se_betas)

  out <- data.table::data.table(
    predictor = predictors,
    beta = betas,
    se_beta = se_betas,
    p_value = p_values,
    or = ors,
    or_lower = or_lowers,
    or_upper = or_uppers
  )
  if (intercept == FALSE) out <- out[!(predictors %in% c("Intercept", "(Intercept)")), ]
  return(out)
}

#' Extract model information from logistf models
#' @param logistf_model A fitted glm model object
#' @param intercept Whether to include the intercept in the output
#' @importFrom data.table data.table
#' @importFrom stats qnorm vcov
#' @return A data.table with the following columns:
#' \itemize{
#'  \item \code{predictor}: The name of the predictor
#' \item \code{beta}: The beta coefficient
#' \item \code{se_beta}: The standard error of the beta coefficient
#' \item \code{p_value}: The p-value of the beta coefficient
#' \item \code{or}: The odds ratio
#' \item \code{or_lower}: The lower bound of the odds ratio
#' \item \code{or_upper}: The upper bound of the odds ratio
#' }
#' @export
betas_from_logistf <- function(logistf_model, intercept = FALSE) {
  predictors <- logistf_model$terms

  betas <- logistf_model$coefficients
  se_betas <- sqrt(diag(vcov(logistf_model)))
  p_values <- logistf_model$prob
  ors <- exp(betas)
  logistf_logistf_mod_confint <- suppressMessages({
    confint(logistf_model)
  })
  or_lowers <- exp(logistf_logistf_mod_confint[, 1])
  or_uppers <- exp(logistf_logistf_mod_confint[, 2])

  out <- data.table::data.table(
    predictor = predictors,
    beta = betas,
    se_beta = se_betas,
    p_value = p_values,
    or = ors,
    or_lower = or_lowers,
    or_upper = or_uppers
  )
  if (intercept == FALSE) out <- out[!(predictors %in% c("Intercept", "(Intercept)")), ]
  return(out)
}


#' Extract beta coefficients from fitted cv.glmnet models
#' @param cv.glmnet_model A fitted cv.glmnet model object
#' @param intercept Whether to include the intercept in the output
#' @param lambda The lambda value to use for the beta coefficients (default is "lambda.min")
#' @importFrom data.table data.table
#' @return A data.table with the following columns:
#' \itemize{
#'  \item \code{predictor}: The name of the predictor
#' \item \code{beta}: The beta coefficient
#' }
#' @export
betas_from_cv.glmnet <- function(cv.glmnet_model, intercept = FALSE, lambda = "lambda.min") {
  out <- predict(cv.glmnet_model, type = "coef", s = lambda) |>
    as.matrix() |>
    (\(x) {
      data.table::data.table(
        predictor = rownames(x),
        beta = x[, 1]
      )
    })()
  if (intercept == FALSE) out <- out[!(predictor %in% c("Intercept", "(Intercept)")), ]
  return(out)
}

#' Extract model information from glm or logistf models
#' @param model A fitted (svy)glm or logistf model object
#' @return A data.table with the following columns:
#' \itemize{
#'  \item \code{predictor}: The name of the predictor
#' \item \code{beta}: The beta coefficient
#' \item \code{se_beta}: The standard error of the beta coefficient
#' \item \code{p_value}: The p-value of the beta coefficient
#' \item \code{or}: The odds ratio
#' \item \code{or_lower}: The lower bound of the odds ratio
#' \item \code{or_upper}: The upper bound of the odds ratio
#' }
#' @export
betas_from_mod <- function(model, intercept = FALSE, lambda = "lambda.min") {
  model_class <- class(model)
  if ("glm" %in% model_class) {
    model_class <- "glm"
  } else if ("logistf" %in% model_class) {
    model_class <- "logistf"
  } else if ("cv.glmnet" %in% model_class) {
    model_class <- "cv.glmnet"
  } else {
    stop("only glm, logistf, and cv.glmnet models supported")
  }
  out <- switch(
    model_class,
    glm     = betas_from_glm(model, intercept = intercept),
    logistf = betas_from_logistf(model, intercept = intercept),
    cv.glmnet = betas_from_cv.glmnet(model, intercept = intercept, lambda = lambda)
  )
  out[["model"]] <- class(model)[1]
  return(out)
}

#' Calculate predicted values from a model
#' @param data_table A data.table with the data to predict on
#' @param beta_table A data.table with the beta coefficients
#' @param expit_out Whether to also return the predicted values on the probability scale
#' @param simplify Whether to return only predicted values (TRUE) or the full data.table with predicted values (FALSE)
#' @importFrom data.table data.table copy set
#' @return A data.table with the following columns:
#' \itemize{
#'  \item \code{pred}: The predicted value
#' \item \code{pred_expit}: The predicted value on the logit scale
#' }
#' @export
sum_beta_weights <- function(data_table, beta_table, expit_out = TRUE, simplify = TRUE) {
    pred_cols <- beta_table[["predictor"]]
    pred_cols <- pred_cols[!(pred_cols %in% c("Intercept", "(Intercept)"))]
    not_in_data <- pred_cols[!(pred_cols %in% colnames(data_table))]
    if (length(not_in_data) > 0) {
        message(paste0("The following predictors are not in the data: ", paste(not_in_data, collapse = ", ")))
    }
    pred_cols <- pred_cols[!(pred_cols %in% not_in_data)]
    data_out <- data.table::copy(data.table::as.data.table(data_table))

    if ("(Intercept)" %in% beta_table[["predictor"]]) {
        data_out[["pred"]] <- beta_table[beta_table[["predictor"]] == "(Intercept)", ][["beta"]]
        expit <- function(x) 1 / (1 + exp(-x))
    } else {
        data_out[["pred"]] <- 0
    }

    for (i in pred_cols) {
        data.table::set(
            data_out,
            j = "pred",
            value = data_out[["pred"]] + (beta_table[beta_table[["predictor"]] == i, ][["beta"]] * data_out[[i]])
        )
    }
    if (simplify) {
        data_out <- data_out |> dplyr::select(id, pred)
    }
    if ("(Intercept)" %in% beta_table[["predictor"]] & expit_out == TRUE) {
        data_out <- data_out |> dplyr::mutate(pred_expit = expit(pred))
    }
    data_out
}
