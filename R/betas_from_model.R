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
betas_from_mod <- function(model, intercept = FALSE) {
  model_class <- class(model)
  if ("glm" %in% model_class) {
    model_class <- "glm"
  } else if ("logistf" %in% model_class) {
    model_class <- "logistf"
  } else {
    stop("only glm and logistf models supported")
  }
  switch(
    model_class,
    glm = betas_from_glm(model, intercept = intercept),
    logistf = betas_from_logistf(model, intercept = intercept)
  )
}