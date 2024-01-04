#' Helper: modified rcompanion::nagelkerke
#' @param fit A fitted model object
#' @param null A fitted null model object
#' @param restrictNobs Whether to restrict the number of observations to the number of observations in the fitted model (default: FALSE)
#' @export
nagelkerke_r2 <- function(
  fit, null = NULL, restrictNobs = FALSE
) {
  if(!restrictNobs & is.null(null)) {
    null = update(fit, ~ 1)
  }
  if(restrictNobs  & is.null(null)) {
    null = update(fit, ~ 1, data = fit$model)
  }
  if(restrictNobs  & !is.null(null)) {
    null = update(null, data = fit$model)
  }

  N = nobs(fit)

  m = suppressWarnings(logLik(fit, REML=FALSE))[1]
  n = suppressWarnings(logLik(null, REML=FALSE))[1]
  mf = 1 - m/n
  cs = 1 - exp(-2/N * (m - n))
  nk = cs/(1 - exp(2/N * n))

  data.table(
    "mcfadden_r2"   = mf,
    "cox_snell_r2"  = cs,
    "nagelkerke_r2" = nk
  )     
}
