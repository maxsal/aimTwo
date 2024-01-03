#' Helper: modified rcompanion::nagelkerke
#' @param fit A fitted model object
#' @param null A fitted null model object
#' @param restrictNobs Whether to restrict the number of observations to the number of observations in the fitted model (default: FALSE)
#' @export
nagelkerke_r2 <- function(fit, null=NULL, restrictNobs=FALSE)
{
   TOGGLE =   (class(fit)[1]=="lm"
             | class(fit)[1]=="gls"
             | class(fit)[1]=="lme"
             | class(fit)[1]=="glm"
             | class(fit)[1]=="negbin"
             | class(fit)[1]=="zeroinfl"
             | class(fit)[1]=="clm"
             | class(fit)[1]=="vglm"
             | class(fit)[1]=="betareg"
             | class(fit)[1]=="rq"
             | class(fit)[1]=="glmmTMB")
   NOGGLE = is.null(null)
   ERROR  = "Note: For models fit with REML, these statistics are based on refitting with ML"
   ERROR2 = "None"
   
  if(!restrictNobs & NOGGLE  & TOGGLE){null = update(fit, ~ 1)}
  if(restrictNobs  & NOGGLE  & TOGGLE){null = update(fit, ~ 1, data=fit$model)}
  if(restrictNobs  & !NOGGLE){null = update(null, data=fit$model)}
    
  Y = matrix(rep(NA,2),
            ncol=1)
  colnames(Y) = ""
  rownames(Y) = c("Model:", "Null:")
  
  Z = matrix(rep(NA, 3),
             ncol=1)
  colnames(Z) = c("Pseudo.R.squared")
  rownames(Z) = c("McFadden", "Cox and Snell (ML)", 
                  "Nagelkerke (Cragg and Uhler)")
  
  X = matrix(rep(NA,4),
             ncol=4)
  colnames(X) = c("Df.diff","LogLik.diff","Chisq","p.value")
  rownames(X) = ""
  
  U = matrix(rep(NA,2),
            ncol=1)
  colnames(U) = ""
  rownames(U) = c("Model:", "Null:")
  
  Y[1]= toString(fit$call)

  if(TOGGLE | (!NOGGLE)){
  Y[2]= toString(null$call)
  N = nobs(fit)
  U[1,1]= nobs(fit); U[2,1]= nobs(null)
  }

  if(U[1,1] != U[2,1]){
    ERROR2 = "WARNING: Fitted and null models have different numbers of observations"}
  
  m = suppressWarnings(logLik(fit, REML=FALSE))[1]
  n = suppressWarnings(logLik(null, REML=FALSE))[1]
  mf = 1 - m/n
  Z[1,] = signif(mf, digits=6)
  cs = 1 - exp(-2/N * (m - n))
  Z[2,] = signif(cs, digits=6)
  nk = cs/(1 - exp(2/N * n))
  Z[3,] = signif(nk, digits=6)
  
  o = n - m
  dfm = attr(logLik(fit),"df")
  dfn = attr(logLik(null),"df")
  if(class(fit)[1]=="vglm"){dfm=df.residual(fit)}
  if(class(fit)[1]=="vglm"){dfn=df.residual(null)}
  dff = dfn - dfm
  CHI = 2 * (m - n)
  P = pchisq(CHI, abs(dff), lower.tail = FALSE)
  
  X [1,1] = dff
  X [1,2] = signif(o, digits=5)             
  X [1,3] = signif(CHI, digits=5)
  X [1,4] = signif(P, digits=5)     
  
  W=ERROR
  
  WW=ERROR2
  
  V = list(Y, Z, X, U, W, WW) 
  names(V) = c("Models", "Pseudo.R.squared.for.model.vs.null", 
               "Likelihood.ratio.test", "Number.of.observations",
               "Messages", "Warnings")
  return(V)            
}