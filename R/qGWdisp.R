#' @name  qGWdisp
#' 
#' @title qGWdisp
#' 
#' @description As qGW, but the GW scale parameter "scale" is replaced by the logarithm of the disperson coefficient 
#'  
#' 
#' @param logdisp      logarithm of the dispersion coefficient scale/q0 (double(1))
#' @param logdispStd   (optional) standard deviation of logdisp
#' 
#' @usage Value <- QGW(p, p0, q0, logdisp, theta, q0Std, logdispStd, thetaStd)
#' 
#' @return list containing quantile(s) and optionally, their standard deviations 
#
#' @author Cees de Valk \email{ceesfdevalk@gmail.com}
#'
#' @references
#' De Valk, C. and Cai, J.J. (2018), A high quantile estimator based on 
#' the log-generalized Weibull tail limit. Econometrics and Statistics 6, 107-128, see
#' \url{https://doi.org/10.1016/j.ecosta.2017.03.001}
#' 
#' @export
qGWdisp <- function(p, p0, q0, logdisp, theta, q0Std, logdispStd, thetaStd) {
  if (missing(q0Std)) {q0Std <- NA}
  if (missing(logdispStd)) {logdispStd <- NA}
  if (missing(thetaStd)) {thetaStd <- NA}
  
  # quantile
  lambda <- log(p)/log(p0)
  ha <- h(theta,lambda)
  disp <- exp(logdisp)
  q <- q0*(1 + disp*ha)

  # derivative of ha to theta
  dha <- (1/theta)*(lambda^theta*log(lambda)-ha)
  id <- abs(theta)< .Machine$double.eps
  if (any(id)) {dha[id] <- 0.5*(log(lambda))^2}
  
  # the following asymptotic expression is fairly accurate
  # (the last term can normally be ignored but with given, precise,
  # theta and logdisp estimates, it may not be negligible)
  varq0= q0Std^2*(q/q0)^2
  varlogdisp= (ha*disp*q0)^2*logdispStd^2
  vartheta= (disp*q0)^2*dha^2*thetaStd^2
  qStd <- sqrt(varq0 + varlogdisp + vartheta)
  res <- list(quantile= q, quantileStd= qStd, varlocation= varq0, 
              varlogdisp= varlogdisp, vartailindex= vartheta)
}