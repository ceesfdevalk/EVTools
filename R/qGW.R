#' @name  qGW
#' 
#' @title qGW
#' 
#' @description quantile of the GW (Generalised Weibull) tail and (optionally) its standard deviation 
#' 
#' @param p            probability(ies) of exceedance at which to evaluate the quantile(s) (double(np))
#' @param p0           reference probability for which the quantile q0 is given (double(1))
#' @param q0           see under p0
#' @param scale        GW scale parameter at p0 (double(1))
#' @param theta        GW tail index (double(1))
#' @param q0Std        (optional) standard deviation of q0 (double(1))
#' @param scaleStd     (optional) standard deviation of scale (double(1))
#' @param thetaStd     (optional) standard deviation of theta (double(1))
#' 
#' @usage Value <- QGW(p, p0, q0, scale, theta, q0Std, scaleStd, thetaStd)
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
qGW <- function(p, p0, q0, scale, theta, q0Std, scaleStd, thetaStd) {

  if (missing(q0Std)) {q0Std <- NA}
  if (missing(scaleStd)) {scaleStd <- NA}
  if (missing(thetaStd)) {thetaStd <- NA}
  
  # quantile
  lambda <- log(p)/log(p0)
  ha <- h(theta,lambda)
  q <- q0 + scale*ha

  # derivative of ha to theta
  dha <- (1/theta)*(lambda^theta*log(lambda)-ha)
  id <- abs(theta)< .Machine$double.eps
  if (any(id)) {dha[id] <- 0.5*(log(lambda))^2}
  
  # the following asymptotic expression is fairly accurate
  # (the last term can normally be ignored but with given, precise,
  # theta and logdisp estimates, it may not be negligible)
  qStd <- sqrt(ha^2*scaleStd^2 + scale^2*dha^2*thetaStd^2 + q0Std^2)
  
  res <- list(quantile= q, quantileStd= qStd)
}