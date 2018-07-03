#' @name  qckGW
#' 
#' @title qckGW
#' 
#' @description quantile of the GW (Generalised Weibull) tail and (optionally) its standard deviation 
#' 
#' @param p            probability(ies) of exceedance at which to evaluate the quantile(s) (double(np))
#' @param es           fitted GW model (see e.g. FitGW_iHill.R)
#' @param l            (optional) threshold rank: no. of order statistics exceeding the threshold for 
#'                     the scale parameter (see ref. for notation) 
#'
#' @usage Value <- qckGW(p, es, l)
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
qGW <- function(p, es, l) {
 if (missing(l)) {
    if (length(es$l)> 1) {
      stop('Argument l must be supplied: es contains values for 
           multiple thresholds.')
    } else if (is.numeric(es$l) & es$l> 0) {
      l <- es$l> 0
    } else {
      stop('es$l not a positive number.')
    }
  i0 <- which(es$l== l) 
  if (length(i0)< 1) {
    stop('l is not in es$l.')
  } else {
    y <- es$y
    q0 <- es$q0
    scale <- es$scale
    theta <- es$theta
    q0Std <- es$q0Std
    scaleStd <- es$scaleStd
    thetaStd <- es$thetaStd
  }
  
  # quantile
  lambda <- -log(p)/es$y
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
  
  