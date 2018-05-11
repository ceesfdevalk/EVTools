#' @name  selectThresholdBT
#' 
#' @title selectThresholdBT
#' 
#' @description Automatic threshold selection for GW (Generalised Weibull) tail 
#'              estimation, based on fluctuation statistics of the GW tail index
#' 
#' @param tailindex    GW tail index estimates from FitGW_iHill.R (double(nl))
#' @param tailindexStd Standard deviation of tail index from FitGW_iHill.R (double(nl))
#' @param l            number of order statistics above the threshold for scale
#'                     and location estimation from FitGW_iHill.R (double(nl))
#' @param l0           (optional) lower bound of range of l to check
#' 
#' @usage Value <- selectThresholdBT(tailindex, tailindexStd, l, l0) 
#' 
#' @return i     index in the vectors tailindex and l, representing the selected threshold}    
#'
#' @author Cees de Valk \email{ceesfdevalk@gmail.com}
#' 
#' @details Works together with FitGW_iHill.M. 
#'        Based on ref. De Valk and Cai(2018), in particular on eq.(27).
#         Idea of simple estimate is from Boucheron & Thomas (2015); see also  
#         Drees and Kaufmann (1998).
#         Phi in this equation (see (22)) can be approximated by a rescaled Brownian
#         motion in the large-sample limit; then we use the BM statistics as in
#         Lemma 1 of de Valk & Cai (2018)
#'           
#' @references
#' De Valk, C. and Cai, J.J. (2018), A high quantile estimator based on 
#' the log-generalized Weibull tail limit. Econometrics and Statistics 6, 107-128, see
#' \url{https://doi.org/10.1016/j.ecosta.2017.03.001}
#' 
#' @export
selectThresholdBT <- function(theta, thetaStd, l, l0) {
  if (missing(l0)) {l0 <- min(l)}
  nl <-  length(l)
  mdif <- rep(NA, nl)
 
  j0 <- max(l0-l[1],1)
  for (j in (j0:(nl-1))) {
    lf <- sqrt(2*log(log(max(l[j+1],3))))
    mdif[j+1] <- max(abs(theta[j0:j]-theta[j+1])/(thetaStd[j0:j]*lf))
  }
  mdif[nl] <- 0
  i <- max(which(mdif < 1))
}

