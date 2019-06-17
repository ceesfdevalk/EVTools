#' @name  selectThresholdP0
#' 
#' @title selectThresholdP0
#' 
#' @description Automatic threshold selection for GW (Generalised Weibull) tail 
#'              estimation, based on fluctuation statistics of the GW tail index
#' 
#' @param tailindex    GW tail index estimates from FitGW_iHill.R (double(nl))
#' @param tailindexStd Standard deviation of tail index from FitGW_iHill.R (double(nl))
#' @param l            number of order statistics above the threshold for scale
#'                     and location estimation from FitGW_iHill.R (double(nl))
#' @param rthresh      (optional) ratio of the probability of nonexceedance of fluctuation 
#'                     size to its maximum, for threshold to be accepted. Default is 0.4.             
#'  
#' @usage Value <- selectThresholdP0(tailindex, tailindexStd, l, rthresh) 
#' 
#' @return list containing the elements
#'   \item{i}{index in the vectors tailindex and l, representing the selected threshold}    
#'   \item{P}{p-values of observed fluctuation statistic}    
#'   \item{bias}{estimate of bias in tail index based on fluctuation statistic}    
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
#         Serial dependence is accounted for in a simplistic way using the EI.
#'           
#' @references
#' De Valk, C. and Cai, J.J. (2018), A high quantile estimator based on 
#' the log-generalized Weibull tail limit. Econometrics and Statistics 6, 107-128, see
#' \url{https://doi.org/10.1016/j.ecosta.2017.03.001}
#' 
#' @export
selectThresholdP0 <- function(theta, thetaStd, l, rthresh) {
  if (missing(rthresh)) {rthresh <- 0.}
  # parameter: fluctuation probability threshold
  # l <- thetaStd^(-2) # overwrite
  ind <- thetaStd< 1
  l <- l[ind]
  l <- l/min(l)
  theta <- theta[ind]
  thetaStd <- thetaStd[ind]
  
  nl <-  length(l)
  ll <- log(pmax(l,3))
  lf <- log(ll)
  slf <- sqrt(2*lf)
  
  alpha <- bias <- rep(0, nl) 
  
  pb <- txtProgressBar(1, nl)
  for (j in (1:(nl-1))) {
    alpha[j+1] <- max(abs(theta[1:j]-theta[j+1])/thetaStd[1:j])/slf[j+1]
    bias[j+1] <- max(abs(theta[1:j]-theta[j+1])-thetaStd[1:j]*slf[j+1])
    setTxtProgressBar(pb, j)
  }
  bias <- pmax(bias,0)                # this is unbiased if nonzero! not soft-clipping
  # a <- ll^(2*(alpha-1))*sqrt(4*pi/lf)
  a <- ll^(2*(alpha-1))*sqrt(pi/lf)   # this is for abs. value; see Darling & Erdos
  P <- 1-exp(-1/a)                    # probability of being outside alpha
  P[nl] <- 0                          # otherwise threshold chouce may not be defined
  
  # Pc <- rev(cummax(rev(P)))           # cumul. maximum of P from the right
  # Pc[length(Pc)] <- 0                 # in case Pc are all 1
  # Pj <- -c(diff(Pc),0)                # probability jumps, upward fron the right
  # # A big jump indicates likely good threshold where 
  # # fluctuation statistics match asymptotic theory
  # Pj <- Pj*(Pc> 0.1)
  # 
  # threshold choice a la Boucheron-Thomas (but more delicate)
  i <- max(which(P> max(P)*rthresh))              # threshold selection based on max of P
  # i <- max(which(P> quantile(P, 0.99)*rthresh)) # threshold selection based on high quantile of P
  res <- list("i"= i, "P"= P, "bias"= bias)
}

