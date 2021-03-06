#' @name  selectThresholdP0
#' 
#' @title selectThresholdP0
#' 
#' @description Automatic threshold selection for GW (Generalised Weibull) tail 
#'              estimation, based on fluctuation statistics of the GW tail index
#' 
#' @param tailindex    GW tail index estimates from FitGW_iHill.R (double(nl))
#' @param tailindexStd Standard deviation of tail index from FitGW_iHill.R (double(nl))
#' @param k            number of order statistics above the threshold for 
#'                     estimation of tail index from FitGW_iHill.R (double(nl))
#' @param rthresh      (optional) ratio of the probability of nonexceedance of fluctuation 
#'                     size to its maximum, for threshold to be accepted. Default is 0.5.             
#'  
#' @usage Value <- selectThresholdP0(tailindex, tailindexStd, k, rthresh) 
#' 
#' @return list containing the elements
#'   \item{i}{index in the vectors tailindex and l, representing the selected threshold}  
#'   \item{k}{the elements of k for which P and bias are computed}      
#'   \item{P}{p-values of observed fluctuation statistic}    
#'   \item{bias}{estimate of bias in tail index based on fluctuation statistic}    
#'
#' @author Cees de Valk \email{ceesfdevalk@gmail.com}
#' 
#' @details Works together with FitGW_iHill.R and similar methods.  
#'        Based on ref. De Valk and Cai(2018), in particular on eq.(27).
#'        Idea of simple estimate is from Boucheron & Thomas (2015); see also  
#'        Drees and Kaufmann (1998).
#'        Phi in this equation (see (22)) can be approximated by a rescaled Brownian
#'        motion in the large-sample limit; then we use the BM statistics as in
#'        Lemma 1 of de Valk & Cai (2018).
#'        Serial dependence is accounted for in a simplistic way implicitly through 
#'        tailindexStd.
#'           
#' @references
#' De Valk, C. and Cai, J.J. (2018), A high quantile estimator based on 
#' the log-generalized Weibull tail limit. Econometrics and Statistics 6, 107-128.
#' Boucheron, S., Thomas, M. (2015), Tail index estimation, concentration and adaptivity. 
#' Electron. J. Stat. 9, 2751–2792.
#' Drees, H., Kaufmann, E. (1998), Selecting the optimal sample fraction in univariate 
#' extreme value estimation. Stoch. Process. Appl. 75, 149–172.
#' 
#' @export
selectThresholdP0 <- function(theta, thetaStd, k, rthresh) {
  if (missing(rthresh)) {rthresh <- 0.5}
  # parameter: fluctuation probability threshold
  # l <- thetaStd^(-2) # overwrite
  
  # the highly variable estimates of the index are excluded (not of interest, and noise spoils
  # the statistics of the fluctuations)
  ind <- thetaStd< 1
  k <- k[ind]
  theta <- theta[ind]
  thetaStd <- thetaStd[ind]
  # then k needs to be normalized in order to fit into the framework of Darling-Erdos 
  # The normalization is based on the stationarity of the Ornstein-Uhlenbeck process
  l <- k/min(k)
  
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
  # threshold choice a la Boucheron-Thomas, but more delicate
  i <- max(which(P> max(P)*rthresh))              # threshold selection based on max of P
  # i <- max(which(P> quantile(P, 0.99)*rthresh)) # threshold selection based on high quantile of P
  res <- list("i"= i, "k"= k, "P"= P, "bias"= bias)
}

