#' @name  selectThresholdP1
#' 
#' @title selectThresholdP1
#' 
#' @description Automatic threshold selection for GW (Generalised Weibull) tail 
#'              estimation, based on fluctuation statistics of the GW tail index. More efficient 
#'              than selectThresholdP0
#' 
#' @param tailindex    GW tail index estimates from FitGW_iHill.R (double(n))
#' @param tailindexStd Standard deviation of tail index from FitGW_iHill.R (double(n))
#' @param k            number of order statistics above the threshold for 
#'                     estimation of tail index from FitGW_iHill.R (double(n))
#' @param rthresh      (optional) ratio(s) of the probability of nonexceedance of fluctuation 
#'                     size to its maximum, for threshold(s) to be accepted. Default is 0.5 (double(lr))            
#' @param kmin         tailindex at k< kmin will be skipped from the analysis (double(1))
#' @usage Value <- selectThresholdP1(tailindex, tailindexStd, k, rthresh, kmin) 
#' 
#' @return list containing the elements
#'   \item{i}{index(indices) in the array P (see below), representing the selected threshold(s)}  
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
selectThresholdP1 <- function(theta, thetaStd, k, rthresh, kmin) {
  if (missing(rthresh)) {rthresh <- 0.5}
  if (missing(kmin)) {kmin <- min(k)}  
  # parameter: fluctuation probability threshold
  # l <- thetaStd^(-2) # overwrite
  
  # the highly variable estimates of the index are excluded (not of interest, and noise spoils
  # the statistics of the fluctuations)
  ind <- k>= kmin  # differenced estimates above this will be taken into account
  k <- k[ind]
  theta <- theta[ind]
  thetaStd <- thetaStd[ind]
  # then k needs to be normalized in order to fit into the framework of Darling-Erdos 
  # The normalization is based on the stationarity of the Ornstein-Uhlenbeck process
  l <- k/min(k)
  nl <-  length(l)
  
  # This is where we compute the fluctuation statistic (to save
  # time, it is not done everywhere)
  id <- unique(round(exp((0:1.e4)*log(nl)*1.e-4)))
  id <- id[id< nl & id>1 & l[id]>3]   # make sure l[id]>3
  kid <- k[id]
  lid <- l[id]
  ll <- log(pmax(lid,3))
  lf <- log(ll)
  b <- sqrt(2*lf)
  nid <- length(id)
  alpha <- bias <- rep(0, nid) 
  
  pb <- txtProgressBar(1, nid)
  # for (j in (1:(nl-1))) {
  for (jj in (1:nid)) {   
    j= id[jj]-1
    # alpha[jj] <- max((abs(theta[1:j]-theta[j+1]))/thetaStd[1:j])/slf[jj]
    # following is a bit conservative (not quite till the end; difference does not matter in the limit)
    # alpha[jj] <- max((abs(theta[1:j]-theta[j+1])+thetaStd[j+1]*2)/thetaStd[1:j])
    # alpha[jj] <- max((abs(theta[1:j]-theta[j+1]))/thetaStd[1:j])
    bias[jj] <- max(abs(theta[1:j]-theta[j+1])-thetaStd[1:j]*b[jj])
    setTxtProgressBar(pb, j)
  }
  bias <- pmax(bias,0)                # this is unbiased if nonzero! not soft-clipping
  # a <- ll^(2*(alpha-1))*sqrt(4*pi/lf)
  # a <- ll^(2*(alpha-1))*sqrt(pi/lf)   # this is for abs. value; see Darling & Erdos
  a <- 2*exp(b*alpha)*sqrt(2*pi)/ll^2/b
  P <- 1-exp(-2/a)                    # probability of being outside alpha
  # P[nl] <- 0                        # otherwise threshold chouce may not be defined
  
  # Pc <- rev(cummax(rev(P)))           # cumul. maximum of P from the right
  # Pc[length(Pc)] <- 0                 # in case Pc are all 1
  # Pj <- -c(diff(Pc),0)                # probability jumps, upward fron the right
  # # A big jump indicates likely good threshold where 
  # # fluctuation statistics match asymptotic theory
  # Pj <- Pj*(Pc> 0.1)
  # 
  # threshold choice a la Boucheron-Thomas, but more delicate
  lr <- length(rthresh)
  i <- rep(NA, lr)
  for (j in 1:lr) {
    i[j] <- max(which(P> max(P)*rthresh[j]))
  }
  # threshold selection based on max of P
  # i <- max(which(P> quantile(P, 0.99)*rthresh)) # threshold selection based on high quantile of P
  res <- list("i"= i, "k"= kid, "P"= P, "bias"= bias)
}

