selectThresholdP0 <- function(theta, thetaStd, l, l0) {
  #
  # module: selectThresholdP0.R
  #
  # purpose: Select a good threshold for estimation of (log)-GW tail index
  #
  # usage:  i <- selectThresholdP0(theta, thetaStd, l, l0= 50)
  #
  # theta       double(nl)  (log)-GW tail index estimates
  # thetaStd    double(nl)  standard dev. of (log)-GW tail index estimates
  # l           integer(nl) no. of upper order statistics used to estimate scale 
  #                         (and indirectly, also the tail index)
  # l0          integer     lower bound of range of l to check
  # output: 
  #    i            index into the vectors l and theta (and all similar output
  #                 vectors of FitGW.M): l(i) is the selected no. of upper 
  #                 order statistics
  #
  # remark:     works together with FitGW_iHill.M
  #
  # method: Based on De Valk, C. and Cai, J.J. (2018), A high quantile estimator
  #         based on the log-generalized Weibull tail limit. Econometrics and 
  #         Statistics (in press). In particular: on eq.(27)
  #         Idea of simple estimate is from Boucheron & Thomas (2015); see also  
  #         Drees and Kaufmann (1998)
  #         Phi in this equation (see (22)) can be approximated by a rescaled Brownian
  #         motion (in the large-sample limit); then we use the BM statistics as in
  #         Lemma 1 of de Valk & Cai (2018)
  #         Serial dependence is accounted for in a simplistic way using the EI.
  #
  pthresh <- 0.3 # parameter: fluctuation probability threshold

  nl <-  length(l)
  ll <- log(pmax(l,3))
  lf <- log(ll)
  slf <- sqrt(2*lf)
  
  alpha <- bias <- rep(0, nl) 
  
  for (j in (1:(nl-1))) {
    alpha[j+1] <- max(abs(theta[1:j]-theta[j+1])/thetaStd[1:j])/slf[j+1]
    bias[j+1] <- max(abs(theta[1:j]-theta[j+1])-thetaStd[1:j]*slf[j+1])
  }
  bias <- pmax(bias,0)                # this is unbiased if nonzero! not soft-clipping
  a <- ll^(2*(alpha-1))*sqrt((4*pi)/lf)
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
  # i <- max(which(P> max(P)*pthresh))        # threshold selection based on max of P
  i <- max(which(P> quantile(P, 0.99)*pthresh))        # threshold selection based on high quantile of P
  res <- list("i"= i, "P"= P, "bias"= bias)
}

