selectThresholdBT <- function(theta, thetaStd, l, l0) {
  #
  # module: selectThresholdBT.R
  #
  # purpose: Select a good threshold for estimation of (log)-GW tail index
  #
  # usage:  i <- selectThresholdBT(theta, thetaStd, l, l0= 50)
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

