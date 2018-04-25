#' @title Fit GW Hill
#' @description Fit a Generalised Weibull (GW) tail to the sample X and
#' @param X double(n)   data sample
#' @param p probabilities of exceedance of quantiles to be estimated
#'
#' @author Cees de Valk, KNMI
#' @export
FitGW_iHill <- function(X, p, N, EI, theta0, logdisp0, l0, XId) {
  #
  # module: FitGW_iHill.R
  # 
  # purpose: Fit a Generalised Weibull (GW) tail to the sample X and
  # (optionally) estimate quantiles for probabilities of exceedance in p
  #
  # usage:  estimates <- FitGW_iHill(X, p, N, EI= 1, theta0= NULL, logdisp0= NULL, l0= NULL, XId= '')
  #
  # X           double(n)   data sample
  # p           double(np)  (optional) probabilities of exceedance of quantiles to be estimated
  # N           integer(1)  (optional) (effective) sample size (in case X is not complete but
  #                         contains only (peak) values above some threshold)
  # EI          double(1)   (optional) extremal index (default: 1): reciprokal of cluster length in time-steps  
  #                         (e.g. estimated by ????) (EI <= 1!!!)
  # theta0      double(1)   (optional) value of theta0 in case it is imposed 
  # logdisp0    double(1)   (optional) value of log of dispersion coeff. in case it is imposed 
  #                         (dispersion coeff. is the raio of scale par. to location par.)
  # l0          integer(0)  (optional) value of l (no. of order stats used) in case it is imposed 
  # XId     character       (optional) data identifier to store with output for traceability
  # estimates   list with members 
  #   l                   no. of order statistics used for scale and quantile estimation    
  #   k                   no. of order statistics used for tail index estimation 
  #   sigma               = 1: fixed algorithm parameter (see ref. eq. (30))
  #   theta               estimates or imposed value of log-GW tail index 
  #   tailindexStd        standard deviations of tail index estimates
  #   logdisp             estimates or imposed value of log of dispersion coeff.  
  #   logdispStd          standard deviations of log of dispersion coeff. estimates
  #   scale               estimates of log-GW scale parameter
  #   locationStd         standard deviation of order statistic
  #   lambda              ratio of logarithms of probabilities of exceedance of quantile and threshold  
  #   p                   probabilities of exceedance of quantiles to be estimated 
  #   quantile            quantile estimates
  #   quantileStd         standard deviations of quantile estimates
  #   orderstats          data X sorted (decreasing)
  #   df                  = "GW": fitted distribution function tail (Generalised Weibull
  #   estimator           = "iteratedHill": see "method" below
  #
  # remark: In case quantiles are to be estimated for given frequencies mu and
  #         (a) if X contains all values (possibly above some threshold), then set 
  #               p <- mu*d/EI, 
  #               N <- T/d (if X contains all observed values over time T, then T/d = n)
  #               EI <- EI
  #         with T the total observation period, and d the time step. 
  #         Note that frequency and time step are defined with reference to the same unit of time!! 
  #         (b) if X contains only the (approx. Poisson) peak values above some threshold (so you want 
  #         to do a POT analysis), then 
  #               p <- mu*d/EI, 
  #               N <- (T/d)*EI,
  #               EI <- 1
  #         (note that d/EI is mean duration of an event, EI/d is mean number of events per unit of time)
  #
  # method: See De Valk, C. and Cai, J.J. (2018), A high quantile estimator
  #         based on the log-generalized Weibull tail limit. Econometrics and 
  #         Statistics (in press)
  #
  # Fixed parameter (controls the relationship between l and k)
  sigma2 <- 1
   
  # Handle arguments
  if (missing(p)) {p <- NULL}
  if (missing(N)) {N <- 0} 
  if (missing(EI)) {EI <- 1}
  if (missing(theta0)) {theta0 <- NULL}
  if (missing(logdisp0)) {logdisp0 <- NULL}
  if (missing(l0)) {l0 <- NULL}
  if (missing(XId)) {XId <- ''}
  X <- c(X)
  n <- length(X)  
  estimates <- NULL
  
  if (n > 0) {
    
    # Correct N if invalid or smaller than n 
    if (N < n) {N <- n}  # this is actually very wrong if POT is done
    
    # Harmonic numbers (H) and expectations of exponential order statistics (see (15)) 
    L <- 1:n 
    H <- c(0, cumsum(1/L))
    th <- sum(1/(1:N))-H
    th <- th[1:n]
    H <- H[2:(n+1)]
    
    # Defined as in proof of Theorem 2
    u <- cumsum(log(th[2:(n-1)]))/(1:(n-2))-log(th[3:n])
    
    #
    # k is found by iteration from l (see (30))
    #
    if (length(l0)== 0) {
      l <- 10:n
    } else {
      l <- l0
    }
    nl <- length(l)
    k <- rep(NA, nl)
    for (jj in 1:length(l)) {
      k[jj] <- which.min(abs((L[1:n]-1)/th^2-l[jj]/sigma2))
    }
    k <- pmax(k, l)
    k <- pmin(pmax(1, k), n+1);
    
    # Order statistics, decreasing
    X0 <- -sort(-X)
    
    # Adjust l and k based on requirements 
    ind <- which(k> 2 & X0[k]> -Inf & l> 0 & k> l*2 & k< n)
    if (length(ind)> 0) {
      k <- k[ind]
      l <- l[ind]
    }
    nl <- length(k)
    mk <- max(k)
    
    theta <- NULL
    thetaStd <- NULL
    g <- NULL
    q <- NULL
    qStd <- NULL   
    Pfluctuation <- NULL
    bias <- NULL
    logdisp <- NULL
    logdispStd <- NULL
    
    if (length(k)> 0) {
      # Hill estimator
      hill0 <- cumsum(X0[1:(mk-1)])/L[1:(mk-1)]-X0[2:mk]
      id <- (hill0== 0)
      if (any(id)) {hill0[id] <- min(hill0[!id])}
      
      if (length(theta0)== 0) {
        # The mean excess and the mean excesss of its logarithm
        hill1 <- cumsum(log(hill0[1:(mk-2)]))/L[1:(mk-2)]-log(hill0[2:(mk-1)])
        
        # Simple estimator of GW index (nondimensional curvature)
        theta <- 1+hill1[k-2]/u[k-2]
        
        # Asymptotic standard deviation of theta
        thetaStd= sqrt(sigma2/(l*EI))
        sigma2m <- sigma2 #for use in estimation of stand. dev. of quantile
      } else {
        theta <- rep(theta0[1],nl)
        thetaStd <- rep(NA,nl)
        sigma2m <- 0 #this will ensure a correct stand. dev. of the quantile
      }
      
      # Scale estimator
      if (length(logdisp0)== 0) {
        normg <- rep(0,nl)
        for (i in 1:nl) {
          normg[i] <- mean(h(theta[i],th[1:(l[i]-1)]/th[l[i]]))
        }
        g <- hill0[l-1]/normg     
        logdisp <- log(g/pmax(X0[l], .01))  # Log of dispersion coefficient
        logdispStd <- sqrt(1/(EI*l))
      } else {
        g <- X0[l]*logdisp0
        logdisp <- rep(logdisp0, nl)
        logdispStd <- rep(NA, nl)
      }
      
      # Standard deviation of X0[l] as estimator of location q(th[l])
        X0lStd <- hill0[l-1]/sqrt(l*EI)

      # Quantile estimation
      lp= length(p)
      if (lp> 0) {
        q <- matrix(NA, nl, lp)
        qStd <- q
        for (i in 1:lp) {
          lambda <- -log(p[i])/th[l] # factor of -logs of probabilities
          
          # Quantiles
          q[, i] <- X0[l]+g*h(theta, lambda)
          
          # Asymptotic standard deviations of quantiles (ADJUST!!!!)
          ha <- h(theta,lambda)
          dha <- (1/theta)*(lambda^theta*log(lambda)-ha)
          qStd[, i]= g*sqrt(ha^2+sigma2m*dha^2)/sqrt(l*EI)
        }
      }
    }
    
    #
    # compute Boucheron-Thomas (cf. Drees-Kaufmann) type simple choice of 
    # threshold l
    # 
    estimatesBT <- NULL  
    estimatesP0 <- NULL  
    P <- NULL
    bias <- NULL  

    if (length(l0)== 0) {
      if (length(theta0)== 0) { 
        vtest <- theta
        vtestStd <- thetaStd    
      } else if (length(logdisp0)== 0) {
        vtest <- logdisp
        vtestStd <- logdispStd  
      } else {
        vtest <- X0[l]
        vtestStd <- X0lStd              
      }
      
      i <- selectThresholdBT(vtest, vtestStd, l, 50)
      estimatesBT <- list("k"= k[i], "l"= l[i], "y"= th[l[i]], 
                          "N"= N, "sigma"= sqrt(sigma2), "EI"= EI,
                          "tailindex"= theta[i], "scale"= g[i], 
                          "location"= X0[l[i]], "locationStd"= X0lStd[i],
                          "logdisp"= logdisp[i], "logdispStd"= logdispStd[i],
                          "p"= p, "quantile"= q[i, ], 
                          "tailindexStd"= thetaStd[i], "quantileStd"= qStd[i, ], 
                          "df"= "GW", 
                          "estimator"= "iteratedHill", "XId"= XId)
      
      #
      # compute refined estimate including estimate of bias
      #
      rP0 <- selectThresholdP0(vtest, vtestStd, l, 50)
      i <- rP0$i
      Pfluctuation <- rP0$P
      bias <- rP0$bias
      estimatesP0 <- list("k"= k[i], "l"= l[i], "y"= th[l[i]], 
                          "N"= N, "sigma"= sqrt(sigma2), "EI"= EI,
                          "tailindex"= theta[i], "scale"= g[i], 
                          "logdisp"= logdisp[i], "logdispStd"= logdispStd[i],
                          "location"= X0[l[i]], "locationStd"= X0lStd[i],
                          "p"= p, "quantile"= q[i, ], 
                          "tailindexStd"= thetaStd[i], "quantileStd"= qStd[i, ], 
                          "df"= "GW", 
                          "estimator"= "iteratedHill", "XId"= XId)
    }
    
    # # Average weighted with inverse of variance times jump probability
    # weight= Pj/thetaStd^2
    # thetaP <- sum(theta*weight)/sum(weight) 
    # # thetaStd= sqrt(sigma2/(l*EI))
    # temp <- sqrt(sum(Pj^2*l+2*Pj*cumsum(Pj*l)))/sum(Pj*l)
    # thetaStdP <- temp*sqrt(sigma2/EI) # from B.M. statistics
    # qP <- qStdP <- p # as dummy
    # for (j in (1:lp)) {                     # same for the quantiles
    #   qP[j] <- sum(q[, j]*weight)/sum(weight) 
    #   # qStdP[j] <- sqrt(sum(qStd[, j]^2*weight)/sum(weight))
    #   # from B.M. statistics: 
    #   Q <- qStd[, j]*sqrt(l)
    #   qStdP[j] <- sqrt(sum(Pj^2*l*Q^2+2*Pj*Q*cumsum(Pj*Q*l)))/sum(Pj*l)
    # }
    # 
    # # As above, but weighted median estimates (robust)
    # thetaPP <- weightedMedian(theta,weight) 
    # thetaStdPP <- thetaStdP*sqrt(pi/2)  # maybe wrong but not very wrong
    # qPP <- qStdPP <- p                           
    # for (j in (1:lp)) {qPP[j] <- weightedMedian(q[, j],weight)} 
    # qStdPP <- qStdP*sqrt(pi/2)          # maybe wrong but not very wrong
    
    # the estimator below is suboptimal (relatively hight variance)
    # sP <- sum(P)
    # lP <- sum(P*l)/sP    
    # thetaP <- sum(P*theta)/sP
    # thetaStdP <- sqrt(sum(P*thetaStd^2)/sP)
    # Pm <- matrix(P, nl, lp)
    # sPm <- colSums(Pm)
    # qP <- colSums(Pm*q)/sPm
    # qStdP <- sqrt(colSums(Pm*qStd^2)/sPm)
    #
    # # weighted median estimates (robust)
    # thetaPP <- weightedMedian(theta,P) 
    # qPP <- qP
    # for (j in (1:lp)) {qPP[j] <- weightedMedian(q[, j],P)} 
    # lPP <- weightedMedian(l,P) 
    #
    # estimatesP <- list("l"= l[Pj> 0], "sigma"= sqrt(sigma2), "EI"= EI, "N"= N,
    #                    "tailindex"= thetaP, "tailindexStd"= thetaStdP, 
    #                    "scale"= gP, 
    #                    "p"= p, "quantile"= qP, "quantileStd"= qStdP, 
    #                    "df"= "GW", 
    #                    "estimator"= "iteratedHill", "XId"= XId)
    # 
    # estimatesPP <- list("l"= l[Pj> 0], "sigma"= sqrt(sigma2), "EI"= EI, "N"= N,
    #                    "tailindex"= thetaPP, "tailindexStd"= thetaStdPP, 
    #                    "scale"= gPP, 
    #                    "p"= p, "quantile"= qPP, "quantileStd"= qStdPP, 
    #                    "df"= "GW", 
    #                    "estimator"= "iteratedHill", "XId"= XId)
    
    # plot(l, P, type= 'l', log= 'x', ylim= c(0,1)); grid()
    
  } # if (n > 0)
  
  estimates <- list("k"= k, "l"= l, "y"= th[l], 
                    "N"= N, "sigma"= sqrt(sigma2), "EI"= EI,
                    "tailindex"= theta, "tailindexStd"= thetaStd, 
                    "scale"= g, "logdisp"= logdisp, "logdispStd"= logdispStd,
                    "location"= X0[l], "locationStd"= X0lStd,
                    "p"= p, "quantile"= q, "quantileStd"= qStd, 
                    "orderstats"= X0, "df"= "GW", 
                    "estimator"= "iteratedHill", "XId"= XId, 
                    "estimatesBT"= estimatesBT,  # Boucheron-Thomas estimate
                    "Pfluctuation"= Pfluctuation,# fluctuation size p-value
                    "bias"= bias,                # order of magnitude of bias
                    "estimatesP0"= estimatesP0)  # single-threshold-estimates (a la B-T)
                    # "estimatesP"= estimatesP,    # multi-threshold mean (jump points)
                    # "estimatesPP"= estimatesPP)  # multi-threshold median (jump points)
  
  return(estimates)
}

