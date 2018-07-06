#' @name  FitGW_ADGG
#' 
#' @title FitGW_ADGG
#' 
#' @description Fit a Generalised Weibull (GW) upper tail to the sample X and estimate quantiles, using
#' the moment estimator of Albert et al. (2018) 
#' 
#' @param X data sample (double(n))
#' @param p probabilities of exceedance of the quantiles to be estimated (double(np))  
#' @param N (optional) (effective) sample size, in case X is not complete but contains only (peak) values above some threshold (integer(1))
#' @param r11 (optional) factor to increase estimator variance by, to account for serial dependence (default: 1) (double(1))
#' @param fixedpar (optional): fixed model parameters not to be estimated, and their standard errors (list; see Details)
#' @param l0 (optional) value of l (no. of order stats used) in case it is imposed (integer(0))
#' @param XId (optional) data identifier to store with output for traceability (character)
#' 
#' @usage Value <- FitGW_iHill(X, p, N= 0, r11= 1, fixedpar= NULL, l0= NULL, XId= '')
#' 
#' @return A list, with members: 
#'   \item{l}{no. of order statistics used for scale and quantile estimation}    
#'   \item{k}{no. of order statistics used for tail index estimation} 
#'   \item{sigma}{= 1: fixed algorithm parameter (see ref. eq. (30))}
#'   \item{tailindex}{estimates or imposed value of log-GW tail index} 
#'   \item{tailindexStd}{standard deviations of tail index estimates}
#'   \item{logdisp}{estimates or imposed value of log of dispersion coeff.}  
#'   \item{logdispStd}{standard deviations of log of dispersion coeff. estimates}
#'   \item{scale}{estimates of log-GW scale parameter}
#'   \item{locationStd}{standard deviation of order statistic}
#'   \item{lambda}{ratio of logarithms of probabilities of exceedance of quantile and threshold}  
#'   \item{p}{probabilities of exceedance of quantiles to be estimated} 
#'   \item{quantile}{quantile estimates}
#'   \item{quantileStd}{standard deviations of quantile estimates}
#'   \item{orderstats}{data X sorted (decreasing)}
#'   \item{df}{= "GW": fitted distribution function tail (Generalised Weibull}
#'   \item{estimator}{= "iteratedHill": see "method" below}
#' 
#' @details
#'  
#'  This version differs from FitGW_iHill: for small values of the index, it is less biased.
#'  But it performs less well if the index equals one.
#'  
#'  Pre-determined model parameters are to be supplied in the list fixedpar (see above):
#'  \itemize{
#'   \item{$theta0: (optional) value of tailindex in case it is imposed (double(1))}
#'   \item{$theta0Std: (optional) its standard deviation (double(1))}
#'   \item{$logdisp0: (optional) value of log of dispersion coeff. in case it is imposed (dispersion coeff. is the raio of scale par. to location par.) (double(1))}
#'   \item{$logdisp0Std: (optional) its standard deviation (double(1))}        
#'   }
#'   
#'   In case a quantile is to be estimated for a \emph{frequency}, say f, and 
#'   \enumerate{
#'   \item{if X contains all values (possibly above some threshold), then with
#'   EI an estimate of the Extremal Index, set
#'   p = f*d/EI and N = T/d, with T the length of the observation period and d the time step. 
#'         Note that f and d are defined with reference to the same unit of time!! In this case,
#'         r11 needs to be estimated.
#'             }
#'   \item{if X contains only the (approximately Poisson) peak values above some threshold 
#'         (i.e., you want to do a PoT analysis), then set p = f*d/EI, N = (T/d)*EI, r11= 1 
#'         (note that d/EI is the mean duration of an "event" and EI/d is the 
#'         mean number of "events" per unit of time).
#'              }
#' }   
#'           
#' @references
#' De Valk, C. and Cai, J.J. (2018), A high quantile estimator based on 
#' the log-generalized Weibull tail limit. Econometrics and Statistics 6, 107-128, see
#' \url{https://doi.org/10.1016/j.ecosta.2017.03.001}
#
#' @author Cees de Valk \email{ceesfdevalk@gmail.com}
#' 
#' @export
FitGW_iHill <- function(X, p, N, r11, fixedpar, l0, XId) {
  # fixed parameter 
  sigma2 <- 0
  
  # Handle arguments
  if (missing(p)) {p <- NULL}
  if (missing(N)) {N <- 0} 
  if (missing(r11)) {r11 <- 1}
  if (missing(fixedpar)) {fixedpar <- NULL}
  if (missing(l0)) {l0 <- NULL}
  if (missing(XId)) {XId <- ''}
  
  theta0 <- fixedpar$theta0
  theta0Std <- fixedpar$theta0Std
  logdisp0 <- fixedpar$logdisp0
  logdisp0Std <- fixedpar$logdisp0Std
  
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
    
    #
    # k is found by iteration from l (see (30))
    #
    if (length(l0)== 0) {
      l <- 10:n
    } else {
      l <- l0
    }
    nl <- length(l)
    k <- l  # Newton iteration for k if sigma2> 0
    if (sigma2> 0){
      for (jj in 1:10) {
        k <- k-(k-l*th[k]^2/sigma2)/(1+2/th[k])
      }
    }
    k <- pmax(k, l)
    k <- pmin(pmax(1, k), n+1);
    
    # Order statistics of logarithms, decreasing
    X <- -sort(-X)
    X0 <- log(pmax(X,0))
    
    # Adjust l and k based on requirements 
    ind <- which(k> 2 & X0[k]> -Inf & l> 0 & k< n)
    if (length(ind)> 0) {
      k <- k[ind]
      l <- l[ind]
    }
    nl <- length(k)
    mk <- max(k)
    ml <- max(l)
    
    theta <- NULL
    thetaStd <- NULL
    g <- NULL
    q <- NULL
    qStd <- NULL   
    Pfluctuation <- NULL
    bias <- NULL
    logdisp <- NULL
    logdispStd <- NULL
    
    thetagrid <- seq(-1.49, 1.49, .02)  # not zero, for convenience
    lg <- length(thetagrid)
    
    if (length(k)> 0) {
      # Hill estimator
      mom1 <- cumsum(X0[1:(mk-1)])/L[1:(mk-1)]-X0[2:mk]
      mom2 <- cumsum(X0[1:(mk-1)]^2)/L[1:(mk-1)]-X0[2:mk]^2
      id <- pmin(mom1, mom2)<= 0
      if (any(id)) {
        mom1[id] <- min(mom1[!id])
        mom2[id] <- min(mom2[!id])
      }
      momrat <- mom1^2/mom2
      
      if (length(theta0)== 0) {
        thetaneg <- rep(NA, mk-1)
        err <- rep(Inf, mk-1)
        mu1t <- rep(Inf, mk-1)
        for (j in (1:lg)){
          tj <- thetagrid[j]
          temp <- cumsum(h(tj, th[1:mk-1]))/L[1:mk-1]
          mu1 <- th[2:mk]^(-tj)*temp + h(tj, 1/th[2:mk])
          tj2 <- tj*2
          temp <- cumsum(h(tj2, th[1:mk]))/L[1:mk]
          temp <- th[1:mk-1]^(-tj2)*temp + h(tj2, 1/th[2:mk])
          mu2 <- (temp-mu1)*(2/tj)
          murat <- mu1^2/mu2
          err1 <- abs(log(momrat/murat))
          id <- err1< err
          thetaneg[id] <- tj
          err[id] <- err1[id]
          mu1t[id] <- mu1[id] 
        }
 
        # now for thetapos abd combine
        mu1_0 <- cumsum(log(th[1:(mk-1)]))/L(1:(mk-1))-log(th[2:mk])
        thetapos <- mom1/mu1_0
        theta <- thetaneg + thetapos
        theta <- theta[k]
          
        # Asymptotic standard deviation of theta
        thetaStd= sqrt(r11/k)*th[k]
        # sigma2m <- sigma2 #for use in estimation of stand. dev. of quantile
      } else {
        theta <- rep(theta0[1], nl)
        if (length(theta0Std)> 0){
          thetaStd <- rep(theta0Std[1], nl)
        } else {
          thetaStd <- rep(0, nl)
        }
      }
      
      # Scale estimator
      if (length(logdisp0)== 0) {
        mm <- mom1[l]/mu1t[l]
        g <- mm*X[l]
        
        logdisp <- log(mm)  # Log of dispersion coefficient
        logdispStd <- sqrt(r11/l)
      } else {
        g <- X[l]*exp(logdisp0[1])
        logdisp <- rep(logdisp0[1], nl)
        if (length(logdisp0Std)> 0){
          logdispStd <- rep(logdisp0Std[1], nl)
        } else {
          logdispStd <- rep(0 ,nl)
        }        
      }
      
      # Standard deviation of X0[l] as estimator of location q(th[l])
      X0lStd <- g*sqrt(r11/l)/th(l)
      
      # Quantile estimation
      lp= length(p)
      if (lp> 0) {
        q <- matrix(NA, nl, lp)
        qStd <- q
        for (i in 1:lp) {
          lambda <- -log(p[i])/th[l] # factor of -logs of probabilities
          
          # Quantiles
          q[, i] <- X[l]+g*h(theta, lambda)
          
          # Asymptotic standard deviations of quantiles
          ha <- h(theta,lambda)
          dha <- (1/theta)*(lambda^theta*log(lambda)-ha)
          id <- abs(theta)< .Machine$double.eps
          if (any(id)) {dha[id] <- 0.5*(log(lambda))^2}
          # the following asymptotic expression is pretty accurate
          # (the last term can normally be ignored but with given, precise,
          # theta and logdisp estimates, it may not be negligible)
          var <- g^2*(ha^2*logdispStd^2+dha^2*thetaStd^2) + X0lStd^2
          qStd[, i]= sqrt(var)
        }
      }
    }
    
    #
    # compute Boucheron-Thomas (cf. Drees-Kaufmann) type simple choice of 
    # threshold l
    # 
    #    estimatesBT <- NULL  
    # estimatesP0 <- NULL  
    # P <- NULL
    # bias <- NULL  
    # 
    # if (length(l0)== 0) {
    #   if (length(theta0)== 0) { 
    #     vtest <- theta
    #     vtestStd <- thetaStd    
    #   } else if (length(logdisp0)== 0) {
    #     vtest <- logdisp
    #     vtestStd <- logdispStd  
    #   } else {
    #     vtest <- X0[l]
    #     vtestStd <- X0lStd              
    #   }
    #   
    # i <- selectThresholdBT(vtest, vtestStd, l, 50)
    # estimatesBT <- list("k"= k[i], "l"= l[i], "y"= th[l[i]], 
    #                     "N"= N, "sigma"= sqrt(sigma2),"r11"= r11,
    #                     "tailindex"= theta[i], "scale"= g[i], 
    #                     "location"= X0[l[i]], "locationStd"= X0lStd[i],
    #                     "logdisp"= logdisp[i], "logdispStd"= logdispStd[i],
    #                     "p"= p, "quantile"= q[i, ], 
    #                     "tailindexStd"= thetaStd[i], "quantileStd"= qStd[i, ], 
    #                     "df"= "GW", 
    #                     "estimator"= "iteratedHill", "XId"= XId)
    #
    # compute refined estimate including estimate of bias
    #
    # rP0 <- selectThresholdP0(vtest, vtestStd, l)
    # i <- rP0$i
    # Pfluctuation <- rP0$P
    # bias <- rP0$bias
    # estimatesP0 <- list("k"= k[i], "l"= l[i], "y"= th[l[i]], 
    #                     "N"= N, "sigma"= sqrt(sigma2), "r11"= r11,
    #                     "tailindex"= theta[i], "scale"= g[i], 
    #                     "logdisp"= logdisp[i], "logdispStd"= logdispStd[i],
    #                     "location"= X0[l[i]], "locationStd"= X0lStd[i],
    #                     "p"= p, "quantile"= q[i, ], 
    #                     "tailindexStd"= thetaStd[i], "quantileStd"= qStd[i, ], 
    #                     "df"= "GW", 
    #                     "estimator"= "iteratedHill", "XId"= XId)
    # }
    # Below is all outdated stuff, don't use it
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
    
    estimates <- list("k"= k, "l"= l, "y"= th[l], 
                      "N"= N, "sigma"= sqrt(sigma2), "r11"= r11,
                      "tailindex"= theta, "tailindexStd"= thetaStd, 
                      "scale"= g, "logdisp"= logdisp, "logdispStd"= logdispStd,
                      "location"= X0[l], "locationStd"= X0lStd,
                      "p"= p, "quantile"= q, "quantileStd"= qStd, 
                      "orderstats"= X0, "df"= "GW", 
                      "estimator"= "iteratedHill", "XId"= XId)
                      # "estimatesBT"= estimatesBT,  # Boucheron-Thomas estimate
                      # "Pfluctuation"= Pfluctuation,# fluctuation size p-value
                      # "bias"= bias,                # order of magnitude of bias
                      # "estimatesP0"= estimatesP0)  # single-threshold-estimates (a la B-T)
                      # "estimatesP"= estimatesP,    # multi-threshold mean (jump points)
                      # "estimatesPP"= estimatesPP)  # multi-threshold median (jump points)
                      
  } # if (n > 0)
  
  return(estimates)
}

