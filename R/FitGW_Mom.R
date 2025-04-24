#' @name  FitGW_Mom
#' 
#' @title FitGW_Mom
#' 
#' @description Fit a Generalised Weibull (GW) upper tail to the sample X and estimate quantiles, using
#' the moment estimator of Albert et al. (2018).
#' 
#' @param X data sample (double(n))
#' @param p probabilities of exceedance of the quantiles to be estimated (double(np))  
#' @param N (optional) (effective) sample size, in case X is not complete but contains only (peak) values above some threshold (integer(1))
#' @param r11 (optional) factor to increase estimator variance by, to account for serial dependence (default: 1) (double(1) or list, see Details)
#' @param fixedpar (optional): fixed model parameters not to be estimated, and their standard errors (list; see Details)
#' @param l0 (optional) value of l (no. of order stats used) in case it is imposed (integer(0))
#' @param sigma (optional) fixed algorithm parameter (see de Valk & Cai (2018) eq. (30)) (double(1)) 
#' @param metadata (optional) information about the variable and, if applicable, the time-series (list; see Details)
#' 
#' @usage Value <- FitGW_Mom(X, p, N= 0, r11= 1, fixedpar= NULL, l0= NULL, sigma= Inf, metadata= NULL)
#' 
#' @return A list, with members: 
#'   \item{l}{no. of order statistics used for scale and quantile estimation}    
#'   \item{k}{no. of order statistics used for tail index estimation} 
#'   \item{sigma}{fixed algorithm parameter (see ref. eq. (30))}
#'   \item{tailindex}{estimates or imposed value of GW tail index} 
#'   \item{tailindexStd}{standard deviations of tail index estimates}
#'   \item{logdisp}{estimates or imposed value of log of dispersion coeff.}  
#'   \item{logdispStd}{standard deviations of log of dispersion coeff. estimates}
#'   \item{scale}{estimates of GW scale parameter}
#'   \item{locationStd}{standard deviation of order statistic}
#'   \item{lambda}{ratio of logarithms of probabilities of exceedance of quantile and threshold}  
#'   \item{p}{probabilities of exceedance of quantiles to be estimated} 
#'   \item{quantile}{quantile estimates}
#'   \item{quantileStd}{standard deviations of quantile estimates}
#'   \item{tailindexraw}{estimates of GW tail index over all possible thresholds} 
#'   \item{tailindexrawStd}{standard deviation of tailindexraw}
#'   \item{kraw}{no. of order statistics used for estimation of tailindexraw} 
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
#'   The serial dependence coefficient r11 can be a positive number, or a list 
#'   produced by R11.R. 
#'   
#'   In case a quantile is to be estimated for a \emph{frequency} f, and 
#'   \enumerate{
#'   \item{if X contains all values (possibly above some threshold), then with
#'   EI an estimate of the Extremal Index from EI.R, set
#'   p = f*d/EI and N = T/d, with T the length of the observation period and d the time step. 
#'         Note that f and d are defined with reference to the same unit of time!! In this case,
#'         r11 needs to be estimated.
#'        }
#'   \item{if X contains only the n (approximately Poisson) peak values above some threshold 
#'         (in a PoT analysis),  it is recommended to set r11= 1 and take p = f*d and 
#'         N = T/d*EI. EI need to be estimated (see above). In this case, EI can also be 
#'         estimated also as EI= n*d/Tt= n/nt with Tt the time spent above the threshold and 
#'         nt the number of time-series values above the threshold. 
#'        } 
#' }   
#'  metadata may contain the following fields (in addition to your own meta data):
#'  \itemize{
#'   \item{$varname: variable name}
#'   \item{$varunit: physical unit of variable}
#'   \item{$timeunit: time unit (e.g. year)}
#'   \item{$timestep: time step in units of timeunit}
#'   \item{$timelength: length of time covered by time-series, in units of timeunit} 
#'   \item{$EI: extremal index (see above)}
#'   \item{$nexcess (for PoT only): no. of data values (as opposed to peak values) exceeding the threshold}
#'  }   
#'           
#' @references
#' Albert, C., Dutfoy,A., Gardes, L., Girard, S. (2018), An extreme quantile estimator 
#' for the log-generalized Weibull-tail model. See: https://hal.inria.fr/hal-01783929v1.
#' De Valk, C. and Cai, J.J. (2018), A high quantile estimator based on 
#' the log-generalized Weibull tail limit. Econometrics and Statistics 6, 107-128, see
#' \url{https://doi.org/10.1016/j.ecosta.2017.03.001}
#
#' @author Cees de Valk \email{ceesfdevalk@gmail.com}
#' 
#' @export
FitGW_Mom <- function(X, p, N, r11, fixedpar, l0, sigma, metadata) {
library(gsl)
  
  # Handle arguments
  if (missing(p)) {p <- NULL}
  if (missing(N)) {N <- 0} 
  if (missing(r11)) {r11 <- 1}
  if (missing(fixedpar)) {fixedpar <- NULL}
  if (missing(l0)) {l0 <- NULL}
  if (missing(sigma)) {sigma <- Inf}
  if (missing(metadata)) {metadata <-NULL}
  
  # fixed parameter 
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
    # H <- c(0, cumsum(1/L))
    # th <- sum(1/(1:N))-H
    # th <- th[1:n]
    # H <- H[2:(n+1)]
    th <- log(N/(1:n))
    
    #
    # k is found by iteration from l (see (30))
    #
    if (length(l0)== 0) {
      l <- 10:(n-1)
    } else {
      l <- l0[l0< n]
    }
    nl <- length(l)
    k <- l  # Newton iteration for k 
    if (sigma< Inf){
      for (jj in 1:10) {
        thk <- log(N/k)
        k <- k-(k-l*(thk/sigma+1)^2)/(1+2*l*(thk/sigma+1)/(k*sigma))
      }
    }
    k <- round(k)
    k <- pmin(pmax(1, k), n);
    k <- pmax(k, l)      # experimental
    
    # Order statistics of logarithms, decreasing
    X0 <- -sort(-X)
    X00 <- log(pmax(X0,0))
    
    # Adjust l and k based on requirements 
    ind <- which(k> 2 & X00[k]> -Inf & l> 0 & k>= l & k< n)
    if (length(ind)> 0) {
      k <- k[ind]
      l <- l[ind]
    } else {
      k <- NULL
      l <- NULL
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
    thetaraw <- NULL
    thetarawStd <- NULL
    kraw <- NULL    
    
    thetagrid <- seq(-0.99, 1.99, .02)  # not zero, for convenience; range can be extended to above .99
                                       # but not below -.99 (at least not without modifying the formula
                                       # for mu1b and possibly for mu1 as well, see below)
    lg <- length(thetagrid)
    
    if (ml> 0) {
      # Empirical moments 
      mom1 <- cumsum(X00[1:(mk-1)])/L[1:(mk-1)]-X00[2:mk]
      mom2 <- cumsum(X00[1:(mk-1)]^2)/L[1:(mk-1)] - X00[2:mk]^2 - 2*mom1*X00[2:mk]
      id <- pmin(mom1, mom2)<= 0
      if (any(id)) {
        mom1[id] <- min(mom1[!id])
        mom2[id] <- min(mom2[!id])
      }
      momrat <- mom1^2/mom2
      
      minmurat <- Inf
      maxmurat <- -Inf
      if (length(theta0)== 0) {
        thetaneg <- rep(NA, mk-1)
        err <- rep(Inf, mk-1)
        mu1t <- rep(Inf, mk-1)
        for (j in (1:lg)){
          tj <- thetagrid[j]
          
          # Discrete approximation (in same spirit as expressions in de Valk & Cai)
          # temp <- cumsum(h(tj, th[1:(mk-1)]))/L[1:(mk-1)]
          # hi <- h(tj, 1/th[2:mk])
          # mu1 <- th[2:mk]^(-tj)*temp + hi
          # temp1 <- cumsum(h(tj, th[1:(mk-1)])^2)/L[1:(mk-1)]
          # mu2 <- th[2:mk]^(-tj*2)*temp1 - hi^2 + 2*hi*mu1
          # mu1_0 <- cumsum(log(th[1:(mk-1)]))/L[1:(mk-1)]-log(th[2:mk])
          
          # Analytical expressions (implementing the integrals in Albert et al.)
          mu1 <- (exp(th[2:mk])*th[2:mk]^(-tj)*gamma_inc(tj+1, th[2:mk]) - 1)/tj
          mu1b <- (exp(th[2:mk])*th[2:mk]^(-tj*2)*gamma_inc(tj*2+1, th[2:mk]) - 1)/(tj*2)
          mu2 <- (mu1b-mu1)*(2/tj)
            
          murat <-  mu1^2/mu2
          
          minmurat <- pmin(minmurat, murat)
          maxmurat <- pmax(maxmurat, murat)
          
          # muratNum <- rep(NA, mk-1)
          # for (jk in 1:mk-1) {muratNum[jk] <- PsiNum(tj, th[jk+1])}
            
          err1 <- abs(log(momrat/murat)) # is this how the inversion is done?
          # cat(length(err))
          # cat(length(err1))
          id <- (err1< err)
          if (any(is.na(id))) 
            {cat("id contains NA")}
          else if (any(id)) {
            thetaneg[id] <- tj
            err[id] <- err1[id]
            mu1t[id] <- mu1[id] 
          }
        }
        
        # now for thetapos abd combine
        # mu1_0 <- exp(th[2:mk])*gamma_inc(0, th[2:mk])
        mu1_0 <- expint_E1(th[2:mk])*exp(th[2:mk])
        thetapos <- mom1/mu1_0
        thetaraw <- thetaneg + thetapos
        kraw <- (1:(mk-1))+1
        # thetaraw[mk] <- thetaraw[mk-1]  # last one may be NA
        theta <- thetaraw[k-1]
  
        # Asymptotic standard deviation of theta
        if (is.list(r11)) {
          r11value <- approx(r11$p, r11$r, kraw/N, rule= 2)$y 
        } else {
          r11value <- r11
        }
        # thetaStd= th[k]*sqrt(r11value/k) replaced 
        thetarawStd= th[kraw]*sqrt(r11value/kraw)
        thetaStd= thetarawStd[k-2]
        thetaStd <- rev(cummax(rev(thetaStd)))  # to avoid unrealistic small values
        

      } else {
        theta <- rep(theta0[1], nl)
        if (length(theta0Std)> 0){
          thetaStd <- rep(theta0Std[1], nl)
        } else {
          thetaStd <- rep(0, nl)
        }
      }
      
      if (is.list(r11)) {
        r11value <- approx(r11$p, r11$r, l/N, rule= 2)$y 
      } else {
        r11value <- r11
      }
      
      # Scale estimator
      if (length(logdisp0)== 0) {
        mm <- mom1[l]/mu1t[l]
        g <- mm*X0[l]
        
        logdisp <- log(mm)  # Log of dispersion coefficient
        logdispStd <- sqrt(r11value/l)
      } else {
        g <- X0[l]*exp(logdisp0[1])
        logdisp <- rep(logdisp0[1], nl)
        if (length(logdisp0Std)> 0){
          logdispStd <- rep(logdisp0Std[1], nl)
        } else {
          logdispStd <- rep(0 ,nl)
        }        
      }
      
      # Standard deviation of X0[l] as estimator of location q(th[l])
      X0lStd <- g*sqrt(r11value/l)/th[l]
      
      # Quantile estimation
      lp= length(p)
      if (lp> 0) {
        q <- matrix(NA, nl, lp)
        qStd <- q
        for (i in 1:lp) {
          lambda <- -log(p[i])/th[l] # factor of -logs of probabilities
          
          # Quantiles
          q[, i] <- X0[l]+g*h(theta, lambda)
          
          # Asymptotic standard deviations of quantiles
          ha <- h(theta,lambda)
          dha <- (1/theta)*(lambda^theta*log(lambda)-ha)
          id <- abs(theta)< 1.e-10
          if (any(id)) {dha[id] <- 0.5*(log(lambda))^2}
          # the following asymptotic expression is pretty accurate
          # (the last term can normally be ignored but with given, precise,
          # theta and logdisp estimates, it may not be negligible)
          var <- g^2*(ha^2*logdispStd^2+dha^2*thetaStd^2) + X0lStd^2
          qStd[, i]= sqrt(var)
          qStd[, i] <- rev(cummax(rev(qStd[, i])))  # to avoid unrealistic small values     
        }
      }
    }
    
    estimates <- list("k"= k, "l"= l, "y"= th[l], 
                      "N"= N, "sigma"= sigma, "r11"= r11,
                      "tailindex"= theta, "tailindexStd"= thetaStd, 
                      "scale"= g, "logdisp"= logdisp, "logdispStd"= logdispStd,
                      "location"= X0[l], "locationStd"= X0lStd,
                      "p"= p, "quantile"= q, "quantileStd"= qStd, 
                      "tailindexraw"= thetaraw, "tailindexrawStd"= thetarawStd, "kraw"= kraw,
                      "orderstats"= X0, "df"= "GW", 
                      "method"= "FitGW_Mom", "metadata"= metadata)
                      # "estimatesBT"= estimatesBT,  # Boucheron-Thomas estimate
                      # "Pfluctuation"= Pfluctuation,# fluctuation size p-value
                      # "bias"= bias,                # order of magnitude of bias
                      # "estimatesP0"= estimatesP0)  # single-threshold-estimates (a la B-T)
                      # "estimatesP"= estimatesP,    # multi-threshold mean (jump points)
                      # "estimatesPP"= estimatesPP)  # multi-threshold median (jump points)
                      
  } # if (n > 0)
  
  return(estimates)
}

