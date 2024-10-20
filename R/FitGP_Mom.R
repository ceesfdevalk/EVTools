#' @title FitGP_Mom
#' 
#' @description Fit GP tail using moment method (see reference)
#' 
#' @param X data sample (double(n))
#' @param p probabilities of exceedance of the quantiles to be estimated (double(np))  
#' @param N (optional) (effective) sample size, in case X is not complete but contains only (peak) values above some threshold (integer(1))
#' @param r11 (optional) factor to increase estimator variance by, to account for serial dependence (default: 1) (double(1) or list, see Details)
#' @param fixedpar (optional): fixed model parameters not to be estimated, and their standard errors (double(1) or list, see Details)
#' @param l0 (optional) value of l (no. of order stats used) in case it is imposed (integer(0))
#' @param sigma (optional) determines the ratio of k to l, see below (double(1))
#' @param metadata (optional) data identifier to store with output for traceability (character)
#' 
#' @usage Value <- FitGP_Mom(X, p, N= 0, r11= 1, fixedpar= NULL, l0= NULL, sigma= Inf,metadata= NULL)
#' 
#' @return A list, with members: 
#'   \item{l}{no. of order statistics used for scale and quantile estimation}    
#'   \item{k}{no. of order statistics used for tail index estimation} 
#'   \item{tailindex}{estimates or imposed value of GP tail index} 
#'   \item{tailindexStd}{standard deviations of tail index estimates}
#'   \item{logdisp}{estimates or imposed value of log of dispersion coeff.}  
#'   \item{logdispStd}{standard deviations of log of dispersion coeff. estimates}
#'   \item{scale}{estimates of GP scale parameter}
#'   \item{locationStd}{standard deviation of order statistic}
#'   \item{lambda}{ratio of logarithms of probabilities of exceedance of quantile and threshold}  
#'   \item{p}{probabilities of exceedance of quantiles to be estimated} 
#'   \item{quantile}{quantile estimates}
#'   \item{quantileStd}{standard deviations of quantile estimates}
#'   \item{orderstats}{data X sorted (decreasing)}
#'   \item{df}{= "GP": fitted distribution function tail (Generalised Pareto)}
#'   \item{estimator}{= "maximum likelihood": see "method" below}
#' 
#' @details
#'  
#'  Pre-determined model parameters are to be supplied in the list fixedpar (see above):
#'  \itemize{
#'   \item{$gamma0: (optional) value of tailindex in case it is imposed (double(1))}
#'   \item{$gamma0Std: (optional) its standard deviation (double(1))}
#'   \item{$logdisp0: (optional) value of log of dispersion coeff. in case it is imposed (dispersion coeff. is the raio of scale par. to location par.) (double(1))}
#'   \item{$logdisp0Std: (optional) its standard deviation (double(1))}        
#'   }
#'   
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
#'   The serial dependence coefficient r11 can be a positive number, or a list 
#'   produced by R11.R. 
#'   
#'   In case a quantile is to be estimated for a \emph{frequency}, say f, and 
#'   \enumerate{
#'   \item{if X contains all values (possibly above some threshold), then with
#'   EI an estimate of the Extremal Index from EI.R, set
#'   p = f*d/EI and N = T/d, with T the length of the observation period and d the time step. 
#'         Note that f and d are defined with reference to the same unit of time!! In this case,
#'         r11 needs to be estimated.
#'        }
#'   \item{if X contains only the n (approximately Poisson) peak values above some threshold 
#'         (in a PoT analysis),  it is recommended to set r11= 1 and take p = f*d/EI and 
#'         N = T/d*EI; in this case (for GP), EI can be any value; e.g. take p= fT/n and N= n.
#'        } 
#' }   
#'           
#' @references
#' De Haan, L. and A. Ferreira (2006), Extreme Value Theory - An Introduction. Springer. 
#
#' @author Cees de Valk \email{ceesfdevalk@gmail.com}
#' 
#' @export
FitGP_Mom <- function(X, p, N, r11, fixedpar, l0, sigma, metadata) {
  
  # Handle arguments
  if (missing(p)) {p <- NULL}
  if (missing(N)) {N <- 0} 
  if (missing(r11)) {r11 <- 1}
  if (missing(fixedpar)) {fixedpar <- NULL}
  if (missing(l0)) {l0 <- NULL}
  if (missing(sigma)) {sigma <- Inf}
  if (missing(metadata)) {metadata= NULL}
  
  # fixed parameter 
  gamma0 <- fixedpar$gamma0
  gamma0Std <- fixedpar$gamma0Std
  logdisp0 <- fixedpar$logdisp0
  logdisp0Std <- fixedpar$logdisp0Std
  
  X <- c(X)
  n <- length(X)  
  estimates <- NULL
  
  if (n > 0) {
    
    # Correct N if invalid or smaller than n 
    if (N < n) {N <- n}  # this is actually very wrong if POT is done
    
    # Harmonic numbers (H) and expectations of exponential order statistics (see (15)) 
    # L <- 1:n 
    # H <- c(0, cumsum(1/L))
    # th <- sum(1/(1:N))-H
    # th <- th[1:(n*5)]
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
    ind <- which(k> 2 & X0[k]> -Inf & l> 0 & k>= l & k< n)
    if (length(ind)> 0) {
      k <- k[ind]
      l <- l[ind]
    } else {
      k <- NULL
      l <- NULL
    }
    p0 <- l/N
    y <- log(N/l)
    
    #
    # to do: possibly further reduce size of arrays l and k here
    # ....
    
    nl <- length(k)
    mk <- max(k)
    ml <- max(l)
    
    gamma <- rep(NA, nl)
    g <- rep(1, nl)    # needed as starting value
    
    #
    # main loop
    #
    if (nl> 0) {
      # Empirical moments 
      mom1 <- cumsum(X00[1:(mk-1)])/(1:(mk-1))-X00[2:mk]
      mom2 <- cumsum(X00[1:(mk-1)]^2)/(1:(mk-1)) - X00[2:mk]^2 - 2*mom1*X00[2:mk]
      id <- pmin(mom1, mom2)<= 0
      if (any(id)) {
        mom1[id] <- min(mom1[!id])
        mom2[id] <- min(mom2[!id])
      }
      momrat <- mom1^2/mom2
      gammas <- mom1 + 1 - 0.5/(1-momrat)
      gamma <- gammas[k-1]
      
      # standard deviations etc.
      # NOTE: if k= l, the asymptotic covariance matrix of the parameters 
      # possibly has nonzero diagonal elements 
      
      # Asymptotic standard deviation of gamma ?????? ()
      if (is.list(r11)) {
        r11value <- approx(r11$p, r11$r, k/N, rule= 2)$y 
      } else {
        r11value <- rep(r11, length(k))
      }
      gammaStd <- sqrt(r11value/k*(1+gamma^2))  
      id <- gamma<0
      gd <- gamma[id]
      gammaStd[id] <- (1-gd)*sqrt(r11value[id]/k[id]*(1-2*gd)*(1-gd+6*gd^2)/(1-3*gd)/(1-4*gd))
      gammaStd <- rev(cummax(rev(gammaStd)))  # to avoid unrealistic small values
      
      if (length(gamma0)> 0){
        gamma <- rep(gamma0[1], nl)
        if (length(gamma0Std)> 0){
          gammaStd <- rep(gamma0Std[1], nl)
        } else {
          gammaStd <- rep(0, nl)
        }      
      }
      
      # scale
      g <- X0[l]*mom1[l-1]*(1-pmin(0, gamma))
       
       
      if (is.list(r11)) {
        r11value <- approx(r11$p, r11$r, l/N, rule= 2)$y 
      } else {
        r11value <- r11
      }
      
      # Scale estimator
      if (length(logdisp0)== 0) {
        logdisp <- log(g/X0[l])  # Log of dispersion coefficient
        logdispStd <- sqrt(r11value/l)*sqrt(2+gamma^2) # from de Haan & Ferreira
        id <- gamma<0
        gd <- gamma[id]
        logdispStd[id] <- sqrt(r11value[id]/l[id]*(2-16*gd+51*gd^2-69*gd^3+50*gd^4-24*gd^5)/
                                 (1-2*gd)/(1-3*gd)/(1-4*gd))
        
        # logdispStd <- sqrt(r11value/l)
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
      X0lStd <- g*sqrt(r11value/l)
      
      # Quantile estimation
      lp= length(p)
      if (lp> 0) {
        q <- matrix(NA, nl, lp)
        qStd <- q
        for (i in 1:lp) {
          lambda <- p0/p[i] # factor of -logs of probabilities
          
          # Quantiles
          q[, i] <- X0[l]+g*h(gamma, lambda)
          
          # Asymptotic standard deviations of quantiles  All Questionable!!!
          ha <- h(gamma, lambda)
          dha <- (1/gamma)*(lambda^gamma*log(lambda)-ha)
          id <- abs(gamma)< 1.e-10
          if (any(id)) {dha[id] <- 0.5*(log(lambda))^2}
          # the following asymptotic expression is pretty accurate
          # if l<< k (sigma< Inf)
          # (the last term can normally be ignored but with given, precise,
          # gamma and logdisp estimates, it may not be negligible)
          var <- g^2*(ha^2*logdispStd^2 + dha^2*gammaStd^2) + X0lStd^2
          ind <- (k== l)
          if (sum(ind)>0) {
            # dependence term, specifically for MLE if k=l (from de Haan & Ferreira)
            depterm <- 2*g^2*ha*dha*r11value/l*(gamma-1)  
            id <- gamma<0
            gd <- gamma[id]
            depterm[id] <- 2*g[id]^2*ha[id]*dha[id]*r11value[id]/l[id]*
              (1-gd)^2*(-1+4*gd-12*gd^2)/((1-3*gd)*(1-4*gd))
            var[ind] <- var[ind] + depterm[ind]
          }
          qStd[, i]= sqrt(var)
          qStd[, i] <- rev(cummax(rev(qStd[, i])))  # to avoid unrealistic small values     
        }
      }
      
      estimates <- list("k"= k, "l"= l, "y"= y, 
                        "N"= N, "sigma"= sigma, "r11"= r11,
                        "tailindex"= gamma, "tailindexStd"= gammaStd, 
                        "scale"= g, "logdisp"= logdisp, "logdispStd"= logdispStd,
                        "location"= X0[l], "locationStd"= X0lStd,
                        "p"= p, "quantile"= q, "quantileStd"= qStd, 
                        "orderstats"= X0, "df"= "GP", 
                        "method"= "FitGP_Mom", "metadata"= metadata)
      # "estimatesBT"= estimatesBT,  # Boucheron-Thomas estimate
      # "Pfluctuation"= Pfluctuation,# fluctuation size p-value
      # "bias"= bias,                # order of magnitude of bias
      # "estimatesP0"= estimatesP0)  # single-threshold-estimates (a la B-T)
      # "estimatesP"= estimatesP,    # multi-threshold mean (jump points)
      # "estimatesPP"= estimatesPP)  # multi-threshold median (jump points)
      
    }  
  }
  return(estimates)
}

