#' @name  FitGW_iHilli
#' 
#' @title FitGW_iHilli
#' 
#' @description Fit a Generalised Weibull (GW) upper tail to the sample X and estimate quantiles: variant
#' of FitGW_iHill in which an implicit joint solution of tail index and scale is compued 
#' 
#' @param X data sample (double(n))
#' @param p probabilities of exceedance of the quantiles to be estimated (double(np))  
#' @param N (optional) (effective) sample size, in case X is not complete but contains only (peak) values above some threshold (integer(1))
#' @param r11 (optional) factor to increase estimator variance by, to account for serial dependence (default: 1) (double(1))
#' @param fixedpar (optional): fixed model parameters not to be estimated, and their standard errors (list; see Details)
#' @param l0 (optional) value of l (no. of order stats used) in case it is imposed (integer(0))
#' @param sigma (optional) determines the ratio of k to l (double(1))
#' @param XId (optional) data identifier to store with output for traceability (character)
#' 
#' @usage Value <- FitGW_iHilli(X, p, N= 0, r11= 1, fixedpar= NULL, l0= NULL, sigma= 1, XId= '')
#' 
#' @return A list, with members: 
#'   \item{l}{no. of order statistics used for scale and quantile estimation}    
#'   \item{k}{no. of order statistics used for tail index estimation} 
#'   \item{sigma}{algorithm parameter (see ref. eq. (30))}
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
FitGW_iHilli <- function(X, p, N, r11, fixedpar, l0, sigma, XId) {

  # Handle arguments
  if (missing(p)) {p <- NULL}
  if (missing(N)) {N <- 0} 
  if (missing(r11)) {r11 <- 1}
  if (missing(fixedpar)) {fixedpar <- NULL}
  if (missing(l0)) {l0 <- NULL}
  if (missing(sigma)) {sigma <- 1}
  if (missing(XId)) {XId <- ''}
  
  # fixed parameters 
  sigma2 <- sigma^2
  
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
    k <- l  # start Newton iteration for k 
    if (sigma2< Inf) {
      for (jj in 1:10) {
        k <- k-(k-l*th[k]^2)/(1+2*l*th[k]/k)
      }
    }
    k <- pmin(pmax(1, k), n);
    
    # Order statistics, decreasing
    X0 <- -sort(-X)
    
    # Adjust l and k based on requirements 
    # ind <- which(k> 2 & X0[k]> -Inf & l> 0 & k> l*2 & k< n)
    ind <- which(k> 2 & X0[k]> -Inf & l> 0  & k<= n & k>= l)
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
    
    if (nl> 0) {
      # Hill estimator
      hill0 <- cumsum(X0[1:(mk-1)])/L[1:(mk-1)]-X0[2:mk]
      id <- (hill0<= 0)
      if (any(id)) {hill0[id] <- min(hill0[!id])}
      
      if (length(theta0)== 0) {
        # The mean excess and the mean excesss of its logarithm
        hill1 <- cumsum(log(hill0[1:(mk-2)]))/L[1:(mk-2)]-log(hill0[2:(mk-1)])
        
        # u is defined as in proof of Theorem 2 of ref. 
        u <- cumsum(log(th[2:(n-1)]))/(1:(n-2))-log(th[3:n])
        
        # Simple estimator of GW index (nondimensional curvature)
        theta <- 1+hill1[k-2]/u[k-2]
        
        # Asymptotic standard deviation of theta
        thetaStd= th[k]*sqrt(r11/k)
        sigma2m <- sigma2 #for use in estimation of stand. dev. of quantile
      } else {
        theta <- rep(theta0[1], nl)
        if (length(theta0Std)> 0){
          thetaStd <- rep(theta0Std[1], nl)
        } else {
          thetaStd <- rep(0, nl)
        }
      }
      
      # Refinement of GW index estimator
      thetagrid <- seq(-1.99, 1.99, .02)   # wide range
      lg <- length(thetagrid)
      
      err <- rep(Inf, mk-1)    # error tracking (lowest value sofar)
      g <- rep(NA, nl)         # scale estimates
      thetaref <- rep(NA, nl)
      for (i in 1:lg) {
        ti <- thetagrid[i]
        temp <- cumsum(h(ti, th[1:(mk-1)]))/L[2:mk]
        w <- th[2:mk]^(-ti)*temp + h(ti, 1/th[2:mk])
        w1 <- cumsum(log(w[1:(mk-2)]))/L[1:(mk-2)]-log(w[2:(mk-1)])
        
        err1 <- abs(ti + 1 - theta + w1[k-2]/u[k-2])
        id <- (err1< err)
        if (!any(is.na(id))) {
          thetaref[id] <- tj
          err[id] <- err1[id]
          # g <- hill0[l-1]/normg
          g[id] <- hill0[l[id]-1]/w[l[id]-1]
        }
      }
      
      # Scale estimator (continued)
      if (length(logdisp0)== 0) {
        logdisp <- log(g/pmax(X0[l], .0001))  # Log of dispersion coefficient
        logdispStd <- sqrt(r11/l)
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
      X0lStd <- hill0[l-1]*sqrt(r11/l)
      
      # Quantile estimation
      lp= length(p)
      if (lp> 0) {
        q <- matrix(NA, nl, lp)
        qStd <- q
        for (i in 1:lp) {
          lambda <- -log(p[i])/th[l] # factor of -logs of probabilities
          
          # Quantiles
          q[, i] <- X0[l]+g*h(thetaref, lambda)
          
          # Asymptotic standard deviations of quantiles
          ha <- h(thetaref, lambda)
          dha <- (1/thetaref)*(lambda^thetaref*log(lambda)-ha)
          id <- abs(thetaref)< .Machine$double.eps
          if (any(id)) {dha[id] <- 0.5*(log(lambda))^2}
          # the following asymptotic expression is pretty accurate
          # (the last term can normally be ignored but with given, precise,
          # theta and logdisp estimates, it may not be negligible)
          var <- g^2*(ha^2*logdispStd^2+dha^2*thetaStd^2) + X0lStd^2
          qStd[, i]= sqrt(var)
        }
      }
    }
    
    estimates <- list("k"= k, "l"= l, "y"= th[l], 
                      "N"= N, "sigma"= sqrt(sigma2), "r11"= r11,
                      "tailindex"= theta, "tailindexref"= thetaref,
                      "tailindexStd"= thetaStd, 
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

