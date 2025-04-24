#' @name  FitGW_MLE
#' 
#' @title FitGW_MLE
#' 
#' @description Fit a Generalised Weibull (GW) upper tail to the sample X and estimate quantiles, using
#' the ML estimator for tail index and scale 
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
#' @usage Value <- FitGW_MLE(X, p, N= 0, r11= 1, fixedpar= NULL, l0= NULL, sigma= Inf, metadata= NULL)
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
#' De Valk, C. and Cai, J.J. (2018), A high quantile estimator based on 
#' the log-generalized Weibull tail limit. Econometrics and Statistics 6, 107-128, see
#' \url{https://doi.org/10.1016/j.ecosta.2017.03.001}
#' Ferreira & Dombry (??) ML for GEV as example
#
#' @author Cees de Valk \email{ceesfdevalk@gmail.com}
#' 
#' @export
FitGW_MLE <- function(X, p, N, r11, fixedpar, l0, sigma, metadata) {
  
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
    
    # Adjust l and k based on requirements 
    ind <- which(k> 2 & X0[k]> -Inf & l> 0 & k>= l & k<= (n-1))
    if (length(ind)> 0) {
      k <- k[ind]
      l <- l[ind]
    } else {
      k <- NULL
      l <- NULL
    }
    y <- log(N/l)
    
    #
    # to do: possibly further reduce size of arrays l and k here
    # ....
    
    nl <- length(k)
    mk <- max(k)
    ml <- max(l)
    
    theta <- rep(NA, nl)
    g <- rep(1, nl)    # needed as starting value
    
    #
    # function to be minimized
    #
    negllGW <- function(par) {
      #
      # negative log-likelihood of conditional GW distribution given the exceedance
      # of min(x), minimised for y= -log(p0) (see header)
      #
      k <- length(xglobal)
      k1 <- k-1
      logg <- par[1]
      g <- exp(logg)
      if (length(par)< 2) {
        theta <- thetaglobal
      } else {
        theta <- par[2]  
      }
      z <- (xglobal[1:k1]-xglobal[k])/g
      y <- log(Nglobal/k)
      if (abs(theta)< 1.e-10) {
        f <- -k1*log(y) + k1*logg - sum(z) + y*sum(exp(z))   
      } else {
        f <- -k1*log(y) + k1*logg - (1/theta-1)*sum(log(pmax(0, 1+theta*z))) +
          y*sum(pmax(0, 1+theta*z)^(1/theta))
      }
      if (min(1+z*theta)<= 0) {f <- Inf}
      f
    }
    
    #
    # main loop
    #
    if (nl> 0) {
      
      pb <- txtProgressBar(1, nl)
      for (j in (1:nl)) {
        lj <- l[j]
        kj <- k[j]
        
        if (length(theta0)== 0) {
          xglobal <- X0[1:kj]
          Nglobal <- N
          par0 <- c(0, 0.1)
          optimout <- try(optim(par0, negllGW, method= "BFGS"), silent=TRUE)
          if (class(optimout)== 'try-error') {
            optimout <- try(optim(par0, negllGW, method= "Nelder-Mead"), silent=TRUE)
          }
          if (class(optimout)!= 'try-error') {   
            par1 <- optimout$par
            theta[j] <- par1[2]
            g[j] <- exp(par1[1])
          }
        } else {
          theta[j] <- theta0[1]
        }
        
        if ((lj< kj)| (length(theta0)> 0)) {           # then estimate scale at a different threshold
          par2 <- log(g[j])
          thetaglobal <- theta[j]
          xglobal <- X0[1:lj]
          optimout <- optim(par2, negllGW, method= "Brent", lower= par2-10, upper= par2+10)
          par3 <- optimout$par
          g[j] <- exp(par3)
        }
        setTxtProgressBar(pb, j)
      }   # for j in ....
      
      # standard deviations etc.
      # NOTE: if k= l, the asymptotic covariance matrix of the parameters 
      # possibly has nonzero diagonal elements 
      
      # Asymptotic standard deviation of theta ?????? ()
      if (is.list(r11)) {
        r11value <- approx(r11$p, r11$r, k/N, rule= 2)$y 
      } else {
        r11value <- r11
      }
      thetaStd= th[k]*sqrt(r11value/k) # ?????? approximation: very similar to FitGW_iHilli.R
      thetaStd <- rev(cummax(rev(thetaStd)))  # to avoid unrealistic small values
      
      if (length(theta0)> 0){
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
        logdisp <- log(g/X0[l])  # Log of dispersion coefficient
        logdispStd <- sqrt(r11value/l)  # ??????
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
          lambda <- -log(p[i])/y # factor of -logs of probabilities
          
          # Quantiles
          q[, i] <- X0[l]+g*h(theta, lambda)
          
          # Asymptotic standard deviations of quantiles  All Questionable!!!
          ha <- h(theta, lambda)
          dha <- (1/theta)*(lambda^theta*log(lambda)-ha)
          id <- abs(theta)< 1.e-10
          if (any(id)) {dha[id] <- 0.5*(log(lambda))^2}
          # the following asymptotic expression is pretty accurate
          # (the last term can normally be ignored but with given, precise,
          # theta and logdisp estimates, it may not be negligible)
          var <- g^2*(ha^2*logdispStd^2 + dha^2*thetaStd^2) + X0lStd^2
          qStd[, i]= sqrt(var)
          qStd[, i] <- rev(cummax(rev(qStd[, i])))  # to avoid unrealistic small values     
        }
      }
      
      estimates <- list("k"= k, "l"= l, "y"= y, 
                        "N"= N, "sigma"= sigma, "r11"= r11,
                        "tailindex"= theta, "tailindexStd"= thetaStd, 
                        "scale"= g, "logdisp"= logdisp, "logdispStd"= logdispStd,
                        "location"= X0[l], "locationStd"= X0lStd,
                        "p"= p, "quantile"= q, "quantileStd"= qStd, 
                        "orderstats"= X0, "df"= "GW", 
                        "method"= "FitGW_MLE", "metadata"= metadata)
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

