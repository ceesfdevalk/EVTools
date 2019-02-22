#' @title FitGP_MLE2
#' 
#' @description As FitGP_MLE: test version with distinct thresholds for scale and shape 
#' 
#' @param X data sample (double(n))
#' @param p probabilities of exceedance of the quantiles to be estimated (double(np))  
#' @param N (optional) (effective) sample size, in case X is not complete but contains only (peak) values above some threshold (integer(1))
#' @param r11 (optional) factor to increase estimator variance by, to account for serial dependence (default: 1) (double(1) or list, see Details)
#' @param fixedpar (optional): fixed model parameters not to be estimated, and their standard errors (double(1) or list, see Details)
#' @param l0 (optional) value of l (no. of order stats used) in case it is imposed (integer(0))
#' @param sigma (optional) determines the ratio of k to l, see below (double(1))
#' @param XId (optional) data identifier to store with output for traceability (character)
#' 
#' @usage Value <- FitGP_MLE2(X, p, N= 0, r11= 1, fixedpar= NULL, l0= NULL, sigma= Inf, XId= '')
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
#'             }
#'   \item{if X contains only the (approximately Poisson) peak values above some threshold 
#'         (i.e., you want to do a PoT analysis), then set p = f*d/EI, N = (T/d)*EI, r11= 1 
#'         (note that d/EI is the mean duration of an "event" and EI/d is the 
#'         mean number of "events" per unit of time estimated by EI.R).
#'              }
#' }   
#'           
#' @references
#' De Haan, L. and A. Ferreira (2006), Extreme Value Theory - An Introduction. Springer. 
#
#' @author Cees de Valk \email{ceesfdevalk@gmail.com}
#' 
#' @export
FitGP_MLE2 <- function(X, p, N, r11, fixedpar, l0, sigma, XId) {
  
  # Handle arguments
  if (missing(p)) {p <- NULL}
  if (missing(N)) {N <- 0} 
  if (missing(r11)) {r11 <- 1}
  if (missing(fixedpar)) {fixedpar <- NULL}
  if (missing(l0)) {l0 <- NULL}
  if (missing(sigma)) {sigma <- 1}
  if (missing(XId)) {XId <- ''}
  
  # fixed parameter 
  sigma2 <- sigma^2
  
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
      l <- 10:n
    } else {
      l <- l0
    }
    nl <- length(l)
    k <- l  # Newton iteration for k if sigma2> 0
    if (sigma2< Inf){
      for (jj in 1:10) {
        thk <- log(N/k)
        k <- k-(k-l*thk^2/sigma2)/(1+2*l*thk/k/sigma2)
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
    # function to be minimized
    #
    negllGP <- function(par) {
      #
      # negative log-likelihood of conditional GP distribution given the exceedance
      # of min(x)
      #
      k <- length(xglobal)
      k1 <- k-1
      logg <- par[1]
      g <- exp(logg)
      if (length(par)< 2) {
        gamma <- gammaglobal
      } else {
        gamma <- par[2]  
      }
      z <- (xglobal[1:k1]-xglobal[k])/g
      if (abs(gamma)< 1.e-10) {
        f <- k1*logg + sum(z)
      } else {
        f <- k1*logg + (1/gamma+1)*sum(log(pmax(0, 1+gamma*z)))
      }
      if (min(1+z*gamma)<= 0) {f <- Inf}
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
        
        if (length(gamma0)== 0) {
          xglobal <- X0[1:kj]
          Nglobal <- N
          par0 <- c(0, 0.1)
          optimout <- try(optim(par0, negllGP, method= "BFGS"), silent=TRUE)
          if (class(optimout)== 'try-error') {
            optimout <- try(optim(par0, negllGP, method= "Nelder-Mead"), silent=TRUE)
          }
          if (class(optimout)!= 'try-error') {   
            par1 <- optimout$par
            gamma[j] <- par1[2]
            g[j] <- exp(par1[1])
          }
        } else {
          gamma[j] <- gamma0[1]
        }
        
        if ((lj< kj)| (length(gamma0)> 0)) {           # then estimate scale at a different threshold
          par2 <- log(g[j])
          gammaglobal <- gamma[j]
          xglobal <- X0[1:lj]
          optimout <- optim(par2, negllGP, method= "Brent", lower= par2-10, upper= par2+10)
          par3 <- optimout$par
          g[j] <- exp(par3)
        }
        setTxtProgressBar(pb, j)
      }   # for j in ....
      
      # standard deviations etc.
      # NOTE: if k= l, the asymptotic covariance matrix of the parameters 
      # possibly has nonzero diagonal elements 
      
      # Asymptotic standard deviation of gamma ?????? ()
      if (is.list(r11)) {
        r11value <- approx(r11$k, r11$r, k, rule= 2)$y 
      } else {
        r11value <- r11
      }
      gammaStd= sqrt(r11value/k)*(1+gamma)  
      
      if (length(gamma0)> 0){
        if (length(gamma0Std)> 0){
          gammaStd <- rep(gamma0Std[1], nl)
        } else {
          gammaStd <- rep(0, nl)
        }      
      }
      
      if (is.list(r11)) {
        r11value <- approx(r11$k, r11$r, l, rule= 2)$y 
      } else {
        r11value <- r11
      }
      
      # Scale estimator
      if (length(logdisp0)== 0) {
        logdisp <- log(g/X0[l])  # Log of dispersion coefficient
        logdispStd <- sqrt(r11value/l)*sqrt(1+(1+gamma)^2) # from de Haan & Ferreira
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
      X0lStd <- g*sqrt(r11value/l)/th[l]
      
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
            depterm <- -2*g^2*ha*dha*(1+gamma)*r11value/l  
            var[ind] <- var[ind] + depterm[ind]
          }
          qStd[, i]= sqrt(var)
        }
      }
      
      estimates <- list("k"= k, "l"= l, "y"= y, 
                        "N"= N, "sigma"= sqrt(sigma2), "r11"= r11,
                        "tailindex"= gamma, "tailindexStd"= gammaStd, 
                        "scale"= g, "logdisp"= logdisp, "logdispStd"= logdispStd,
                        "location"= X0[l], "locationStd"= X0lStd,
                        "p"= p, "quantile"= q, "quantileStd"= qStd, 
                        "orderstats"= X0, "df"= "GP", 
                        "estimator"= "Maximum likelihood", "XId"= XId)
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

