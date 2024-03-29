#' @name  FitTailFromAllData
#' 
#' @title FitTailFromAllData
#' 
#' @description Fit and extrapolate the tail of the frequency distribution of a time-series. Optionally, the distribution of the values coinciding with a covariate in some bin is estimated, for each of a specified set of bins. 
#' 
#' @param X data sample (double(n), optionally with a second column of values of a covariate (double(n, 2)))
#' @param freq frequencies of exceedance of the quantiles to be estimated (double(nf))  
#' @param df distribution function to be fitted to the tail: "GP", "GW", "logGW", "Wbl", or "Exp"
#' @param method (optional) name of R-script to estimate the tail (character)
#' @param options (optional) parameters controlling the estimation (list; see Details)
#' @param metadata (optional) information about the variable and the time-series (list; see Details)
#'
#' @usage Value <- FitTailFromAllData(X, freq= NULL, df= "GW", method= "FitGW_iHilli", options= NULL, metadata= NULL)
#'
#' @return A list, with members: 
#'   \item{l}{no. of order statistics used for scale and quantile estimation}    
#'   \item{k}{no. of order statistics used for tail index estimation} 
#'   \item{sigma}{algorithm parameter (see ref. eq. (30))}
#'   \item{tailindex}{estimates or imposed value of GW tail index} 
#'   \item{tailindexStd}{standard deviations of tail index estimates}
#'   \item{logdisp}{estimates or imposed value of log of dispersion coeff.}  
#'   \item{logdispStd}{standard deviations of log of dispersion coeff. estimates}
#'   \item{scale}{estimates of GW scale parameter}
#'   \item{locationStd}{standard deviation of order statistic}
#'   \item{lambda}{ratio of logarithms of probabilities of exceedance of quantile and threshold}  
#'   \item{p}{probabilities (fractions of time) of exceedance of quantiles to be estimated} 
#'   \item{pbin}{probability (fraction of time) of exceedance of the lower bound of the variable} 
#'   \item{freq}{frequencies of exceedance of quantiles to be estimated}    
#'   \item{quantile}{quantile estimates}
#'   \item{quantileStd}{standard deviations of quantile estimates}
#'   \item{tailindexraw}{raw estimates of GW tail index over all possible thresholds (method: FitGW_iHill.R)} 
#'   \item{tailindexrawStd}{standard deviation of tailindexraw}
#'   \item{kraw}{no. of order statistics used for estimation of tailindexraw} 
#'   \item{orderstats}{data X sorted (decreasing)}
#'   \item{df}{see above}
#'   \item{method}{see above}
#'   \item{ }{In addition, several plots are produced (tailindex, quantile, threshold indicator, fitted tail for selected threshold)}
#'
#' @details
#'   
#'  metadata may contain the following fields (in addition to your own meta data):
#'  \itemize{
#'   \item{caseId: user-chosen identifier for later reference (default: current date/time)}
#'   \item{$varname: variable name}
#'   \item{$varunit: physical unit of variable}
#'   \item{$timeunit: time unit (e.g. year)}
#'   \item{$timestep: time step in units of timeunit}
#'   \item{$timelength: length of time covered by time-series, in units of timeunit} 
#'   \item{$nexcess (for PoT only): no. of data values (as opposed to peak values) exceeding the threshold}
#'  }  
#'  
#'  options may contain the following fields:
#'  \itemize{
#'   \item{$covariate: list, which may contain $lbin and $ubin (the lower resp. upper limits of the bins of covariate values on which to restrict the tail estimates)}
#'   \item{$dither: width of uniform distribution of noise to add to data (double(1))}
#'   \item{$pthreshold: fraction of good values exceeding threshold (double(1))}
#'   \item{$maxpthreshold: upper bound on pthreshold (in case pthreshold is estimated)}
#'   \item{$minpthreshold: lower bound on pthreshold (in case pthreshold is estimated)}
#'   \item{$indexselect: if TRUE, threshold is selected based on tail index estimates (logical, default= TRUE)} 
#'   \item{$kmin: no. of order statistics skipped in determining threshold (integer(1)), default= 20)} 
#'   \item{$sigma: determines the ratio of k to l ( (no. of order stats used for estimation of tail index and quantile) (double(1)}
#'   \item{$bootstrap: list. If exists/nonempty, precision is assessed by a moving block bootstrap. May contain $blocktime (block length in terms of time) and $nsamples (no. of bootstrap samples) or $blocks (random starting indices of the blocks, an array of size ($nblocks, $nsamples))}
#'   \item{$fixedpar: fixed model parameters not to be estimated, and their standard errors (list; see below)}
#'   \item{$plotparams: plotparameters (list) with members: $makeplot (default= TRUE), $pconf (coverage probability of confidence interval), $xlim (plot limits for quantile estimates), $freqlim (plot limits for frequencies), $plim (plot limits for fractions of time)}
#'   \item{$EI: value(s) of extremal index (representing serial dependence) to be used for all covariate bins (double(1)) or for each bin (double(length($covariate$lbin)))}
#'  }   
#'
#'  Pre-determined model parameters are to be supplied in options$fixedpar (see above):
#'  \itemize{
#'   \item{$theta0: (optional) value of tailindex in case it is imposed (double(1))}
#'   \item{$theta0Std: (optional) its standard deviation (double(1))}
#'   \item{$logdisp0: (optional) value of log of dispersion coeff. in case it is imposed (dispersion coeff. is the raio of scale par. to location par.) (double(1))}
#'   \item{$logdisp0Std: (optional) its standard deviation (double(1))}        
#'   }
#'   
#' @author Cees de Valk \email{cees.de.valk@knmi.nl}
#' 
#' @export
FitTailFromAllData <- function(X, freq, df, method, options, metadata) {
  library(purrr)
  
  # Default parameters
  
  fixedpar <- sigma <- pthreshold <-  maxpthreshold <- minpthreshold <- indexselect <- NULL
  
  # Handle missing arguments    
  
  if (missing(X)) {stop("Data X must be specified.")}
  X <- data.matrix(X)      # make sure X has a dimension
  N0 <- dim(X)[1] 
  id <- !is.na(rowSums(X))
  X <- data.matrix(X[id, ])
  
  N <- dim(X)[1] 
  if (N< 20) {
    stop("Time series length must be at least 20.")
  }
  
  if (missing(freq)) {
    freq <- NULL
  } else {
    freq <- sort(freq)
  }
  if (missing(df)) {
    df <- "GW"
    # Following is fair choice for common weather/ocean data  
    warning("Tail not specified; a GW tail is estimated.")
  }
  if (missing(method)) {
    method <- "FitGW_iHilli" 
    if (df== "GP") {
      method <- "FitGP_Mom"
    }
  }
  if (missing(options)) {options <-NULL}
  
  # Ensure that metadata exists and contains a unique identifier
  
  if (missing(metadata)) {
    metadata <- NULL
  }
  if (is.null(metadata$caseId)) {
    metadata$caseId <- Sys.time()   # always set at first call
  }
  # Specify estimator with corresponding options
  
  if (df== "Weibull") {tailfit <- "FitWbl"}
  if (df== "GW") {tailfit <- "FitGW"}  
  if (df== "logGW") {tailfit <- "FitlogGW"} 
  if (df== "GP") {
    tailfit <- "FitGP"
    # this method allows threshold optimization but is not the best ("FitGP_MLE" is better)
    if (method== "FitGW_iHilli") {method <- "FitGP_Mom"}  
  } 
  
  # Specific parameters
  timestep <- 1 
  timestep <- metadata$timestep
  timelength <- metadata$timelength
  if (length(timestep)< 1) {
    if (length(timelength)> 0) {
      timestep <- timelength/N0
    } else {
      timestep <- 1      # This makes the probability equal to the frequency
    }
  }
  
  # Estimator options
  
  lbin <- options$covariate$lbin
  ubin <- options$covariate$ubin
  dither <- options$dither
  pthreshold <- options$pthreshold
  maxpthreshold <- options$maxpthreshold
  minpthreshold <- options$minpthreshold
  if (length(pthreshold)> 0) {
    maxpthreshold <- minpthreshold <- pthreshold
  } else {
    if (length(maxpthreshold)< 1) {maxpthreshold <- 0.5} 
    if (length(minpthreshold)< 1) {minpthreshold <- 0}
    if (minpthreshold> maxpthreshold) {
      stop("options$minpthreshold larger or equal to options$maxpthreshold")
    } else if (minpthreshold== maxpthreshold) {
      pthreshold <- maxpthreshold
    }
  }
  if (grepl("ML", method)) {
    warning("Chosen method may take a long time")
  }
  
  # Only sigma= Inf to keep it simple
  # sigma <- options$sigma
  # if (length(sigma)< 1) {sigma <- Inf} # to keep behaviour simple to non-expert
  sigma <- Inf
  fixedpar <- options$fixedpar
  indexselect <- options$indexselect
  if (length(indexselect)< 1) {indexselect <- TRUE} # index is better!
  # kmin <- options$kmin
  # if (length(kmin)< 1) {
  #   kmin <- 20       # default value
  # }
  EI0 <- options$EI
  
  #
  # Prepare bootstrap
  #
  nb <- options$bootstrap$nsamples
  bt <- options$bootstrap$blocktime
  blocks <- options$bootstrap$blocks
  
  if (length(nb)< 1) {nb <- 0}
  
  nblocks <- 0
  if (length(blocks)> 0) {
    blocks <- data.matrix(blocks)
    nb <- dim(blocks)[2]
    nblocks <- dim(blocks)[1]
  } else if (nb> 0) {
    nb <- max(50, nb) # no. of bootstrap samples not too small!
  }
  
  if (!(length(bt)> 0)) {
    bt <- 1  # block length is 1 time unit, normally a year
    if (nb> 0) {
      warning("Bootstrap block length set to 1 time unit.")
    }
  }
  lb <- ceiling(bt/timestep) # block length measured in time steps (rounded upward)
  
  if ((nblocks== 0) & (nb> 0)) {
    nblocks <- ceiling(N/lb)
    blocks <- sample(N-lb+1, size= nblocks*nb, replace= TRUE)
    dim(blocks) <- c(nblocks, nb)
  }
  
  # Determine quantization and dither data if needed
  
  if (length(dither)< 1) {
    sX <- -sort(-X[, 1])
    dX <- -diff(sX)
    dither <- min(dX[dX> 0]) 
    if (max(dX%%dither)/dither> 0.1) {
      dither <- NULL
    }
  }
  if (length(dither)> 0) {
    ind <- X[, 1]> 0.5*dither
    X[ind, 1] <- X[ind, 1] + (runif(sum(ind))-0.5)*dither
  }
  
  #
  # Covariate: determine categories
  #
  cat <- data.matrix(rep(TRUE, N))
  lcats <- 1
  if (dim(X)[2]> 1) {
    cat <- assigncat(X[, 2], lbin, ubin)
    X <- X[, 1]
    lcats <- length(lbin)
  }
  
  if (length(EI0) %in% c(0, 1, lcats)) {
    if (length(EI0)== 1) {
      EI0 <- rep(EI0, lcats)
    }
  } else {  
    stop("Length of options$EI invalid.")
  }
  
  # certain inputs may be different for different bins
  corrlength <- function(x, ll) {
    if (ll>1 & length(x)== 1) {
      x <- rep(x[1], ll)  
    } else if (length(x)> 0 & length(x)!= ll) {
      stop("length of an array supplied in options is incompatible with no. of covariate bins.")
    }
    return(x)
  }
  pick <- function(x, i) {
    if (length(x)> 1) {
      x <- x[i]  
    }
    return(x)
  }
  pthreshold <- corrlength(pthreshold, lcats)
  fixedpar$theta0 <- corrlength(fixedpar$theta0, lcats)
  fixedpar$theta0Std <- corrlength(fixedpar$theta0Std, lcats)
  fixedpar$logdisp0 <- corrlength(fixedpar$logdisp0, lcats)
  fixedpar$logdisp0Std <- corrlength(fixedpar$logdisp0Std, lcats)
  
  #
  # Loop over covariate bins (by recursion)
  #
  estimates <- vector(mode= "list", length= lcats)
  for (i in 1:lcats) {
    
    if (lcats> 1) {    
      print(" ")
      print(c(lbin[i], ubin[i]))
      print(" ")
    }
    ind <- cat[, i]
    
    if (length(EI0)== 0) {
      EIes <- EI(X[ind], makeplot= FALSE)
      EIvalue <- max(EIes$EIFS[1:3]) 
    } else {
      EIvalue <- EI0[i]
    }
    r11es <- r11(X[ind], makeplot= FALSE)

    
    Xmin <- min(X[ind])
    N <- sum(X[ind]> Xmin)+1 
    pbin <- N/length(X)   # fraction of time that X is in bin and above its
    # minimal value
    
    if (N> 20) {
      
      #  p is fraction of "the time that X is above its minimum and X in bin"
      p <- freq*timestep/EIvalue/pbin
      
      kmin <- max(10, round(N*max(minpthreshold*0.5, 1.e-4)))
      
      sX <- -sort(-X[ind])
      n <- round(min(N*maxpthreshold*4, 5.e4*sqrt(kmin)))
      if (length(pthreshold[i])<1) {
        l0 <- seq(kmin, n, kmin)             # only the values needed for threshold optimisation  
      } else {
        l0 <- round(N*pick(pthreshold, i))
      } 
      
      fixedpar0 <- fixedpar
      fixedpar0$theta0 <- pick(fixedpar$theta0, i)
      fixedpar0$theta0Std <- pick(fixedpar$theta0Std, i)
      fixedpar0$logdisp0 <- pick(fixedpar$logdisp0, i)
      fixedpar0$logdispStd0 <- pick(fixedpar$logdispStd0, i)
      
      es <- get(tailfit)(X=sX[1:n], method, p=p, N=N, r11=r11es, fixedpar= fixedpar0, 
                         l0= l0, sigma= sigma, metadata= metadata)
      
      es$lbin <- lbin[i]
      es$ubin <- ubin[i]
      es$pbin <- pbin  # fraction of time that X is above its minimum and in bin
      es$p <- p*pbin # real fraction of time
      es$freq <- freq  # frequency
      es$EIvalue <- EIvalue
      
      if (nb> 0) {   # do bootstrap
        
        be <- vector("list", length = nb)  
        for (j in 1:nb) {
          is <- rep(blocks[, j], each= lb) + rep(0:(lb-1), times= nblocks)
          Xb <- X[is[1:length(X)]]
          catb <- cat[is[1:length(X)], ]
          
          ind <- catb[i, ]
          sXb <- -sort(-Xb[ind])
          be[[j]] <- get(tailfit)(X=sXb[1:n], method, p=p, N=N, r11= 1, fixedpar= fixedpar, 
                                  l0= l0, sigma= sigma, metadata= metadata)
        }
        
        # store bootstrap ensemble of estimates 
        tailindex <- purrr::map(be, "tailindex")
        logdisp <- purrr::map(be, "logdisp")
        scale <- purrr::map(be, "scale")
        location <- purrr::map(be, "location")
        quantile <- purrr::map(be, "quantile")
        k <- purrr::map(be, "k")
        l <- purrr::map(be, "l")
        es$bootstrap <- list(tailindex= tailindex, logdisp= logdisp, scale=scale,
                             location= location, quantile= quantile, k= k, l= l)
        
        # process to precision estimates 
        tt <- unlist(tailindex)
        dim(tt) <- c(length(be[[1]]$tailindex), length(be))
        tStd <- apply(tt, 1, sd)
        es$tailindexStd <- rev(cummax(rev(tStd)))
        
        qq <- unlist(quantile)
        dim(qq) <- c(length(be[[1]]$quantile), length(be))   
        qStd <- apply(qq, 1, sd)
        dim(qStd) <- dim(be[[1]]$quantile)
        qStd <- apply(qStd, 2, f <- function(x) {rev(cummax(rev(x)))})
        es$quantileStd <- qStd
        
        ll <- unlist(logdisp)
        dim(ll) <- c(length(be[[1]]$logdisp), length(be))   
        lStd <- apply(ll, 1, sd)
        es$logdispStd <- rev(cummax(rev(lStd)))
        
        ll <- unlist(location)
        dim(ll) <- c(length(be[[1]]$location), length(be))   
        lStd <- apply(ll, 1, sd)
        es$locationStd <- rev(cummax(rev(lStd)))
        
      } # end bootstrap
      
      # Threshold choice
      
      iselect <- 1
      
      if (length(pthreshold)< 1) {
        if (!indexselect) {
          Pthresh <- selectThresholdP1(es$quantile[, 1], es$quantileStd[, 1], 
                                       es$l, 0.9, kmin= kmin) 
          li <- Pthresh$k[Pthresh$i]
          iselect <- which(es$l== li)
        } else {
          Pthresh <- selectThresholdP1(es$tailindex, es$tailindexStd, es$k, 
                                       0.5, kmin= kmin)
          ki <- Pthresh$k[Pthresh$i]
          iselect <- which(es$k== ki)
        }  
        es$threshold <- Pthresh
        
        # Bound iselect from above by preset limit pmax on probability of exceedance
        lmax <- es$N*maxpthreshold
        lmin <- es$N*minpthreshold
        iselect <- min(iselect, max(which(es$l<= lmax)))
        iselect <- max(iselect, min(which(es$l>= lmin)))
        if (es$l[iselect]< lmin) {
          warning("options$minpthreshold overruled: sample is large, not all data can be processed.")
        }     
      }
      
      # Compute quantiles for selected threshold on a refined frequency grid 
      
      ls <- es$l[iselect]
      lf <- log10(freq)     # Extend frequency array for plotting etc. 
      mlf <- log10(es$l[iselect]/es$N*EIvalue*pbin/timestep)
      freqs <- 10^(min(lf) + (mlf-min(lf))*seq(0, 1, 0.01))
      
      ps <- freqs*timestep/EIvalue/pbin
      ps <- unique(sort(c(ps, p)))
      ps <- pmax(0, pmin(1, ps))
      freqs <- ps*EIvalue*pbin/timestep
      
      esel <- get(tailfit)(X=sX[1:n], method, p=ps, N=N, r11=r11es, fixedpar= fixedpar, 
                           l0= ls, sigma= sigma, metadata= metadata)
      
      esel$pbin <- pbin  # fraction of time that X is above its minimum
      esel$lbin <- lbin[i]
      esel$ubin <- ubin[i]
      
      esel$p <- ps*pbin # fraction of time
      esel$freq <- freqs  # frequency
      esel$EIvalue <- EIvalue
      
      esel$iselect <- iselect
      es$selected <- esel
      
      # Store result
      
      estimates[[i]] <- es
      
    } # if N< 20
    
  } # for (i in 1:lcats)
  
  # Plotting
  
  for (i in 1:lcats) {
    
    es <- estimates[[i]]
    
    fac <- 1.2
    plotparams <- options$plotparams
    if (is.null(plotparams)) {plotparams <- NULL}
    if (is.null(plotparams$makeplot)) {plotparams$makeplot <- FALSE}
    if (plotparams$makeplot) {
      # Plot of tail fit
      genname <- paste(es$df, "-", metadata$varname, "-", lbin[i], "-", ubin[i], "-", metadata$caseId, sep= "")
      
      fname <- paste("Tail-", genname, ".png", sep= "")
      png(filename= fname,units="in", width=7.5*fac, height=7.5*fac, res=72)
      tailplot(es$selected, params= plotparams)
      dev.off()
      
      if (length(es$l)> 1) {
        
        # Plot of tail index estimates vs. l
        fname <- paste("Tailindex-", genname, ".png", sep= "")
        png(filename= fname,units="in", width=5*fac, height=5*fac, res=72)
        tailindexplot(es= es)  
        dev.off()
        
        # Plot of quantile estimates vs. l for the lowest frequency
        fname <- paste("Quantile-", genname, ".png", sep= "")
        png(filename= fname,units="in", width=5*fac, height=5*fac, res=72)
        tailquantileplot(es, params= plotparams) 
        dev.off()
        
        # Plot P-value
        if (length(es$threshold)> 0) {
          fname <- paste("ThresholdP-", genname, ".png", sep= "")
          png(filename= fname,units="in", width=5*fac, height=5*fac, res=72)
          thresholdplot(plotparams, es)
          dev.off()
        }   
      }
    }
    
  } #  for (i in 1:lcats)
  
  # 
  # remove redundant list level
  # 
  if (length(estimates)== 1) {
    estimates <- estimates[[1]]
  }
  
  return(estimates)
}
