#' @name  FitTail_AllDataCov
#' 
#' @title FitTail_AllDataCov
#' 
#' @description Fit a tail to the values of a time-series and estimate upper tail quantiles
#' 
#' @param X data sample (double(n) or (double(n, 2)))
#' @param freq frequencies of exceedance of the quantiles to be estimated (double(nf))  
#' @param df distribution function to be fitted to the tail: "GP", "GW", "logGW", "Wbl", or "Exp"
#' @param method (optional) name of R-script to estimate the tail (character)
#' @param options (optional) parameters controlling the estimation (list; see Details)
#' @param metadata (optional) information about the variable and the time-series (list; see Details)
#' 
#' @usage Value <- FitTail_AllDataCov(X, freq= NULL, df= "GW", method= "FitGW_iHilli", options= NULL, metadata= NULL)
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
#'   \item{p}{probabilities (frwctions of time) of exceedance of quantiles to be estimated} 
#'   \item{p0}{probability of exceedance of the lower bound of the variable} 
#'   \item{freq}{frequencies of exceedance of quantiles to be estimated}    
#'   \item{quantile}{quantile estimates}
#'   \item{quantileStd}{standard deviations of quantile estimates}
#'   \item{tailindexraw}{raw estimates of GW tail index over all possible thresholds (method: FitGW_iHill.R)} 
#'   \item{tailindexrawStd}{standard deviation of tailindexraw}
#'   \item{kraw}{no. of order statistics used for estimation of tailindexraw} 
#'   \item{orderstats}{data X sorted (decreasing)}
#'   \item{df}{see above}
#'   \item{method}{see above}
#' In addition, several plots are produced:
#'   \item{}
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
#'   
#' } 
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
#'   \item{$dither: width of uniform distribution of noise to add to data (double(1))}
#'   \item{$pthreshold: fraction of time that value exceeds threshold (double(1))}
#'   \item{$maxpthreshold: upper bound on pthreshold (in case pthreshold is estimated)}
#'   \item{$minpthreshold: lower bound on pthreshold (in case pthreshold is estimated)}
#'   \item{$indexselect: if TRUE, threshold is selected based on tail index estimates (logical, default= FALSE)} 
#'   \item{$kmin: no. of order statistics skipped in determining threshold (integer(1)), default= 20)} 
#'   \item{$bootstrap: list. If exists/nonempty, precision is assessed by a moving block bootstrap. May 
#'          contain $nsamples (no. of bootstrap samples) and $blocktime (block length in terms of time)
#'   \item{$fixedpar: fixed model parameters not to be estimated, and their standard errors (list; see below)}
#'   \item{$plotparams: plotparameters (list) with members: $makeplot (default= TRUE), $pconf (coverage probability 
#'          of confidence interval), $xlim (plot limits for quantile estimates), $freqlim (plot limits for 
#'          frequencies), $plim (plot limits for fractions of time)}
#'  }                                
#
#' @author Cees de Valk \email{cees.de.valk@knmi.nl}
#' 
#' @export
FitTail_AllDataCov <- function(X, freq, df, method, options, metadata) {
  
  # Default parameters
  
  fixedpar <- sigma <- pthreshold <-  maxpthreshold <- minpthreshold <- indexselect <- NULL
  
  # Handle missing arguments    
  
  if (missing(X)) {stop("Data X must be specified.")}
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
  
  timestep <- metadata$timestep
  if (length(timestep)< 1) {
    timestep <- 1      # This makes the probability equal to the frequency
  }
  
  # Estimator options
  
  nbin <- options$nbin
  dither <- options$dither
  pthreshold <- options$pthreshold
  maxpthreshold <- options$maxpthreshold
  if (length(maxpthreshold)< 1) {maxpthreshold <- 0.5} 
  minpthreshold <- options$minpthreshold
  if (length(minpthreshold)< 1) {minpthreshold <- 0}
  if (minpthreshold> maxpthreshold) {
    stop("options$minpthreshold larger or equal to options$maxpthreshold")
  } else if (minpthreshold== maxpthreshold) {
    minpthreshold <- pthreshold <- maxpthreshold
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
  kmin <- options$kmin
  if (length(kmin)< 1) {kmin <- 20}  # default value
  
  #
  # Prepare bootstrap
  #
  nb <- options$bootstrap$nsamples
  bt <- options$bootstrap$blocktime
  blocks <- options$bootstrap$blocks
  nblocks <- 0
  if (length(nb)< 1) {nb <- 0}
  
  if (length(blocks)> 0) {
    blocks <- data.matrix(blocks)
    nb <- dim(blocks)[2]
    nblocks <- dim(blocks)[1]
  } else if (nb> 0) {
    nb <- max(50, nb) # no. of bootstrap samples not too small!
  }
}

if (!(length(bt)> 0)) {
  bt <- 1  # block length is 1 time unit, normally a year
  warning("Bootstrap block length set to 1 time unit.")
}
lb <- ceiling(bt/timestep) # block length measured in time steps (rounded upward)

if ((nblocks== 0) & (nb> 0)) {
  nblocks <- ceiling(length(X)/lb)
  blocks <- sample(length(X)-lb+1, size= nblocks*nb, replace= TRUE)
  dim(blocks) <- c(nblocks, nb)
}

# Determine quantization and dither data if needed

if (length(dither)< 1) {
  sX <- -sort(-X)
  dX <- -diff(sX)
  dither <- min(dX[dX> 0]) 
  if (max(dX%%dither)/dither> 0.1) {
    dither <- NULL
  }
}
if (length(dither)> 0) {
  ind <- X> 0.5*dither
  X <- X + (runif(sum(ind))-0.5)*dither
}
options$dither <- 0  # when called recursively, nothing is done 

#
# Covariate: determine categories
#
X <- data.matrix(X)
N <- dim(X)[1] 
if (N< 20) {stop("Time series length must be at least 20.")}

cat <- rep(1, N)
if (dim(X)[2]> 1) {
  cat <- X[, 2]
  X <- X[, 1]
  if (length(nbin)> 0) {
    if (nbin>= 2) { 
      binw <- (max(cat)-min(cat))/nbin
      dsa <- diff(sort(cast))
      delta <- min(dsa[dsa> 0])
      binw <- ceil(binw/delta)*delta
      cat <- (((cat+binw/2) %/% binw) %% nbin)*binw
    }
  }
  cats <- sort(unique(cat))
  lcats <- length(cats)
}

#
# Loop over covariate bins (by recursion)
#
estimates <- vector(mode= "list", length= lcats)
for (i in 1:lcats) {
  
  ind <- cat== cats[i]
  
  EIes <- EI(X[ind], makeplot= FALSE)
  r11es <- r11(X[ind], makeplot= FALSE)
  EIvalue <- max(EIes$EIFS[1:3]) 
  
  Xmin <- min(X[ind])
  N <- sum(X[ind]> Xmin)+1 
  pcat <- N/length(X)   # fraction of time that X is in cats[i] and above its
  # minimal value for cats[i]
  
  if (N> 20) {
    
    #  p is fraction of "the time that X is above its minimum and cat in cats[i]"
    p <- freq*timestep/EIvalue/pcat
    
    sX <- -sort(-X[ind])
    n <- round(min(N*maxpthreshold, 5.e4))
    if (length(pthreshold)<1) {
      l0 <- seq(kmin, n, kmin)             # only the values needed for threshold optimisation  
    } else {
      l0 <- round(N*pthreshold)
    }
    
    es <- get(tailfit)(X=sX[1:n], method, p=p, N=N, r11=r11es, fixedpar= fixedpar, 
                       l0= l0, sigma= sigma, metadata= metadata)
    
    es$cat <- cats[i]
    es$pcat <- pcat  # fraction of time that X is above its minimum and cat in cats[i]
    es$p <- p*pcat # real fraction of time
    es$freq <- freq  # frequency
    es$EIvalue <- EIvalue
    
    if (nb> 0) {   # do bootstrap
      
      be <- vector("list", length = nb)  
      for (j in 1:nb) {
        is <- rep(blocks[, j], each= lb) + rep(0:(lb-1), times= nblocks)
        Xb <- X[is[1:length(X)]]
        catb <- cat[is[1:length(X)]]
        
        ind <- catb== cats[i]
        sXb <- -sort(-Xb[ind])
        be[[j]] <- get(tailfit)(X=sXb[1:n], method, p=p, N=N, r11= 1, fixedpar= fixedpar, 
                                l0= l0, sigma= sigma, metadata= metadata)
      }
      
      # store bootstrap ensemble of estimates 
      tailindex <- purrr:map(be, "tailindex")
      logdisp <- purrr:map(be, "logdisp")
      scale <- purrr:map(be, "scale")
      location <- purrr:map(be, "logdisp")
      quantile <- purrr:map(be, "quantile")
      k <- purrr:map(be, "k")
      l <- purrr:map(be, "l")
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
      
      # Compute quantiles for selected threshold on a refined frequency grid 
      
      ls <- es$l[iselect]
      lf <- log10(freq)     # Extend frequency array for plotting etc. 
      mlf <- log10(es$l[iselect]/es$N*EIvalue*p0/timestep)
      freqs <- 10^(min(lf) + (mlf-min(lf))*seq(0, 1, 0.01))
      
      ps <- freqs*timestep/EIvalue/pcat
      ps <- unique(sort(c(ps, p)))
      ps <- pmax(0, pmin(1, ps))
      freqs <- ps*EIvalue*pcat/timestep
      
      esel <- get(tailfit)(X=sX[1:n], method, p=ps, N=N, r11=r11es, fixedpar= fixedpar, 
                           l0= ls, sigma= sigma, metadata= metadata)
      
      esel$pcat <- pcat  # fraction of time that X is above its minimum
      esel$p <- ps*pcat # fraction of time
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
    if (is.null(plotparams$makeplot)) {plotparams$makeplot <- TRUE}
    if (plotparams$makeplot) {
      # Plot of tail fit
      genname <- paste(estimates$df, "-", metadata$varname, "-", metadata$caseId, sep= "")
      
      fname <- paste("Tail-", genname, ".png", sep= "")
      png(filename= fname,units="in", width=7.5*fac, height=7.5*fac, res=72)
      tailplot(plotparams, es$selected)
      dev.off()
      
      # Plot of tail index estimates vs. l
      fname <- paste("Tailindex-", genname, ".png", sep= "")
      png(filename= fname,units="in", width=5*fac, height=5*fac, res=72)
      tailindexplot(es= es)  
      dev.off()
      
      # Plot of quantile estimates vs. l for the lowest freqency
      fname <- paste("Quantile-", genname, ".png", sep= "")
      png(filename= fname,units="in", width=5*fac, height=5*fac, res=72)
      tailquantileplot(plotparams, es) 
      dev.off()
      
      # Plot P-value
      fname <- paste("ThresholdP-", genname, ".png", sep= "")
      png(filename= fname,units="in", width=5*fac, height=5*fac, res=72)
      thresholdplot(plotparams, es)
      dev.off()
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
