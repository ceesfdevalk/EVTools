#' @name  FitTail_AllData
#' 
#' @title FitTail_AllData
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
#' @usage Value <- FitTail_AllData(X, freq= NULL, df= "GW", method= "FitGW_iHilli", options= NULL, metadata= NULL)
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
#'   \item{$sigma: determines the ratio of k to l ( (no. of order stats used for estimation of tail index and quantile) (double(1)}
#'   \item{$bootstrap: list. If exists/nonempty, precision is assessed by a moving block bootstrap. May contain $nsamples (no. of bootstrap samples) and $blocktime (block length in terms of time)
#'   \item{$fixedpar: fixed model parameters not to be estimated, and their standard errors (list; see below)}
#'   \item{$plotparams: plotparameters (list) with members: $makeplot (default= TRUE), $pconf (coverage probability of confidence interval), $xlim (plot limits for quantile estimates), $freqlim (plot limits for frequencies), $plim (plot limits for fractions of time)}
#'  }                                
#
#' @author Cees de Valk \email{cees.de.valk@knmi.nl}
#' 
#' @export
FitTail_AllData <- function(X, freq, df, method, options, metadata) {
  
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
    metadata$caseId <- Sys.time()
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
  if (length(timestep)< 1) {timestep= 1} # This makes the probability equal to the frequency
  
  # Estimator options
  delta <- options$dither
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
  
  sigma <- options$sigma
  if (length(sigma)< 1) {sigma <- Inf} # to keep behaviour simple to non-expert
  fixedpar <- options$fixedpar
  indexselect <- options$indexselect
  if (length(indexselect)< 1) {indexselect <- TRUE} # index is better!
  kmin <- options$kmin
  if (length(kmin)< 1) {kmin <- 20}
  
  nb <- options$bootstrap$nsamples
  bt <- options$bootstrap$blocktime
  if (length(nb)> 0) {
    if (nb< 50) {
      nb <- max(50, nb) # no. of bootstrap samples not too small!
      warning("No. of bootstrap samples adjusted upward to 50.")
    }
    if (!length(bt)> 0) {
      bt <- 1  # block length is 1 time unit, normally a year
      warning("Bootstrap block length set to 1 time unit.")
    }
    lb <- ceiling(bt/timestep) # block length measured in time steps (rounded upward)
  } else {
    nb <- 0
  }

  # Sample size and correction for positive probability of X equal to its lower bound,
  # to prevent fitting of distribution containing an atom at its lowest value, 
  # like with rainfall

  X <- data.matrix(X)
  N <- dim(X)[1] 
  if (N< 20) {stop("Time series length must be at least 20.")}
  
  cat <- rep(1, N)
  if (dim(X)[2]> 1) {
    cat <- X[, 2]
    X <- X[, 1]
    if (any(is.na(cat)) & !is.na(cat[1])) {
      nbin <- cat[1]
      if (nbin> 1) { 
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
  
  Xmin <- min(X)
  N <- sum(X> Xmin)+1 
  if (N< 20) {stop("Time series must have at least 19 values above its minimum.")}
  p0 <- N/length(X) # fraction of time that X is above its minimum
  
  # Determine quantization and dither data if needed
  
  if (length(delta)< 1) {
    sX <- -sort(-X)
    dX <- -diff(sX)
    delta <- min(dX[dX> 0]) 
    if (max(dX%%delta)/delta> 0.1) {
      delta <- NULL
    }
  }
  if (length(delta)> 0) {
    X <- X + (runif(length(X))-0.5)*delta
  }
  X <- pmax(X, Xmin) # to prevent a change of range due to dithering
  
  # Estimate extremal index EI and dependence coefficient r11
  
  EIes <- EI(X, makeplot= FALSE)
  r11es <- r11(X, makeplot= FALSE)
  EIvalue <- max(EIes$EIFS[1:3]) 
  
  # Convert frequency to fraction of time p
  
  p <- freq*timestep/EIvalue/p0 #  fraction of "the time that X is above its minimum"
  p <- p[p>0 & p< 1]
  if (length(p)< 1) {p <- NULL}
  
  # Tail estimation
  
  l0 <- N*pthreshold
  if (length(l0)<1) {
    l0 <- NULL
  } else {
    l0 <- min(l0, round(N*min(pthreshold, p0)))
  }
  
  
  sX <- -sort(-X)
  n <- min(N, 5.e4)
  estimates <- get(tailfit)(X=sX[1:n], method, p=p, N=N, r11=r11es, fixedpar= fixedpar, 
                       l0= l0, sigma= sigma, metadata= metadata)
  estimates$p0 <- p0  # fraction of time that X is above its minimum
  estimates$p <- p*p0 # fraction of time
  estimates$freq <- freq  # frequency
  estimates$EIvalue <- EIvalue
  iselect <- which(estimates$l== round(N*pthreshold))
  
  if (nb> 0) {            
    
    # Moving block bootstrap to overrule asymptotic approximations of precision
    
    nblocks <- ceiling(length(X)/lb)
    istart <- sample(length(X)-lb+1, size= nblocks*nb, replace= TRUE)
    dim(istart) <- c(nblocks, nb)
    be <- vector("list", length = nb)  
    for (j in 1:nb) {
      is <- rep(istart[, j], each= lb) + rep(0:(lb-1), times= nblocks)
      Xb <- X[is[1:length(X)]]
      sXb <- -sort(-Xb)
      be[[j]] <- get(tailfit)(X=sXb[1:n], method, p=p, N=N, r11= 1, fixedpar= fixedpar, 
                              l0= l0, sigma= sigma, metadata= metadata)
    }
    
    # Process estimates from bootstrap samples to obtain precision estimates 
    tt <- unlist(purrr::map(be, "tailindex"))
    dim(tt) <- c(length(be[[1]]$tailindex), length(be))
    tStd <- apply(tt, 1, sd)
    estimates$tailindexStd <- rev(cummax(rev(tStd)))
    
    qq <- unlist(purrr::map(be, "quantile"))
    dim(qq) <- c(length(be[[1]]$quantile), length(be))   
    qStd <- apply(qq, 1, sd)
    dim(qStd) <- dim(be[[1]]$quantile)
    qStd <- apply(qStd, 2, f <- function(x) {rev(cummax(rev(x)))})
    estimates$quantileStd <- qStd
    
    ll <- unlist(purrr::map(be, "logdisp"))
    dim(ll) <- c(length(be[[1]]$logdisp), length(be))   
    lStd <- apply(ll, 1, sd)
    estimates$logdispStd <- rev(cummax(rev(lStd)))
  }
  
  # Threshold choice
  
  if (length(pthreshold) <1) {
    if (!indexselect) {
      Pthresh <- selectThresholdP1(estimates$quantile[, 1], estimates$quantileStd[, 1], 
                                   estimates$l, 0.9, kmin= kmin) 
      li <- Pthresh$k[Pthresh$i]
      iselect <- which(estimates$l== li)
    } else {
      Pthresh <- selectThresholdP1(estimates$tailindex, estimates$tailindexStd, estimates$k, 
                                   0.5, kmin= kmin)
      ki <- Pthresh$k[Pthresh$i]
      iselect <- which(estimates$k== ki)
    }  
    estimates$threshold <- Pthresh
    # Bound iselect from above by preset limit pmax on probability of exceedance
    lmax <- estimates$N*maxpthreshold
    lmin <- estimates$N*minpthreshold
    iselect <- min(iselect, max(which(estimates$l< lmax)))
    iselect <- max(iselect, min(which(estimates$l> lmin)))
    if (estimates$l[iselect]< lmin) {
      warning("options$minpthreshold overruled: sample is large, not all data can be processed.")
    }
  }
  
  # Compute quantiles for selected threshold on refined frequency grid 
  
  ls <- estimates$l[iselect]
  lf <- log10(freq)     # Extend frequency array for plotting etc. 
  mlf <- log10(estimates$l[iselect]/estimates$N*EIvalue*p0/timestep)
  freqs <- 10^(min(lf) + (mlf-min(lf))*seq(0, 1, 0.01))
  ps <- freqs*timestep/EIvalue/p0
  ps <- unique(sort(c(ps, p)))
  ps <- pmax(0, pmin(1, ps))
  freqs <- ps*EIvalue*p0/timestep
  
  es <- get(tailfit)(X=sX[1:n], method, p=ps, N=N, r11=r11es, fixedpar= fixedpar, 
                     l0= ls, sigma= sigma, metadata= metadata)
  
  es$p0 <- p0  # fraction of time that X is above its minimum
  es$p <- ps*p0 # fraction of time
  es$freq <- freqs  # frequency
  es$EIvalue <- EIvalue
  
  estimates$selected <- es
  
  # Plotting
  
  fac <- 1.2
  plotparams <- options$plotparams
  if (is.null(plotparams)) {plotparams <- NULL}
  if (is.null(plotparams$makeplot)) {plotparams$makeplot <- TRUE}
  if (plotparams$makeplot) {
    # Plot of tail fit
    genname <- paste(estimates$df, "-", metadata$varname, "-", metadata$caseId, sep= "")
    
    fname <- paste("Tail-", genname, ".png", sep= "")
    png(filename= fname,units="in", width=7.5*fac, height=7.5*fac, res=72)
    tailplot(plotparams, es)
    dev.off()
    
    # Plot of tail index estimates vs. l
    fname <- paste("Tailindex-", genname, ".png", sep= "")
    png(filename= fname,units="in", width=5*fac, height=5*fac, res=72)
    tailindexplot(es= estimates)  
    dev.off()
    
    # Plot of quantile estimates vs. l for the lowest freqency
    fname <- paste("Quantile-", genname, ".png", sep= "")
    png(filename= fname,units="in", width=5*fac, height=5*fac, res=72)
    tailquantileplot(plotparams, estimates) 
    dev.off()
    
    # Plot P-value
    fname <- paste("ThresholdP-", genname, ".png", sep= "")
    png(filename= fname,units="in", width=5*fac, height=5*fac, res=72)
    thresholdplot(plotparams, estimates)
    dev.off()
  }
  
  return(estimates)
}
