#' @name  FitTail_PoT
#' 
#' @title FitTail_PoT
#' 
#' @description Fit a tail to a PoT (peaks-over-threshold) sample from a time series, and estimate upper tail quantiles
#' 
#' @param X data sample (double(n))
#' @param freq frequencies of exceedance of the quantiles to be estimated (double(nf))  
#' @param df distribution function to be fitted to the tail: "GP", "GW", "logGW", "Wbl", or "Exp"
#' @param method (optional) name of R-script to estimate the tail (character)
#' @param options (optional) parameters controlling the estimation (list; see Details)
#' @param metadata (optional) information about the variable and the time-series (list; see Details)
#' 
#' @usage Value <- FitTail_PoT(X, freq= NULL, df= "GW", method= "FitGW_iHilli", options= NULL, metadata= NULL)
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
#' 
#'  options may contain the following fields:
#'  \itemize{
#'   \item{$pthreshold: fraction of time that value exceeds threshold (double(1))}
#'   \item{$maxpthreshold: upper bound on pthreshold (in case pthreshold is estimated)}
#'   \item{$minpthreshold: lower bound on pthreshold (in case pthreshold is estimated)}
#'   \item{$indexselect: if TRUE, threshold is selected based on tail index estimates (logical, default= FALSE)} 
#'   \item{$kmin: no. of order statistics skipped in determining threshold (integer(1)), default= 20)} 
#'   \item{$sigma: determines the ratio of k to l ( (no. of order stats used for estimation of tail index and quantile) (double(1)}
#'   \item{$fixedpar: fixed model parameters not to be estimated, and their standard errors (list; see below)}
#'   \item{$plotparams: plotparameters (list) with members: $makeplot (default= TRUE), $pconf (coverage probability of confidence interval), $xlim (plot limits for quantile estimates), $freqlim (plot limits for frequencies), $plim (plot limits for fractions of time)}
#'  }                                
#
#' @author Cees de Valk \email{cees.de.valk@knmi.nl}
#' 
#' @export
FitTail_PoT <- function(X, freq, df, method, options, metadata) {
  
  # Default parameters
  fixedpar <- sigma <- pthreshold <-  maxpthreshold <- minpthreshold <- indexselect <- NULL
  
  # Handle missing arguments    
  if (missing(X)) {stop("Data X must be specified.")}
  if (missing(freq)) {
    freq <- NULL
  } else {
    freq <- sort(freq)
  }
  if (missing(df)) {df <- "GW"}
  if (missing(method)) {
    method <- "FitGW_iHilli"  # for now: for GP, the default will be changed
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
  
  # specific parameters
  timestep <- metadata$timestep
  if (length(timestep)< 1) {timestep= 1} # This makes the probability equal to the frequency
  
  # Estimator options
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
  sigma <- options$sigma
  if (length(sigma)< 1) {sigma <- Inf} # to keep behaviour simple to non-expert
  fixedpar <- options$fixedpar
  indexselect <- options$indexselect
  if (length(indexselect)< 1) {indexselect <- TRUE} # index is better!
  kmin <- options$kmin
  if (length(kmin)< 1) {kmin <- 20}
  
  # Sample size and correction for positive probability of X equal to its lower bound,
  # to prevent fitting of distribution containing an atom at its lowest value, 
  # like with rainfall
  N <- length(X) 
  if (N< 20) {stop("Time series length must be at least 20.")}
  Xmin <- min(X)
  N <- sum(X> Xmin)+1 
  if (N< 20) {stop("Time series must have at least 19 values above its minimum.")}
  p0 <- N/length(X) # fraction of time that X is above its minimum
  n <- N
  
  # Determine quantization and dither data if needed
  sX <- -sort(-X)
  dX <- -diff(X)
  deltaX <- min(dX[dX> 0])
  if (max(dX%%deltaX)/deltaX< 0.1) {
    X <- X + (runif(n)-0.5)*deltaX
  }
  X <- pmax(X, Xmin) # to prevent a change of range due to dithering
  
  # Estimate extremal index EI and dependence coefficient r11
  EIes <- EI(X, makeplot= FALSE)
  r11es <- r11(X, makeplot= FALSE)
  
  # Convert frequency to fraction of time p
  EIvalue <- max(EIes$EIFS[1:3]) 
  p <- freq*timestep/EIvalue/p0 #  fraction of "the time that X is above its minimum"
  p <- p[p>0 & p< 1]
  if (length(p)< 1) {p <- NULL}
  
  # Tail estimation
  sX <- -sort(-X)
  n <- min(n, 5.e5)
  l0 <- round(N*pthreshold)
  if (length(l0)<1) {l0 <- NULL}
  estimates <- get(tailfit)(X=sX[1:n], method, p=p, N=N, r11=r11es, fixedpar= fixedpar, 
                       l0= l0, sigma= sigma, metadata= metadata)
  estimates$p0 <- p0  # fraction of time that X is above its minimum
  estimates$p <- p*p0 # fraction of time
  estimates$freq <- freq  # frequency
  estimates$EIvalue <- EIvalue
  iselect <- which(estimates$l== round(N*pthreshold))
  
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
    iselect <- min(iselect, max(which(estimates$l< estimates$N*maxpthreshold)))
    iselect <- max(iselect, min(which(estimates$l> estimates$N*minpthreshold)))
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
  
  if (is.null(plotparams$plot)) {plotparams$makeplot <- TRUE}
  if (plotparams$plot) {
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
