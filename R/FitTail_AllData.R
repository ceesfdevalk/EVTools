#' @name  FitTail_AllData
#' 
#' @title FitTail_AllData
#' 
#' @description Fit a parametric tail to the values of a time-series and estimate quantiles
#' 
#' @param X data sample (double(n))
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
#'       }
#'   \item{if X contains only the n (approximately Poisson) peak values above some threshold 
#'         (in a PoT analysis),  it is recommended to set r11= 1 and take p = f*d/EI and 
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
#'   \item{$nexcess (for PoT only): no. of data values (as opposed to peak values) exceeding the threshold}
#'  }  
#'  
#'  options may contain the following fields:
#'  \itemize{
#'   \item{$pthreshold: fraction of time that value exceeds threshold (double(1))}
#'   \item{$pthresholdmax: upper bound on pthreshold (in case pthreshold is estimated)}
#'   \item{$indexselect: if TRUE, threshold is selected based on tail index estimates (logical, default= FALSE)} 
#'   \item{$kmin: no. of order statistics skipped in determining threshold (integer(1)), default= 20)} 
#'   \item{$sigma: determines the ratio of k to l ( (no. of order stats used for estimation of tail index and quantile) (double(1)}
#'   \item{$fixedpar: fixed model parameters not to be estimated, and their standard errors (list; see below)}
#'   \item{$plotparams: plotparameters (list) with members: $pconf (coverage probability of confidence interval), $xlim (plot limits for quantile estimates), $freqlim (plot limits for frequencies), $plim (plot limits for fractions of time)}
#'  }                                
#'                                        
#' @references
#' De Valk, C. and Cai, J.J. (2018), A high quantile estimator based on 
#' the log-generalized Weibull tail limit. Econometrics and Statistics 6, 107-128, see
#' \url{https://doi.org/10.1016/j.ecosta.2017.03.001}
#
#' @author Cees de Valk \email{ceesfdevalk@gmail.com}
#' 
#' @export
FitTail_AllData <- function(X, freq, df, method, options, metadata) {
  
  # Default parameters
  fixedpar <- sigma <- pthreshold <-  pthresholdmax <- indexselect <- NULL
  
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
    metadata <- list(caseId= Sys.time())
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
  pthresholdmax <- options$pthresholdmax
  if (length(pthresholdmax)< 1) {pthresholdmax <- 0.5} 
  sigma <- options$sigma
  if (length(sigma)< 1) {sigma <- Inf} # to keep behaviour simple to non-expert
  fixedpar <- options$fixedpar
  indexselect <- options$indexselect
  if (length(indexselect)< 1) {indexselect <- FALSE}
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
  EIes <- EI(X, makeplot= TRUE)
  r11es <- r11(X, makeplot= TRUE)
  
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
      Pthresh <- selectThresholdP1(estimates$quantile[, 1], estimates$quantileStd[, 1], estimates$l, 
                                   0.9, kmin= kmin) 
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
    iselect <- min(iselect, max(which(estimates$l< estimates$N*pthresholdmax)))
  }
  
  # Compute quantiles for selected threshold on refined frequency grid 
  ls <- estimates$l[iselect]
  lf <- log10(freq)     # Extend frequency array for plotting etc. 
  mlf <- log10(estimates$l[iselect]/estimates$l*EIvalue*p0/timestep)
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
  plotparams <- options$plotparams
  
  # Plot of tail fit
  tailplot(plotparams, es)
  
  # Plot of tail index estimates vs. l
  # tailindexplot(es= estimates)  
  
  # Plot of quantile estimates vs. l for the lowest freqency
 #  tailquantileplot(plotparams, estimates)   

  return(es)
}
