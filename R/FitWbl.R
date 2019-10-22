#' @name  FitWbl
#' 
#' @title FitWbl
#' 
#' @description Fit a Weibull upper tail to the sample X and estimate quantiles
#' 
#' @param X data sample (double(n))
#' @param method (optional): "FitWbl_MLE" or (as Wbl is special case of log-GW) "FitGW_iHilli", "FitGW_Mom", "FitGW_MLE"
#' @param p (optional) probabilities of exceedance of the quantiles to be estimated (double(np))  
#' @param N (optional) (effective) sample size, in case X is not complete but contains only (peak) values above some threshold (integer(1))
#' @param r11 (optional) factor to increase estimator variance by, to account for serial dependence (default: 1) (double(1) or list, see Details)
#' @param fixedpar (obsolete for Weibull)
#' @param l0 (optional) value of l (no. of order stats used) in case it is imposed (integer(0))
#' @param sigma (obsolete for Weibull)
#' @param metadata (optional) information about the variable and, if applicable, the time-series (list; see Details)
#' 
#' @usage Value <- FitWbl(X, method= "FitGW_iHilli", p= NULL, N= 0, r11= 1, fixedpar= NULL, l0= NULL, sigma= Inf, metadata= NULL)
#' 
#' @return A list, with members: 
#'   \item{l}{no. of order statistics used for scale and quantile estimation}    
#'   \item{k}{no. of order statistics used for tail index estimation} 
#'   \item{sigma}{= Inf}
#'   \item{tailindex}{estimates or imposed value of Weibull tail index} 
#'   \item{tailindexStd}{standard deviations of tail index estimates}
#'   \item{locationStd}{standard deviation of order statistic}
#'   \item{lambda}{ratio of logarithms of probabilities of exceedance of quantile and threshold}  
#'   \item{p}{probabilities of exceedance of quantiles to be estimated} 
#'   \item{quantile}{quantile estimates}
#'   \item{quantileStd}{standard deviations of quantile estimates}
#'   \item{tailindexraw}{raw estimates of Weibull tail index over all possible thresholds (method: FitGW_iHill.R)} 
#'   \item{tailindexrawStd}{standard deviation of tailindexraw}
#'   \item{kraw}{no. of order statistics used for estimation of tailindexraw} 
#'   \item{orderstats}{data X sorted (decreasing)}
#'   \item{df}{= "Weibull": fitted distribution function tail (Generalised Weibull}
#'   \item{estimator}{= "iteratedHill": see "method" below}
#' 
#' @details
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
#'   \item{$EI: extremal index (see above)}
#'   \item{$nexcess (for PoT only): no. of data values (as opposed to peak values) exceeding the threshold}
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
FitWbl <- function(X, method, p, N, r11, fixedpar, l0, sigma, metadata) {
  
  # Handle arguments
  if (missing(method)) {method <- "FitGW_iHilli"}
  if (missing(p)) {p <- NULL}
  if (missing(N)) {N <- 0} 
  if (missing(r11)) {r11 <- 1}
  if (missing(l0)) {l0 <- NULL}
  if (missing(sigma)) {sigma <- Inf}
  if (missing(metadata)) {metadata <-NULL}
  
  known <- grep(method, c("FitWbl_MLE", "FitGW_iHilli", "FitGW_Mom", "FitGW_MLE"))
  if (length(known)< 1) {
    stop("Specified method is unknown to FitWbl.R")
  }
  
  sigma <- Inf #sigma not applicable to 1-parameter Weibull tail estimation (has no effect anyway)
  fixedpar <- NULL
  
  if (method== "FitWbl_MLE") {  
    
    fixedpar <- list(f0= 0, f0Std= 0) # offset fixed to zero for standard Weibull model
    estimates <- FitWbl_MLE(X, p, N, r11, fixedpar, l0, sigma, metadata)
    
  } else {           # estimating Weibull tail as special case of log-GW tail
    
    fixedpar <- list(theta0= 0, theta0Std= 0)
    estimates <- FitlogGW(X, method, p, N, r11, fixedpar, l0, sigma, metadata)
    
    estimates$tailindex <- estimates$scale
    estimates$tailindexStd <- estimates$logdispStd*estimates$tailindex 
    estimates$scale <- NULL
    estimates$logdisp <- NULL
    estimates$logdispStd <- NULL
    # comment: location parameter and quantile of log-GW are already 
    # good (see FitlogGW.R)
  }
  return(estimates)
}

