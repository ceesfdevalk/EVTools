#'  Pre-determined model parameters are to be supplied in the list fixedpar (see above):
#'  \itemize{
#'   \item{$theta0: (optional) value of tailindex in case it is imposed (double(1))}
#'   \item{$theta0Std: (optional) its standard deviation (double(1))}
#'   \item{$logdisp0: (optional) value of log of dispersion coeff. in case it is imposed (dispersion coeff. is the raio of scale par. to location par.) (double(1))}
#'   \item{$logdisp0Std: (optional) its standard deviation (double(1))}        
#'  }
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
#'   \item{$dither: width of uniform distribution of noise to add to data (double(1))}
#'   \item{$pthreshold: fraction of time that value exceeds threshold (double(1))}
#'   \item{$maxpthreshold: upper bound on pthreshold (in case pthreshold is estimated)}
#'   \item{$minpthreshold: lower bound on pthreshold (in case pthreshold is estimated)}
#'   \item{$indexselect: if TRUE, threshold is selected based on tail index estimates (logical, default= FALSE)} 
#'   \item{$kmin: no. of order statistics skipped in determining threshold (integer(1)), default= 20)} 
#'   \item{$sigma: determines the ratio of k to l ( (no. of order stats used for estimation of tail index and quantile) (double(1)}
#'   \item{$bootstrap: list. If exists/nonempty, precision is assessed by a moving block bootstrap. May contain $nsamples (no. of bootstrap samples) and $blocktime (block length in terms of time)}
#'   \item{$fixedpar: fixed model parameters not to be estimated, and their standard errors (list; see below)}
#'   \item{$plotparams: plotparameters (list) with members: $makeplot (default= TRUE), $pconf (coverage probability of confidence interval), $xlim (plot limits for quantile estimates), $freqlim (plot limits for frequencies), $plim (plot limits for fractions of time)}
#'  }            
#'  
#'  