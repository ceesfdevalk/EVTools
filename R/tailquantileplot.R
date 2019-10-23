#' @name  tailquantileplot
#' 
#' @title tailquantileplot
#' 
#' @description # Plot of tail quantile estimates for a single frequency with confidence interval  
#' 
#' @param params (optional) list (see below)
#' @param estimates list containing tail estimates from a single sample 
#' 
#' @usage tailquantileplot(params, estimates)
#' 
#' @return A plot file (.png)
#' 
#' @details The parameter list params may contain:
#'  \itemize{
#'   \item{$pconf: coverage probability of confidence interval (0.9 by default) (double(1))}
#'   }
#
#' @author Cees de Valk \email{ceesfdevalk@gmail.com}
#' 
#' @export  
tailquantileplot <- function(params= NULL, es= NULL) {
  if (length(es)< 1) {stop("Need data to plot.")}
  lwd <- 2
  metadata <- es$metadata
  varname <- as.character(metadata$varname)
  varunit <- as.character(metadata$varunit)
  timeunit <- metadata$timeunit
  if (is.null(timeunit)) {timeunit= "-"}
  caseId <- metadata$caseId
  if (is.null(caseId)) {caseId <- Sys.time()}
  
  pconf <- params$pconf
  if (is.null(pconf)) {pconf <- 0.9}
  qn <- abs(qnorm((1-pconf)/2)) # half width of normal confidence interval
  
  # axis labels 
  ylab <- paste(varname, " [", varunit, "]", sep= "")
  xlab <- paste("sample fraction for quantile estimate")
  title <- paste(es$df, " quantile", " at freq. ", freq[1], "/", timeunit,  
                 ", case: ", caseId, sep= "")
  
  freq <- es$freq[1]
  q <- es$quantile[, 1]
  qStd <- es$quantileStd[, 1]
  id <- es$l< 0.1*es$N
  ylim <- c(quantile(q[id]-qStd[id]*qn, 0.05), quantile(q[id]+qStd[id]*qn, 0.95))
  ylim <- signif(ylim, digits= 2)       
  
  print(ylim)
  xlim <- c(10^floor(log10(min(es$l)/es$N)), 1)
  
  # plot
  par(pty= 's')
  plot(es$l/es$N, q, type= "l", log= "x", 
       xlim= xlim, ylim= ylim,  lwd= lwd, 
       xlab= xlab, ylab= ylab, main= title,
       yaxp= c(ylim, diff(ylim))) #,  tck = 1)
  lines(es$l/es$N, q + qStd*qn, lwd= 1) 
  lines(es$l/es$N, q - qStd*qn, lwd= 1) 
  if (length(es$selected)> 0) {
    lines(es$selected$l/es$selected$N*c(1,1), ylim, lty= 2) 
    # points(es$selected$l/es$selected$N, es$selected$quantile[1], lty= 2)
  }
  
  grid()
} # klaar
