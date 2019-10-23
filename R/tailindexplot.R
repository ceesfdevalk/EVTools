#' @name  tailindexplot
#' 
#' @title tailindexplot
#' 
#' @description # Plot of tail index estimates with confidence interval  
#' 
#' @param params (optional) list (see below)
#' @param es list containing tail estimates from a single sample 
#' 
#' @usage tailindexplot(params, es)
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
tailindexplot <- function(params= NULL, es= NULL) {
  lwd <- 2
  metadata <- es$metadata
  caseId <- metadata$caseId
  if (is.null(caseId)) {caseId <- Sys.time()}
  varname <- as.character(metadata$varname)
  
  pconf <- params$pconf
  if (is.null(pconf)) {pconf <- 0.9}
  qn <- abs(qnorm((1-pconf)/2)) # half width of normal confidence interval
  
  # axis labels 
  ylab <- paste(es$df, "tail index")
  xlab <- paste("sample fraction for quantile estimate")
  title <- paste(es$df, " tail index", ", case: ", caseId, sep= "")
  
  med <- median(es$tailindex[es$l< 0.1*es$N])
  ylim <- 0.5*round(med/0.5)+c(-1, 1)
  xlim <- c(10^floor(log10(min(es$l)/es$N)), 1)
  
  # plot
  par(pty= 's')
  plot(es$l/es$N, es$tailindex, type= "l", log= "x", 
       xlim= xlim, ylim= ylim,  lwd= lwd, 
       xlab= xlab, ylab= ylab, main= title,
       yaxp= c(ylim, diff(ylim)*10)) #, tck = 1)
  lines(es$l/es$N, es$tailindex + es$tailindexStd*qn, lwd= 1) 
  lines(es$l/es$N, es$tailindex - es$tailindexStd*qn, lwd= 1) 
  if (length(es$selected)> 0) {
    lines(es$selected$l/es$selected$N*c(1,1), ylim, lty= 2) 
    # points(es$selected$l/es$selected$N, es$selected$tailindex, lty= 2)
  }
  grid()
} # klaar
