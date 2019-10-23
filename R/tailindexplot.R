#' @name  tailindexplot
#' 
#' @title tailindexplot
#' 
#' @description # Plot of tail index estimates with confidence interval  
#' 
#' @param params (optional) list (see below)
#' @param estimates list containing tail estimates from a single sample 
#' 
#' @usage tailindexplot(params, estimates)
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
tailindexplot <- function(params= NULL, estimates= NULL) {
  lwd <- 2
  metadata <- estimates$metadata
  caseId <- metadata$caseId
  if (is.null(caseId)) {caseId <- Sys.time()}
  varname <- as.character(metadata$varname)
  
  pconf <- params$pconf
  if (is.null(pconf)) {pconf <- 0.9}
  qn <- abs(qnorm((1-pconf)/2)) # half width of normal confidence interval
  
  # axis labels 
  ylab <- paste(estimates$df, "tail index")
  xlab <- paste("sample fraction for quantile estimate")
  title <- paste(estimates$df, " tail index", ", case: ", caseId, sep= "")
  
  med <- median(estimates$tailindex[estimates$l< 0.1*estimates$N])
  ylim <- 0.5*round(med/0.5)+c(-1, 1)
  
  print("ylim:")
  print(ylim)
  
  xlim <- c(10^floor(log10(min(estimates$l)/estimates$N)), 1)
  
  # plot
  par(pty= 's')
  plot(estimates$l/estimates$N, estimates$tailindex, type= "l", log= "x", 
       xlim= xlim, ylim= ylim,  lwd= lwd, 
       xlab= xlab, ylab= ylab, main= title)
  lines(estimates$l/estimates$N, estimates$tailindex + estimates$tailindexStd*qn, lwd= 1) 
  lines(estimates$l/estimates$N, estimates$tailindex - estimates$tailindexStd*qn, lwd= 1) 
  grid()
} # klaar
