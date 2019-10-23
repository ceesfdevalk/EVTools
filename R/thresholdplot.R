#' @name  thresholdplot
#' 
#' @title thresholdplot
#' 
#' @description # Plot of p-value of a fluctuation size statistic, for threshold selection  
#' 
#' @param params (optional) list (see below)
#' @param es list containing tail estimates from a single sample 
#' 
#' @usage thresholdplot(params, es)
#' 
#' @return A plot file (.png)
#
#' @author Cees de Valk \email{ceesfdevalk@gmail.com}
#' 
#' @export  
thresholdplot <- function(params= NULL, es= NULL) {
  lwd <- 2
  metadata <- es$metadata
  caseId <- metadata$caseId
  if (is.null(caseId)) {caseId <- Sys.time()}
  varname <- as.character(metadata$varname)
  
  # axis labels 
  ylab <- paste(es$df, "p-value")
  xlab <- paste("sample fraction for quantile estimation")
  title <- paste(es$df, " threshold p-value", ", case: ", caseId, sep= "")
  
  ylim <- c(0, 1)
  xlim <- c(10^floor(log10(min(es$l)/es$N)), 1)
  lP <- approx(es$threshold$k, es$k, es$l, rule= 1)$y

  # plot
  par(pty= 's')
  plot(lP/es$N, es$threshold$P, type= "l", log= "x", 
       xlim= xlim, ylim= ylim,  lwd= lwd, 
       xlab= xlab, ylab= ylab, main= title,
       yaxp= c(ylim, diff(ylim)*10)) #, tck = 1)
  
  if (length(es$selected)> 0) {
    lines(es$selected$l/es$N*c(1,1), ylim, lty= 2) 
    # points(es$selected$l/es$selected$N, es$selected$tailindex, lty= 2)
  }
  grid()
} # klaar
