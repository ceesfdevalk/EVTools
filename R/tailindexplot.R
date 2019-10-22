# Plot of tail estimates for the two datasets, 
# with confidence interval at return period of 10 years 
# params: order, pconf,  
tailindexplot <- function(params= NULL, es= NULL) {
  lwd <- 2
  metadata <- es$metadata
  caseId <- metadata$caseId
  if (isnull(caseId)) {caseId <- Sys.time()}
  varname <- as.character(metadata$varname)
  
  pconf <- params$pconf
  if (is.null(pconf)) {pconf <- 0.9}
  qn <- abs(qnorm((1-params$pconf)/2)) # half width of normal confidence interval
  
  # axis labels 
  ylab <- paste(es$df, "tail index")
  xlab <- paste("sample fraction for quantile estimate")
  
  med <- median(es$tailindex[es$l< 0.1*es$N])
  ylim <- med+c(-1, 1)
  xlim <- c(10^floor(log10(min(es$l)/es$N)), 1)
  
  # plot
  par(pty= 's')
  plot(es$l/es$N, es$tailindex, type= "l", log= "x", 
       xlim= xlim, ylim= ylim,  lwd= lwd, 
       xlab= xlab, ylab= ylab)
  lines(es$l/es$N, es$tailindex + es$tailindexStd*qn, lwd= 1) 
  lines(es$l/es$N, es$tailindex - es$tailindexStd*qn, lwd= 1) 
  grid()
} # klaar