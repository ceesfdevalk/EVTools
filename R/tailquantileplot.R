# Plot of tail index estimates 
# with confidence interval
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
  qn <- abs(qnorm((1-params$pconf)/2)) # half width of normal confidence interval
  
  # axis labels 
  ylab <- paste(es$df, "tail quantile")
  xlab <- paste("sample fraction for quantile estimate")
  title <- paste(es$df, " quantile of ", varname, " at freq. ", freq[1], "/", timeunit,  
                 ", case: ", caseId, sep= "")
  
  freq <- es$freq[1]
  q <- es$quantile[, 1]
  qStd <- es$quantileStd[, 1]
  id <- es$l< 0.1*es$N
  ylim <- c(quantile(q[id]-qStd[id]*qn, 0.1), quantile(q[id]+qStd[id]*qn, 0.9))
  ylim <- signif(ylim, digits= 2)       
  
  print(ylim)
  xlim <- c(10^floor(log10(min(es$l)/es$N)), 1)
  
  # plot
  par(pty= 's')
  plot(es$l/es$N, q, type= "l", log= "x", 
       xlim= xlim, ylim= ylim,  lwd= lwd, 
       xlab= xlab, ylab= ylab, main= title)
  lines(es$l/es$N, q + qStd*qn, lwd= 1) 
  lines(es$l/es$N, q - qStd*qn, lwd= 1) 
  grid()
} # klaar
