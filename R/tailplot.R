# Plot of tail estimates for the two datasets, 
# with confidence interval at return period of 10 years 
# params: order, pconf,  
tailplot <- function(params, ...) {
  col <- c("black", "blue", "red", "magenta", "cyan", "grey")
  lwd <- 2
  
  # data of estimates to be plotted 
  es <- list(...)
  les <- length(es)
  
  if (is.null(params)) {
    params <- list()
    params$order <- 0
  }
  if (is.null(params$order)) {
    params$order <- 0
  }
  
  if (!(params$order> 0)) {
    
    if (is.null(params$xlim)) {
      xlim <- c(Inf, -Inf)
      for (i in 1:les) {
        xlim[1] <- min(xlim[1], es[[i]]$location)
        xlim[2] <- max(xlim[2], max(es[[i]]$quantile))
        xlim[2] <- max(xlim[2], max(es[[i]]$orderstats))
      }
      xlim[2] <- xlim[1] + (xlim[2]-xlim[1])*1.1
      xlim <- round(xlim)
      params$xlim <- xlim
    }
    
    if (is.null(params$plim)) {    
      plim <- c(1, 0)
      for (i in 1:les) {
        plim[1] <- min(plim[1], min(es[[i]]$p))
        plim[2] <- max(plim[2], max(es[[i]]$p))
        plim[2] <- max(plim[2], max(exp(-es[[i]]$y)))
      }
      params$plim <- plim
      
      if (is.null(params$pconf)) {
        params$pconf <- 0.9 # default: 90% confidence interval
      }
    }
    
    #
    # plot the successive estimates
    params$order <- 1
    tailplot(params, es[[1]])
    for (i in 2:les) {
      params$order <- i
      tailplot(params, es[[i]])
    }
    return()
  }
  
  # From here on, we only deal with a single dataset/tail estimate
  
  es <- es[[1]]  # unpack from list
  timeseries <- FALSE
  metadata <- es$metadata

  if (!is.null(metadata)) {
    
    varname <- as.character(metadata$varname)
    varunit <- as.character(metadata$varunit)    
    timeunit <- metadata$timeunit
    timestep <- metadata$timestep
    timelength <- metadata$timelength
    nexcess <- metadata$nexcess    
    caseId <- metadata$caseId
    if (is.null(caseId)) {caseId <- Sys.time()}
    
    if (!is.null(timestep)) {
      timeseries <- TRUE
    } else {
      timestep <- 1
    }
    
    EI <- es$EIvalue
    if (is.null(EI)) {
      EI <- 1  # possibly to be modified later if pot
    }
    
    p0 <- es$p0
    if (is.null(p0)) {
      p0 <- 1  
    }    
    
    pot <- FALSE
    if (!is.null(timelength) & timeseries) {
      if (abs(es$N*timestep-timelength*p0) > timestep) {
        pot <- TRUE
        if (EI== 1 & !is.null(nexcess)) 
          EI <- n/nexcess
      } 
    }
    
    if (is.null(timelength) & !is.null(timestep)) {  # this is never for POT!
      timelength <- timestep*es$N/p0
    } else if (!is.null(timelength) & is.null(timestep)) {
      timestep <- timelength/es$N*p0
    }
    
    if (is.null(varunit)) {varunit= "-"}
    
  } else {  # case: metadata= NULL: probabilities instead of frequencies
    
    timeseries <- FALSE
    pot <- FALSE
    varname <- ""
    varunit <- "-"    
    timeunit <- "-" 
    timestep <- 1
    timelength <- es$N
    EI= 1
    
  }
  
  xlim <- params$xlim
  plim <- params$plim
  qn <- abs(qnorm((1-params$pconf)/2)) # half width of normal confidence interval
  
  # axis labels 
  if (is.null(timeseries)) {
    ylab <- "probability of exceedance"
  } else {
    ylab <- paste("frequency of exceedance [/", timeunit, "]", sep= "") 
  }
  xlab <- paste(varname, " [" , varunit, "]", sep= "")
  title <- paste(es$df, "tail", ", sample fraction: ", 
                 signif(es$l/es$N, digits= 2), ", case: ", caseId, sep= "")
  
  # convert p to frequency mu and compute quantile curve
  mu <- es$freq
  if (is.null(mu)) {
    mu <- es$p*EI/timestep
  }
  ylim <- params$freqlim
  if (is.null(ylim)) {
    ylim <- plim*EI/timestep
  }
  
  # plot quantile curve
  i <- (params$order-1)%%length(col)+1
  if (i== 1) {
    par(pty= 's')
    plot(es$quantile, mu, type= "l", log= "y", 
         xlim= xlim, ylim= ylim, col= col[i], lwd= lwd, 
         xlab= xlab, ylab= ylab, main= title,
         xaxp= c(xlim, diff(xlim)), yaxp= c(ylim, 3))
    # grid(ny= diff(ylim))
    grid()
  } else {
    lines(es$quantile, mu, col= col[i], lwd= lwd)
  }
  one= 1.
  b <- es$quantile[1]
  points(b, mu[1]*one, col= col[i], type="o", cex= 1)
  lines(c(b-es$quantileStd[1]*qn, b+es$quantileStd[1]*qn), rep(mu[1]*one, 2), 
        col=col[i], type="o", pch="|", cex= 2, lwd= lwd)
  
  # plot sample curve 
  X <- -sort(-es$orderstats)
  muX <- (1:length(X))/timelength
  if (!pot) {muX <- muX*EI}
  points(X, muX, pch= 20, col= col[i])
  lines(X, muX, col= col[i], type= 's')

} # klaar
