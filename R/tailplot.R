#' @name  tailplot
#' 
#' @title tailplot
#' 
#' @description # Plot of tail estimates for one or several samples, with confidence interval  
#' 
#' @param ... one or several lists es1, es2, ..., each containing a tail estimate to be plotted
#' @param params (optional) list (see below)
#' 
#' @usage  tailplot(es1, es2, ..., params= params)
#' 
#' @return A plot file (.png)
#' 
#' @details The parameter list params must be explicitly specified with "params= ..." and may contain:
#'  \itemize{
#'   \item{$pconf: coverage probability of confidence interval (0.9 by default) (double(1))}
#'   \item{$xlim: plot limits for quantile estimates (double(2))}
#'   \item{$freqlim: plot limits for frequency (double(2))}
#'   }
#
#' @author Cees de Valk \email{ceesfdevalk@gmail.com}
#' 
#' @export
tailplot <- function(..., params) {
  col <- c("black", "blue", "red", "magenta", "cyan", "grey")
  lwd <- 2
  
  # data of estimates to be plotted 
  es <- list(...)
  les <- length(es)
  
  for (j in 1:les) {
    if (!is.null(es[[j]]$selected)) {
      es[[j]] <- es[[j]]$selected
    }
  }
  
  if (missing(params)) {
    params <- NULL
  }
  if (length(params)< 1) {
    params <- list()
    params$order <- 0
  }
  if (is.null(params$order)) {
    params$order <- 0
  }
  
  if (params$order<= 0) {
    
    if (is.null(params$xlim)) {
      xlim <- c(Inf, -Inf)
      for (i in 1:les) {
        xlim[1] <- min(xlim[1], es[[i]]$location, na.rm= TRUE)
        xlim[2] <- max(xlim[2], max(es[[i]]$quantile, na.rm= TRUE), na.rm= TRUE)
        xlim[2] <- max(xlim[2], max(es[[i]]$orderstats, na.rm= TRUE), na.rm= TRUE)
      }
      xlim[2] <- xlim[1] + (xlim[2]-xlim[1])*1.1
      xlim <- round(xlim)
      params$xlim <- xlim
    }
    
    if (is.null(params$freqlim)) {    
      freqlim <- c(1, 0)
      for (i in 1:les) {
        freqlim[1] <- min(freqlim[1], min(es[[i]]$freq, na.rm= TRUE), na.rm= TRUE)
        freqlim[2] <- max(freqlim[2], max(es[[i]]$freq, na.rm= TRUE), na.rm= TRUE)
      }
      params$freqlim <- freqlim
    }
    
    if (length(params$pconf)< 1) {
      params$pconf <- 0.9 # default: 90% confidence interval
    }
    
    #
    # plot the successive estimates
    params$order <- 1
    tailplot(es[[1]], params= params)
    if (les> 1) {
      for (i in 2:les) {
        params$order <- i
        tailplot(es[[i]], params= params)
      }
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
    
    pbin<- es$pbin
    if (is.null(pbin)) {
      pbin <- 1  
    }    
    
    pot <- FALSE
    if (!is.null(timelength) & timeseries) {
      if (abs(es$N*timestep-timelength*pbin) > timestep) {
        pot <- TRUE
        if (EI== 1 & !is.null(nexcess)) 
          EI <- n/nexcess
      } 
    }
    
    if (is.null(timelength) & !is.null(timestep)) {  # this is never for POT!
      timelength <- timestep*es$N/pbin
    } else if (!is.null(timelength) & is.null(timestep)) {
      timestep <- timelength/es$N*pbin
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
    caseId <- Sys.time()
    
  }
  
  xlim <- params$xlim
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
  
  # plot quantile curve
  i <- (params$order-1)%%length(col)+1
  if (i== 1) {
    par(pty= 's')
    plot(es$quantile, mu, type= "l", log= "y", 
         xlim= xlim, ylim= ylim, col= col[i], lwd= lwd, 
         xlab= xlab, ylab= ylab, main= title,
         xaxp= c(xlim, diff(xlim)), yaxp= c(ylim, 3)) #, tck = 1)
    # grid(ny= diff(ylim))
    grid()
  } else {
    lines(es$quantile, mu, col= col[i], lwd= lwd)
  }
  
  one= 1.
  b <- es$quantile[1]
  points(b, mu[1]*one, col= col[i], type="o", cex= 1)
  lines(c(b-es$quantileStd[1]*qn, b+es$quantileStd[1]*qn), rep(mu[1]*one, 2), 
        col=col[i], pch="|", cex= 0.5, lwd= 0.5)
  lines(es$quantile-es$quantileStd*qn, mu, col= col[i], lwd= 0.5)
  lines(es$quantile+es$quantileStd*qn, mu, col= col[i], lwd= 0.5)
  
  # plot sample curve 
  X <- -sort(-es$orderstats)
  lX <- length(X)
  muX <- (1:lX)/timelength
  if (!pot) {muX <- muX*EI}
  id <- unique(round(exp((0:1.e4)*log(lX)*1.e-4)))
  points(X[id], muX[id], pch= 20, col= col[i])
  lines(X[id], muX[id], col= col[i], type= 's')
  
} # klaar
