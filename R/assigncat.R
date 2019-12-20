#
# 
#
assigncat <- function(y, lbin, ubin, cats, binwidth) {
  ly <- lenth(y)
  nbin <- length(lbin)
  
  if (length(binwidth)> 0 & (nbin== 0 & length(ubin)== 0)) {
    nbin <- ceiling((max(y)-min(y))/binwidth)
    cat <- (((y+binwidth/2) %/% binwidth) %% nbin)*binwidth
    cats <- sort(unique(cat))
    lbin <- cats-binwidth/2
    ubin <- cats+binwidth/2
    period <- ubin[nbin]-lbin[1]
    lbin <- lbin %% period    
    ubin <- ubin %% period  
    ind <- sort(lbin, index.return= TRUE)$ix
    lbin <- lbin[ind]
    ubin <- ubin[ind]
    cats <- cats[ind]
    cat0 <- cat
    cat <- rep(FALSE, ly*lbin)
    dim(cat) <- c(ly, lbin)
    cat[cbind(1:ly, cat0)] <- TRUE
  } else if (nbin> 0 & length(ubin)== nbin) {
    if (length(cats)!= nbin) {
       cats <- lbin
       warning("Bins (re)labelled by their lower bounds.")
    }
    if (nbin> 1) {
      ibin <- sort(lbin, index.return= TRUE)$ix
      lbin <- lbin[ibin]
      ubin <- ubin[ibin]
      cats <- cats[ibin]
    } else {
      cats <- NA
    }
    cat <- rep(FALSE, ly*lbin)
    dim(cat) <- c(ly, lbin)
    for (i in 1:nbin) {
       if (ubin[i]> lbin[i]) {
         cat[, i] <- y>= lbin[i] & y< ubin[i]
       } else {
         cat[, i] <- y>= lbin[i] | y< ubin[i]         
       }
    }
  } else {  
    stop("Inconsistencies/omissions in options$covariate found.")
  }
  value <- list(cat= cat, cats= cats, lbin= lbin, ubin= ubin)
}

