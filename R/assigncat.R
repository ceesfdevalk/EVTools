#
# 
#
assigncat <- function(y, lbin, ubin, binwidth) {
  
  nbin <- length(lbin)
  if (length(binwidth)> 0 & (nbin== 0 | length(ubin)== 0 | nbin!= length(ubin))) {
    
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
    
  } else if (nbin> 0 & length(ubin)== nbin) {
    
    if (nbin> 1) {
      if (!all(diff(lbin)> 0)) {
        stop("options$covariate$lbin must be increasing.")
      }
      if (nbin> 2) {
        if (!all(diff(ubin[1:(nbin-1)])> 0)) {
          stop("options$covariate$ubin[1:(nbin-1)] must be increasing.")
        }
      }
    }
    cats <- apply(rbind(lbin, ubin), 2, mean)
    lbin1 <- lbin
    ubin1 <- ubin
    cats1 <- cat
    ind <- which(ubin< lbin)
    if (length(ind)== 1) {
      if (ind!= nbin) {
        stop("If any(ubin< lbin), it must be the last element.")
      }
      if (any(cats< 0)) {
        stop("lbin and ubin must be nonnegative if any(ubin< lbin).")
      }
      cats[nbin] <- 0    # simple choice: apparently it is circular 
      lbin1 <- c(-Inf, lbin)
      ubin1 <- c(ubin[nbin], ubin[1:(nbin-1)], Inf)
      cats1 <- c(cats[nbin], cats)  # simple choice: apparently it is circular 
    } else if (length(ind)> 0) { 
      stop("options$covariate$ubin< options$covariate$lbin for at most one element.")
    }
    # period <- ubin[nbin]-lbin[1]
    b <- rbind(lbin1, ubin1)
    b <- as.vector(b)
    fn <- rbind(cats1, rep(Inf, (nbin+1)))
    fn <- as.vector(fn)
    cat <- approx(b, fn, y, method <- "constant", yleft= Inf,
                  yright= Inf, f= 0, ties= "ordered")$y
    
  }
  value <- list(cat= cat, cats= cats, lbin= lbin, ubin= ubin)
}

