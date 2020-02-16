#
# 
#
assigncat <- function(y, lbin, ubin) {
  ly <- length(y)
  nbin <- length(lbin)
  cats <- NA
  cat <- data.matrix(rep(TRUE, ly))
  if (nbin> 1) {
    if (length(ubin)!= nbin) {
      stop("lbin and ubin have different lengths in $options$covariate.")
    } else {
      cats <- rep(NA, nbin)
      cat <- rep(FALSE, ly*nbin)
      dim(cat) <- c(ly, nbin)
      for (i in 1:nbin) {
        if (ubin[i]> lbin[i]) {
          cat[, i] <- y>= lbin[i] & y< ubin[i]
          cats[i] <- (ubin[i]+lbin[i])/2
        } else {
          cat[, i] <- y>= lbin[i] | y< ubin[i]  
          cats[i] <- (((ubin[i]+lbin[i]-360)/2-.01) %% 360) + .01
        }
      }
    }
  }
  value <- list(cat= cat, cats= cats, lbin= lbin, ubin= ubin)
}

  