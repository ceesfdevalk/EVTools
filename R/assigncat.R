#
# 
#
assigncat <- function(y, lbin, ubin) {
  ly <- length(y)
  nbin <- length(lbin)
  cat <- data.matrix(rep(TRUE, ly))
  if (nbin> 1) {
    if (length(ubin)!= nbin) {
      stop("lbin and ubin have different lengths in $options$covariate.")
    } else {
      cat <- rep(FALSE, ly*nbin)
      dim(cat) <- c(ly, nbin)
      for (i in 1:nbin) {
        if (ubin[i]> lbin[i]) {
          cat[, i] <- y>= lbin[i] & y< ubin[i]
        } else {
          cat[, i] <- y>= lbin[i] | y< ubin[i]  
        }
      }
    }
  }
  value <- list(cat= cat, lbin= lbin, ubin= ubin)
}

  