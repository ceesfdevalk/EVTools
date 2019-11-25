#' @name  CovariateGraph
#' 
#' @title CovariateGraph
#' 
#' @description Conditional probability of a covariate belonging to some bin, given that the primary variate exceeds a threshold
#' 
#' @param X            Data of primary variate X[, 1] and covariate X[, 2] (double(n, 2))
#' @param nbin         (optional) No. of bins: range of covariate is subdivided by nbin bins and covariate is rounded
#'                     to centre of bin (double(1))
#'                     
#' @usage Value <- CovariateGraph(X, nbin)
#' 
#' @return list of length nbin, containing for each bin and each value x of the primary variate 
#'         coinciding with a covariate value in this bin, the number of times that the primary 
#'         variate is larger than or equal to x
#'
#' @author Cees de Valk \email{ceesfdevalk@gmail.com}
#' 
#' @details If no nbin is specified, then the covariate is assumed to be discrete, and the
#'          bins are its range of values 
#'          
#'          For the i-th bin, the graph of its conditional probability p2 as a function of the probability 
#'          p1 of exceedance of the value of the primary variate is given by p2 <- (1:length(l))/l, 
#'          p1 <- l/n with n <- dim(X)[1] and l <- Value[[i]]
#' 
#' @export
CovariateGraph <- function(X, nbin) {
  X <- data.matrix(X)
  dd <- dim(X)
  if (dd[2]< 2) {
    stop("X must have (at least) 2 columns")
  }
  l <- 1:dd[1]
  
  if (!missing(nbin)) {
    if (nbin> 1) { 
      cat <- X[, 2]
      binw <- (max(cat)-min(cat))/nbin
      dsa <- diff(sort(cat))
      delta <- min(dsa[dsa> 0])
      binw <- ceil(binw/delta)*delta
      cat <- (((cat+binw/2) %/% binw) %% nbin)*binw
    }
  }
  
  temp <- sort(-X[, 1], index.return= TRUE) 
  j <- temp$ix
  
  cat <- cat[j]  # category values corresponding to highest X[, 1], highest but one X[, 1], etc.
  cats <- sort(unique(cat))
  lcats <- length(cats)
  
  numbers <- vector(mode= "list", length= lcats)
  for (i in 1:lcats) {
    ind <- which(cat== cats[i])
    if (length(ind)> 0) {
      numbers[[i]] <- l[ind]
    }
  } # for
  
  par(pty= "s")
  plot(1, 1, type= "l", ylim= c(1/dd[1], 1), xlim= c(1/dd[1], 1), log= "xy",
       xlab= "fraction exceeded", ylab= "cond. prob. of category")
  grid()
  
  for (i in 1:lcat) {
    l <- numbers[[i]]
    ll <- length(l)
    i1 <- unique(round(exp((0:4)*max(log(ll))/4)))
    l1 <- l[i1]
    if (ll> 0) {
      lines(l/dd[1], (1:length(l))/l)
      # points(l1/dd[1], i1/l1, i-1, cex= 1.25)
      text(l1/dd[1], i1/l1, as.character(cats[i]), cex= .8, font= 2)
    }
  } # fo 
  
  result <- numbers
}
