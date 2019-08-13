#' @name  r11
#' 
#' @title r11
#' 
#' @description Computes estimates of the extremal dependence coeficient r(1,1) of 
#' Rootzen (1995)
#' 
#' @param X time series (double(nX))
#' @param l (optional, defailt= 50) no of data points above highest level (integer(1))
#' @param ngr (optional, default= 10) number of levels (integer(1))
#' @param makeplot (optional, default= FALSE) (logical(1))
#' 
#' @usage Value <- r11(X, l, ngr, makeplot) 
#' 
#' @return A list, with members: 
#'   \item{k}{no. of data points above level, for each of ngr levels}
#'   \item{p}{sample fraction above level, for each of ngr levels} 
#'   \item{r}{corresponding estimate of r(1,1)} 
#'
#' @references
#'  Rootzen, H. (1995), The tail empirical process for stationary sequences. Technical report,
#'  Department of Mathematics, Chalmers University, Sweden.
#'
#' @author Cees de Valk \email{ceesfdevalk@gmail.com}
#' 
#' @export
r11 <- function(X, l= 50, ngr= 20, makeplot= FALSE) {
  sX <- -sort(-X)
  n <- length(X)
  gr <- seq(log(l), log(n), length.out= ngr)
  pl <- matrix(0, nrow= l, ncol= ngr)
  k <- ceiling(exp(gr))
  p= k/n
  ll <- min(l, n/k) # l must never be too large (actually o(n/k))
  for (i in 1:ngr) {
    cat(i)
    s <- sX[k[i]] 
    for (j in 1:ll[i]) {
      id2 <- X[1:(n-l)]> s & X[(j+1):(n-l+j)]> s
      id1 <- X[1:(n-l)]> s
      pl[j, i] <- sum(id2)/sum(id1) - sum(id1)/(n-l) # last term is mean; small in tail
    }
    pl[, i] <- pl[, i]*(1-(1:ll[i])/ll[i])
  }
  # plot((1:l),cumsum((1-(1:l)/l)*pl))
  r <- 1+2*colSums(pl)
  
  if (makeplot) {
    par(pty = "s")
    plot(p, 1/r, log= 'x', ylim= c(0, 1), 
         main= "Serial tail dependence", xlab= "sample fraction", ylab= '1/r(1,1)')
    grid()
  }
  
  res <- list(k= k, p= p, r= r)
}