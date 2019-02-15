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
r11 <- function(X, l= 50, ngr= 10, makeplot= FALSE) {
  sX <- -sort(-X)
  n <- length(X)
  gr <- seq(log(l), log(n/5), length.out= ngr)
  pl <- matrix(NA, nrow= l, ncol= ngr)
  k <- ceiling(exp(gr))
  p= k/n
  for (i in 1:ngr) {
    cat(i)
    s <- sX[k[i]] 
    for (j in 1:l) {
      id2 <- X[1:(n-l)]> s & X[(j+1):(n-l+j)]> s
      id1 <- X[1:(n-l)]> s
      pl[j, i] <- sum(id2)/sum(id1)
    }
    pl[, i] <- pl[, i]*(1-(1:l)/l)
  }
  # plot((1:l),cumsum((1-(1:l)/l)*pl))
  r <- 1+2*colSums(pl)
  
  if (makeplot) {
    par(pty = "s")
    plot(p, 1/r, log= 'x', ylim= c(0, 1), 
         main= "Serial tail dependence", xlab= "rank no.", ylab= '1/r(1,1)')
    grid()
  }
  
  res <- list(k= k, p= p, r= r)
}