#' @name  EI
#' 
#' @title EI
#' 
#' @description Computes estimates of the Extremal Index
#' 
#' 
#' @param X time series (double(nX))
#' @param l (optional, defailt= 50) no of data points above highest level (integer(1))
#' @param ngr (optional, default= 10) number of levels (integer(1))
#' @param makeplot (optional, default= FALSE) (logical(1))
#' 
#' @usage Value <- EI(X, l, ngr, makeplot) 
#' 
#' @return A list, with members: 
#'   \item{k}{no. of data points above level, for each of ngr levels}    
#'   \item{p}{sample fraction above level, for each of ngr levels}    
#'   \item{EIupcross}{estimate of EI based on upcrossings} 
#'   \item{EIFS}{estimate of EI from Ferro & Segers (2003)} 
#'
#' @references
#'  M. R. Leadbetter, M.R., Lindgren, G. and Rootzen, H. (1983), 
#'  Extremes and Related Properties of Random Sequences and Processes. 
#'  Springer.
#'  Ferro, A.T. and J. Segers (2003), Inference for clusters 
#'  of extreme values. J. R. Statist. Soc. B (2003) 65, 545–556.
#'
#' @author Cees de Valk \email{ceesfdevalk@gmail.com}
#' 
#' @export
EI <- function(X, l= 50, ngr= 20, makeplot= FALSE) {
  sX <- -sort(-X)
  n <- length(X)
  gr <- seq(log(l), log(n), length.out= ngr)
  k <- ceiling(exp(gr)) 
  p= k/n

  EIupcross <- rep(NA, ngr)
  EIFS <- rep(NA, ngr)
  
  for (i in 1:ngr) {
    # cat(i)
    s <- sX[k[i]] 

    # simple upcrossing count    
    id2 <- X[1:(n-1)]<= s & X[2:n]> s
    id1 <- X[2:n]> s
    EIupcross[i] <- sum(id2)/sum(id1)
    
    # F&S estimator
    S <- which(X> s)
    if (length(S)> 1) {
      T <- diff(S)
      sT <- sum(T)
      lT <- length(T)
      if (max(T)> 2) {
        theta <- 2*(sT-lT)^2/(sum((T-1)*(T-2))*lT)
      } else {
        theta <- 2*sT^2/(sum(T^2)*lT)
      }
      EIFS[i] <- pmin(1, theta)
    }
  }
  
  if (makeplot) {
    par(pty = "s")
    title <- "Extremal Index naive: circle; Ferro-Segers: triangle"
    plot(p, EIupcross, log= 'x', ylim= c(0, 1), 
         main= title, xlab= "sample fraction", ylab= "Extremal Index")
    points(p, EIFS, pch= 2)
    grid()
  }
  
  res <- list(k= k, p= p, EIupcross= EIupcross, EIFS= EIFS)
}