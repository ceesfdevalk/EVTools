#' @name  PoTselect
#' 
#' @title PoTselect
#' 
#' @description Select a Peak over Threshold (PoT) sample from a univariate time series 
#' 
#' @param X            time series (sequence) (double(n))
#' @param p            fraction of elements of X exceeding the threshold (double(1))
#' @param separation   minimum distance between selected peaks in the sequence (integer(1))  
#'  
#' @usage XPoT <- PoTselect(X, p, separation)
#' 
#' @return list containing the elements
#'   \item{pot}{selected peaks} 
#'   \item{ind}{indices in the sequence of the selected peaks in pot}  
#'   \item{p}{same as input}    
#'
#' @author Cees de Valk \email{ceesfdevalk@gmail.com}
#' 
#' @export
PoTselect <- function(X, p, separation) {
  
  n <- length(X)
  sX <- -sort(-X)
  q10 <- sX[round(n*p)]
  ind <- 1+which((X[2:(n-1)]>= X[1:(n-2)]) & (X[2:(n-1)]>= X[3:n]) & (X[2:(n-1)]> q10))
  pot <- X[ind]
  count <- 0
  repeat{
    di <- diff(ind)
    if (min(di)>= separation) {
      break
    }
    ldi <- length(di)
    imins <- which((di< separation) & (diff(pot)> 0))
    iplus <- which((di< separation) & (diff(pot)<= 0))  
    pot[imins] <- 0
    pot[iplus+1] <- 0
    id <- pot>0
    ind <- ind[id]
    pot <- X[ind]
    count <- count+1
  }
  potdata <- list("ind"= ind, "pot"= pot, "p"=p)
}



