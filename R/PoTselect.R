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
#' @author Alexandre Tribut, Cees de Valk \email{ceesfdevalk@gmail.com}
#' 
#' 
#' @export
PoTselect <- function(data, p=NULL, separation=1, threshold=NULL) {
  if (any(is.na(data))) {
    stop("NA's in data.")
  }
  n <- length(data)
  sdata <- -sort(-data)
  
  if(!is.null(p)){
    threshold <- sdata[round(n * p)]
    # print(threshold)
  } 
  df <- diff(c(0, data > threshold, 0))
  
  # Find the start and end indices of each excursion
  excursion_starts <- which(df== 1)
  excursion_ends <- which(df== -1)-1
  
  # Ensure they match up
  if (length(excursion_starts) != length(excursion_ends)) {
    stop("Mismatch in excursion start and end points")
  }
  
  
  # test_time = now()
  # Find the start and end indices of each excursion
  ind <- sapply(1:length(excursion_starts), function(i) {
    excursion_starts[i]+which.max(data[excursion_starts[i]:excursion_ends[i]])-1
  })
  
  # Initial potential excursions above threshold
  pot <- data[ind]
  count <- 0
  
  # Remove excursions that are too close
  repeat {
    di <- diff(ind)
    if (min(di) >= separation) {
      break
    }
    
    ldi <- length(di)
    imins <- which((di < separation) & (diff(pot) > 0))
    iplus <- which((di < separation) & (diff(pot) <= 0))
    
    pot[imins] <- 0
    pot[iplus + 1] <- 0
    
    id <- pot > 0
    ind <- ind[id]
    pot <- data[ind]
    
    count <- count + 1
  }
  
  potdata <- list("ind" = ind, "pot" = pot, "p" = p)
  
  return(potdata)
}



find_max_window <- function(data, pot_ind, window_size) {
  best_window <- NULL
  counts <- rep(0, length(data))
  for (i in seq(1,length(data),max(1,round(window_size/10)))) {
    cat("\r",i,'/',length(data))
    flush.console()
    start <- max(1, i - window_size + 1)
    end <- min(length(data), i)
    counts[i] <- sum(pot_ind>start & pot_ind < end)
    
    max_count <- max(counts)
    best_window <- c(which.max(counts),which.max(counts)+window_size)
  }
  
  
  return(list(best_window = best_window, max_count = max_count))
}