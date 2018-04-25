h <- function(theta,lambda) {
  #
  # module: h.r
  #
  # purpose: calculate integral of y^theta, y from 1 to lambda
  #
  # usage:  x <- h(theta, lambda)
  #
  # theta       double(nt)
  # lambda      double(nl)  must be positive
  # x           double
  #
  # remark: If nt= nl, then x has length nt. Else, x has size (nt,nl).
  #
  # method: analytic
  #
  eps <- 5.e-16
  nt <- length(theta)
  nl <- length(lambda)
  x <- NA

  if (length(theta)== 1) {
    if (abs(theta)> eps) {
      x <- (lambda^theta-1)/theta
    } else {
      x <- log(lambda)
    }
  } else if (length(lambda)== length(theta)) {
    x <- log(lambda)
    id <- which(abs(theta)> eps)
    x[id] <- (lambda[id]^theta[id]-1)/theta[id];
  } else if (length(lambda)== 1 ) {
    x <- rep(log(lambda),nt)
    dim(x) <- dim(theta)
    id <- which(abs(theta)> eps)
    x[id] <- (lambda^theta[id]-1)/theta[id];
  } else {
    theta <- rep(theta,nl)
    lambda <- rep(lambda,nt)
    x <- h(theta, lambda)
  }
  return(x)
}
