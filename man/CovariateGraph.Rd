% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CovariateGraph.R
\name{CovariateGraph}
\alias{CovariateGraph}
\title{CovariateGraph}
\usage{
Value <- CovariateGraph(X, nbin)
}
\arguments{
\item{X}{Data of primary variate X[, 1] and covariate X[, 2] (double(n, 2))}

\item{nbin}{(optional) No. of bins: range of covariate is subdivided by nbin bins and covariate is rounded
to centre of bin (double(1))}
}
\value{
list of length nbin, containing for each bin and each value x of the primary variate 
        coinciding with a covariate value in this bin, the number of times that the primary 
        variate is larger than or equal to x
}
\description{
Conditional probability of a covariate belonging to some bin, given that the primary variate exceeds a threshold
}
\details{
If no nbin is specified, then the covariate is assumed to be discrete, and the
         bins are its range of values 
         
         For the i-th bin, the graph of its conditional probability p2 as a function of the probability 
         p1 of exceedance of the value of the primary variate is given by p2 <- (1:length(l))/l, 
         p1 <- l/n with n <- dim(X)[1] and l <- Value[[i]]
}
\author{
Cees de Valk \email{ceesfdevalk@gmail.com}
}
