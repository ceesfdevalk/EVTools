% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EI.R
\name{EI}
\alias{EI}
\title{EI}
\usage{
Value <- EI(X, l, ngr, makeplot)
}
\arguments{
\item{X}{time series (double(nX))}

\item{l}{(optional, defailt= 50) no of data points above highest level (integer(1))}

\item{ngr}{(optional, default= 10) number of levels (integer(1))}

\item{makeplot}{(optional, default= FALSE) (logical(1))}
}
\value{
A list, with members: 
  \item{k}{no. of data points above level, for each of ngr levels}    
  \item{p}{sample fraction above level, for each of ngr levels}    
  \item{EIupcross}{estimate of EI based on upcrossings} 
  \item{EIFS}{estimate of EI from Ferro & Segers (2003)}
}
\description{
Computes estimates of the Extremal Index
}
\references{
M. R. Leadbetter, M.R., Lindgren, G. and Rootzen, H. (1983), 
 Extremes and Related Properties of Random Sequences and Processes. 
 Springer.
 Ferro, A.T. and J. Segers (2003), Inference for clusters 
 of extreme values. J. R. Statist. Soc. B (2003) 65, 545–556.
}
\author{
Cees de Valk \email{ceesfdevalk@gmail.com}
}
