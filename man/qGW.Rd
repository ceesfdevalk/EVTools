% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qGW.R
\name{qGW}
\alias{qGW}
\title{qGW}
\usage{
Value <- qGW(p, p0, q0, scale, theta, q0Std, scaleStd, thetaStd)
}
\arguments{
\item{p}{probability(ies) of exceedance at which to evaluate the quantile(s) (double(np))}

\item{p0}{reference probability for which the quantile q0 is given (double(1))}

\item{q0}{see under p0}

\item{scale}{GW scale parameter at p0 (double(1))}

\item{theta}{GW tail index (double(1))}

\item{q0Std}{(optional) standard deviation of q0 (double(1))}

\item{scaleStd}{(optional) standard deviation of scale (double(1))}

\item{thetaStd}{(optional) standard deviation of theta (double(1))}
}
\value{
list containing quantile(s) and optionally, their standard deviations
}
\description{
quantile of the GW (Generalised Weibull) tail and (optionally) its standard deviation
}
\references{
De Valk, C. and Cai, J.J. (2018), A high quantile estimator based on 
the log-generalized Weibull tail limit. Econometrics and Statistics 6, 107-128, see
\url{https://doi.org/10.1016/j.ecosta.2017.03.001}
}
\author{
Cees de Valk \email{ceesfdevalk@gmail.com}
}
