% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/h.R
\name{h}
\alias{h}
\title{h}
\usage{
Value <- h(theta, lambda)
}
\arguments{
\item{theta}{the exponent(s) (double(nt))}

\item{lambda}{the positive number(s) at which the integral is evaluated (double(nl))}
}
\value{
real number(s): if nt= nl or nt= 1 or nl= 1, then Value has length max(nt, nl). 
Otherwise, Value is an nlxnt matrix.
}
\description{
h(theta, lambda) is the integral of x^theta for x from 1 to lambda.
}
\author{
Cees de Valk \email{ceesfdevalk@gmail.com}
}
