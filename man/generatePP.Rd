% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generator.R
\name{generatePP}
\alias{generatePP}
\title{Poisson process generator}
\usage{
generatePP(intens, Time, max.intens)
}
\arguments{
\item{intens}{Intensity function defined on [0,Time] (needs to be positive).}

\item{Time}{Positive real number. [0,Time] is the total time interval of observation.}

\item{max.intens}{Positive real number. Upper bound of intensity on [0,Time].}
}
\value{
Vector with the values of one realization of the PP.
}
\description{
Generates one realization of an inhomogeneous Poisson process (PP) with given intensity function (using the thinning method).
}
\examples{
# Generate a Poisson Process on time interval [0,30] with intensity function
# intens= function(x) 100*x*exp(-8*x)
# using max.intens = 5

intens <- function(x) 100*x*exp(-8*x)

generatePP(intens, Time=30, max.intens=5)

}
