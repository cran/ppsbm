% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generator.R
\name{generateDynppsbmConst}
\alias{generateDynppsbmConst}
\title{Data under dynppsbm with piecewise constant intensities}
\usage{
generateDynppsbmConst(intens, Time, n, prop.groups, directed = TRUE)
}
\arguments{
\item{intens}{Matrix with piecewise constant intensities \eqn{\alpha^{(q,l)}}. Each row gives the constant values of the piecewise constant intensity for a group pair \eqn{(q,l)} on a regular partition of the time interval [0,Time].}

\item{Time}{Positive real number. [0,Time] is the total time interval of observation.}

\item{n}{Total number of nodes,  \eqn{1\le i \le n}.}

\item{prop.groups}{Vector of group proportions, should be of length \eqn{Q}.}

\item{directed}{Boolean for directed (TRUE) or undirected (FALSE) case.

If directed then \code{intens} should be of length \eqn{Q^2}, else of length \eqn{Q*(Q+1)/2}.}
}
\description{
Generate data under the Dynamic Poisson Process Stochastic Blockmodel (dynppsbm) with piecewise constant intensity functions.
}
\examples{
# Define 2 different piecewise constant intensity functions
# on a 3 parts regular partition of time interval [0,Time]
intens1 <- c(1,3,8)
intens2 <- c(2,3,6)

intens <- matrix(c(intens1,intens2,intens1,intens2),4,3)

Time <- 10
n <- 20
prop.groups <- c(0.2,0.8)
obs <- generateDynppsbmConst(intens,Time,n,prop.groups,directed=TRUE)

}
\references{
MATIAS, C., REBAFKA, T. & VILLERS, F. (2018).  A semiparametric extension of the stochastic block model for longitudinal networks. Biometrika. 105(3): 665-680.
}
