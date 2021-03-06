% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/initialization.R
\name{tauDown_Q}
\alias{tauDown_Q}
\title{Construct initial \eqn{\tau} from \eqn{Q+1}}
\usage{
tauDown_Q(tau, n_perturb = 1)
}
\arguments{
\item{tau}{\eqn{\tau}}

\item{n_perturb}{Number of different perturbations on k-means result}
}
\value{
List of matrixes of initial values for \eqn{\tau} for \eqn{Q} groups from value obtained at \eqn{Q+1}
}
\description{
Construct initial \eqn{\tau} with \eqn{Q} groups from value obtained at \eqn{Q+1} groups
}
\examples{
# Generate first initial tau for generated_Q3 data

n <- 50
Dmax <- 2^3
Q <- 3
d_part <- 1 # less than 3 (owing to Dmax)
n_perturb <- 2
perc_perturb <- 0.2
n_random <- 1
directed <- FALSE

data <- list(Nijk = statistics(generated_Q3$data, n, Dmax, directed = FALSE))

tau <- tauInitial(data,n,Q,d_part,n_perturb,perc_perturb,n_random,directed)

tau.list <- tauDown_Q(tau[[1]],1)

}
