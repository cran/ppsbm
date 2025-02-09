% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{generated_Q3}
\alias{generated_Q3}
\title{Example dataset}
\format{
A list of 3 components:
\describe{
  \item{\code{data}}{Observed data is itself a list of 3 components:
   \itemize{
     \item \code{time.seq} - Vector containing the times (in [0,1]) of the events (length M).
     \item \code{type.seq} - Vector containing the types in \eqn{\{1,\dots, N\}} of the events (length M). Here, \eqn{N=n(n-1)/2} (undirected). 
     \item \code{Time} - Positive real number. [0,Time] is the total time interval of observation.
   }
  }
  \item{\code{z}}{Latent variables. A matrix with size \eqn{Q\times n} and entries 1 (cluster q contains node i) or 0 (else).}
  \item{\code{intens}}{Intensities used to simulate data. A list of \eqn{Q(Q+1)/2} intensity functions. Each one is given as a list of 2 components:
   \itemize{
   \item \code{intens} - a positive function. The intensity function  \eqn{\alpha^{(q,l)}}
   \item \code{max} - positive real number. An upper bound on function \eqn{\alpha^{(q,l)}}
   }
  }
}
}
\usage{
generated_Q3
}
\description{
Example of undirected dataset with \eqn{n=50} individuals in \eqn{Q=3} clusters and final observation Time=1.
}
\details{
This random datatset was obtained using the following code
\preformatted{

intens <- list(NULL)
intens[[1]] <- list(intens=function(x) return (rep(4,length(x))), max=4.1)
intens[[2]] <- list(intens=function(x){
   y <- rep(0,length(x))
   y[x<.25] <- 4
   y[x>.75] <- 10
   return(y)
 }, max=10.1)
intens[[3]] <- list(intens=function(x) return(8*(1/2-abs(x-1/2))), max=4.1)
intens[[4]] <- list(intens=function(x) return(100*x*exp(-8*x)), max=4.698493)
intens[[5]] <- list(intens=function(x) return(exp(3*x)*(sin(6*pi*x-pi/2)+1)/2), max=12.59369)
intens[[6]] <- list(intens=function(x) return(8.1*(exp(-6*abs(x-1/2))-.049)), max=7.8031)

generated_Q3 <- generateDynppsbm(intens,Time=1,n=50,prop.groups=rep(1/3,3),directed=F)
}
}
\references{
MATIAS, C., REBAFKA, T. & VILLERS, F. (2018).  A semiparametric extension of the stochastic block model for longitudinal networks. Biometrika. 105(3): 665-680.
}
\keyword{datasets}
