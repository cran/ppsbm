% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/output.R
\name{bootstrap_and_CI}
\alias{bootstrap_and_CI}
\title{Bootstrap and Confidence Bands}
\usage{
bootstrap_and_CI(
  sol,
  Time,
  R,
  alpha = 0.05,
  nbcores = 1,
  d_part = 5,
  n_perturb = 10,
  perc_perturb = 0.2,
  directed,
  filename = NULL
)
}
\arguments{
\item{sol}{One list (for one value of \eqn{Q}) output by \link[ppsbm]{mainVEM} with \code{hist} method.}

\item{Time}{Positive real number. [0,Time] is the total time interval of observation.}

\item{R}{Number of bootstrap samples.}

\item{alpha}{Level of confidence: \eqn{1- \alpha}.}

\item{nbcores}{Number of cores for parallel execution.

If set to 1 it does sequential execution.

Beware: parallelization with fork (multicore): doesn't work on Windows!}

\item{d_part}{Maximal level for finest partitions of time interval [0,T], used for kmeans initializations on the bootstrap samples
\itemize{
  \item Algorithm takes partition up to depth \eqn{2^d} with \eqn{d=1,...,d_{part}}
  \item Explore partitions \eqn{[0,T], [0,T/2], [T/2,T], ... [0,T/2^d], ...[(2^d-1)T/2^d,T]}
  \item Total number of partitions \eqn{npart= 2^{(d_{part} +1)} - 1}
}}

\item{n_perturb}{Number of different perturbations on k-means result on the bootstrap samples.}

\item{perc_perturb}{Percentage of labels that are to be perturbed (= randomly switched)  on the bootstrap samples.}

\item{directed}{Boolean for directed (TRUE) or undirected (FALSE) case.}

\item{filename}{Name of the file where to save the results.}
}
\description{
Plots confidence bands for estimated intensities between pairs of groups obtained by bootstrap.
}
\details{
Not for sparse models and only for histogram method.
}
\examples{

# data of a synthetic graph with 50 individuals and 3 clusters

n <- 50
Q <- 3

Time <- generated_Q3$data$Time
data <- generated_Q3$data
z <- generated_Q3$z

K <- 2^3

# VEM-algo hist:
sol.hist <- mainVEM(list(Nijk=statistics(data,n,K,directed=FALSE),Time=Time),
n,Qmin=3,directed=FALSE,method='hist',d_part=1,n_perturb=0)[[1]]

# compute bootstrap confidence bands
boot <- bootstrap_and_CI(sol.hist,Time,R=10,alpha=0.1,nbcores=1,d_part=1,n_perturb=0,
     directed=FALSE)

# plot confidence bands
alpha.hat <- exp(sol.hist$logintensities.ql)
vec.x <- (0:K)*Time/K
ind.ql <- 0
par(mfrow=c(2,3))
for (q in 1:Q){
  for (l in q:Q){
    ind.ql <- ind.ql+1
    ymax <- max(c(boot$CI.limits[ind.ql,2,],alpha.hat[ind.ql,]))
    plot(vec.x,c(alpha.hat[ind.ql,],alpha.hat[ind.ql,K]),type='s',col='black',
        ylab='Intensity',xaxt='n',xlab= paste('(',q,',',l,')',sep=""),
        cex.axis=1.5,cex.lab=1.5,ylim=c(0,ymax),main='Confidence bands')
    lines(vec.x,c(boot$CI.limits[ind.ql,1,],boot$CI.limits[ind.ql,1,K]),col='blue',
        type='s',lty=3)
    lines(vec.x,c(boot$CI.limits[ind.ql,2,],boot$CI.limits[ind.ql,2,K]),col='blue',
        type='s',lty=3)
  }
}

}
