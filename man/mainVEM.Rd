% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vem_dynppsbm.R
\name{mainVEM}
\alias{mainVEM}
\title{Adaptive VEM algorithm}
\usage{
mainVEM(
  data,
  n,
  Qmin,
  Qmax = Qmin,
  directed = TRUE,
  sparse = FALSE,
  method = c("hist", "kernel"),
  init.tau = NULL,
  cores = 1,
  d_part = 5,
  n_perturb = 10,
  perc_perturb = 0.2,
  n_random = 0,
  nb.iter = 50,
  fix.iter = 10,
  epsilon = 1e-06,
  filename = NULL
)
}
\arguments{
\item{data}{Data format depends on the estimation method used!!
\enumerate{
  \item Data with \strong{hist} method - List with 2 components:
    \describe{
      \item{data$Time}{Positive real number. [0,data$Time] is the total time interval of observation}

      \item{data$Nijk}{Data matrix with counts per process \eqn{N_{ij}} and sub-intervals ; matrix of size \eqn{N*K}  where \eqn{N = n(n-1)} or \eqn{n(n-1)/2} is the number of possible node pairs in the graph and \eqn{K = 2^{dmax}} is the size of the finest partition in the  histogram approach.

      Counts are pre-computed - Obtained through function \link[ppsbm]{statistics} on data with second format and using a number of subintervals K as a power of 2.}
    }

  \item Data with \strong{kernel} method - List with 3 components:
    \describe{
      \item{data$time.seq}{Vector of the time points of the events (size M).}

      \item{data$type.seq}{Vector of the corresponding node pair indexes in format output by \link[ppsbm]{convertNodePair} of the events (same size M).}

    \item{data$Time}{[0,data$Time] is the total time interval of observation}
    }
}}

\item{n}{Total number of nodes,  \eqn{1\le i \le n}}

\item{Qmin}{Minimum number of groups}

\item{Qmax}{Maximum number of groups}

\item{directed}{Boolean for directed (TRUE) or undirected (FALSE) case}

\item{sparse}{Boolean for sparse (TRUE) or not sparse (FALSE) case}

\item{method}{Either \code{hist} for histogram method or \code{kernel} for kernel method

Beware: \code{hist} is recommended (much faster). You can obtain smooth estimated intensities by using \link[ppsbm]{kernelIntensities} on the output of the \code{hist} method.}

\item{init.tau}{List of initial values of \eqn{\tau} - all tau's are matrices with size \eqn{Q\times n} (might be with different values of Q)}

\item{cores}{Number of cores for parallel execution

If set to 1 it does sequential execution

Beware: parallelization with fork (multicore) : doesn't work on Windows!}

\item{d_part}{Maximal level for finest partition of time interval [0,T] used for k-means initializations.
\itemize{
  \item Algorithm takes partition up to depth \eqn{2^d} with \eqn{d=1,...,d_{part}}
  \item Explore partitions \eqn{[0,T], [0,T/2], [T/2,T], ... [0,T/2^d], ...[(2^d-1)T/2^d,T]}
  \item Total number of partitions \eqn{npart= 2^{(d_{part} +1)} -1}
}}

\item{n_perturb}{Number of different perturbations on k-means result

      When \eqn{Qmin < Qmax}, number of perturbations on the result with \eqn{Q-1} or \eqn{Q+1} groups}

\item{perc_perturb}{Percentage of labels that are to be perturbed (= randomly switched)}

\item{n_random}{Number of completely random initial points. The total number of initializations for the VEM is \eqn{npart*(1+n_{perturb}) +n_{random}}}

\item{nb.iter}{Number of iterations of the VEM algorithm}

\item{fix.iter}{Maximum number of iterations of the fixed point into the VE step}

\item{epsilon}{Threshold for the stopping criterion of VEM and fixed point iterations}

\item{filename}{Name of the file where to save the results along the computation (increasing steps for \eqn{Q}, these are the longest).

   The file will contain a list of 'best' results.}
}
\value{
The function outputs a list of Qmax-Qmin+1 components. Each component is the solution obtained for a number of clusters Q, with \eqn{Qmin\le Q \le Qmax} and is a list of 8 elements:
 \itemize{
   \item {\code{tau}} - Matrix with size \eqn{Q\times n} containing the estimated values in \eqn{(0,1)} that cluster q contains node i.
   \item {\code{rho}} - When method=\code{hist} only. Either 1 (non sparse method) or a vector with length \eqn{Q(Q+1)/2} (undirected case) or \eqn{Q^2} (directed case) with estimated values for the sparsity parameters \eqn{\rho^{(q,l)}}. See Section S6 in the supplementary material paper of Matias et al. (Biometrika, 2018) for more details.
   \item {\code{beta}} - When method=\code{hist} only. Vector with length \eqn{Q(Q+1)/2} (undirected case) or \eqn{Q^2} (directed case) with estimated values for the sparsity parameters \eqn{\beta^{(q,l)}}. See Section S6 in the supplementary material paper Matias et al. (Biometrika, 2018) for more details.
   \item {\code{logintensities.ql}} - When method=\code{hist} only. Matrix with size \eqn{Q(Q+1)/2\times K} (undirected case) or \eqn{Q^2\times K} (directed case). Each row contains estimated values of the log intensity function \eqn{\log(\alpha^{(q,l)})} on a regular partition with K parts of the time interval [0,Time].
   \item {\code{best.d}} - When method=\code{hist} only. Vector with length \eqn{Q(Q+1)/2} (undirected case) or \eqn{Q^2} (directed case) with estimated value for the exponent of the best partition to estimate intensity \eqn{\alpha^{(q,l)}}. The best number of parts is \eqn{K=2^d}.
   \item {\code{J}} - Estimated value of the ELBO.
   \item {\code{run}} - Which run of the algorithm gave the best solution. A run relies on a specific initialization of the algorithm. A negative value maybe obtained in the decreasing phase (for Q) of the algorithm.
   \item {\code{converged}} - Boolean. If TRUE, the algorithm stopped at convergence. Otherwise it stopped at the maximal number of iterations.
 }
}
\description{
Principal adaptive VEM algorithm for histogram (with model selection) or for kernel method.
}
\details{
The sparse version works only for the histogram approach.
}
\examples{
# load data of a synthetic graph with 50 individuals and 3 clusters
n <- 20
Q <- 3

Time <- generated_Q3_n20$data$Time
data <- generated_Q3_n20$data
z <- generated_Q3_n20$z

step <- .001
x0 <- seq(0,Time,by=step)
intens <-  generated_Q3_n20$intens

# VEM-algo kernel
sol.kernel <- mainVEM(data,n,Q,directed=FALSE,method='kernel', d_part=0,
    n_perturb=0)[[1]]
# compute smooth intensity estimators
sol.kernel.intensities <- kernelIntensities(data,sol.kernel$tau,Q,n,directed=FALSE)
# eliminate label switching
intensities.kernel <- sortIntensities(sol.kernel.intensities,z,sol.kernel$tau,
    directed=FALSE)

# VEM-algo hist
# compute data matrix with precision d_max=3 (ie nb of parts K=2^{d_max}=8).
K <- 2^3
Nijk <- statistics(data,n,K,directed=FALSE)
sol.hist <- mainVEM(list(Nijk=Nijk,Time=Time),n,Q,directed=FALSE, method='hist',
    d_part=0,n_perturb=0,n_random=0)[[1]]
log.intensities.hist <- sortIntensities(sol.hist$logintensities.ql,z,sol.hist$tau,
     directed=FALSE)

# plot estimators
par(mfrow=c(2,3))
ind.ql <- 0
for (q in 1:Q){
  for (l in q:Q){
    ind.ql <- ind.ql + 1
    true.val <- intens[[ind.ql]]$intens(x0)
    values <- c(intensities.kernel[ind.ql,],exp(log.intensities.hist[ind.ql,]),true.val)
    plot(x0,true.val,type='l',xlab=paste0("(q,l)=(",q,",",l,")"),ylab='',
        ylim=c(0,max(values)+.1))
    lines(seq(0,1,by=1/K),c(exp(log.intensities.hist[ind.ql,]),
        exp(log.intensities.hist[ind.ql,K])),type='s',col=2,lty=2)
    lines(seq(0,1,by=.001),intensities.kernel[ind.ql,],col=4,lty=3)
  }
}

}
\references{
DAUDIN, J.-J., PICARD, F. & ROBIN, S. (2008). A mixture model for random graphs. Statist. Comput. 18, 173–183.

DEMPSTER, A. P., LAIRD, N. M. & RUBIN, D. B. (1977). Maximum likelihood from incomplete data via the EM algorithm. J. Roy. Statist. Soc. Ser. B 39, 1–38.

JORDAN, M., GHAHRAMANI, Z., JAAKKOLA, T. & SAUL, L. (1999). An introduction to variational methods for graphical models. Mach. Learn. 37, 183–233.

MATIAS, C., REBAFKA, T. & VILLERS, F. (2018).  A semiparametric extension of the stochastic block model for longitudinal networks. Biometrika. 105(3): 665-680.

MATIAS, C. & ROBIN, S. (2014). Modeling heterogeneity in random graphs through latent space models: a selective review. Esaim Proc. & Surveys 47, 55–74.
}
