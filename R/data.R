#' Example dataset 
#' 
#' Example of undirected dataset with \eqn{n=50} individuals in \eqn{Q=3} clusters and final observation Time=1.
#'
#' @format A list of 3 components:
#' \describe{
#'   \item{\code{data}}{Observed data is itself a list of 3 components:
#'    \itemize{
#'      \item \code{time.seq} - Vector containing the times (in [0,1]) of the events (length M).
#'      \item \code{type.seq} - Vector containing the types in \eqn{\{1,\dots, N\}} of the events (length M). Here, \eqn{N=n(n-1)/2} (undirected). 
#'      \item \code{Time} - Positive real number. [0,Time] is the total time interval of observation.
#'    }
#'   }
#'   \item{\code{z}}{Latent variables. A matrix with size \eqn{Q\times n} and entries 1 (cluster q contains node i) or 0 (else).}
#'   \item{\code{intens}}{Intensities used to simulate data. A list of \eqn{Q(Q+1)/2} intensity functions. Each one is given as a list of 2 components:
#'    \itemize{
#'    \item \code{intens} - a positive function. The intensity function  \eqn{\alpha^{(q,l)}}
#'    \item \code{max} - positive real number. An upper bound on function \eqn{\alpha^{(q,l)}}
#'    }
#'   }
#' }
#' 
#' @details
#' This random datatset was obtained using the following code
#' \preformatted{
#' 
#' intens <- list(NULL)
#' intens[[1]] <- list(intens=function(x) return (rep(4,length(x))), max=4.1)
#' intens[[2]] <- list(intens=function(x){
#'    y <- rep(0,length(x))
#'    y[x<.25] <- 4
#'    y[x>.75] <- 10
#'    return(y)
#'  }, max=10.1)
#' intens[[3]] <- list(intens=function(x) return(8*(1/2-abs(x-1/2))), max=4.1)
#' intens[[4]] <- list(intens=function(x) return(100*x*exp(-8*x)), max=4.698493)
#' intens[[5]] <- list(intens=function(x) return(exp(3*x)*(sin(6*pi*x-pi/2)+1)/2), max=12.59369)
#' intens[[6]] <- list(intens=function(x) return(8.1*(exp(-6*abs(x-1/2))-.049)), max=7.8031)
#' 
#' generated_Q3 <- generateDynppsbm(intens,Time=1,n=50,prop.groups=rep(1/3,3),directed=F)
#' }
#' 
#' 
#' @references
#' MATIAS, C., REBAFKA, T. & VILLERS, F. (2018).  A semiparametric extension of the stochastic block model for longitudinal networks. Biometrika. 105(3): 665-680. 
#' 
"generated_Q3"


#' Example dataset 
#' 
#' Example of undirected dataset with \eqn{n=20} individuals in \eqn{Q=3} clusters and observation Time=1
#'
#' @format A list of 3 components:
#' \describe{
#'   \item{\code{data}}{Observed data is itself a list of 3 components:
#'    \itemize{
#'      \item \code{time.seq} - Vector containing the times of the events (length M)
#'      \item \code{type.seq} - Vector containing the types of the events (length M) 
#'      \item \code{Time} - Positive real number. [0,Time] is the total time interval of observation
#'    }
#'   }
#'   \item{\code{z}}{Latent variables. A matrix with size \eqn{Q\times n} and entries 1 (cluster q contains node i) or 0 (else).}
#'   \item{\code{intens}}{Intensities used to simulate data. A list of \eqn{Q(Q+1)/2} intensity functions. Each one is given as a list of 2 components:
#'    \itemize{
#'    \item \code{intens} - a positive function. The intensity function  \eqn{\alpha^{(q,l)}}
#'    \item \code{max} - positive real number. An upper bound on function \eqn{\alpha^{(q,l)}}
#'    }
#'   }
#' }
#' 
#' @details
#' This random datatset was obtained using the following code
#' \preformatted{
#' intens <- list(NULL)
#' intens[[1]] <- list(intens=function(x) return (rep(4,length(x))), max=4.1)
#' intens[[2]] <- list(intens=function(x){
#'    y <- rep(0,length(x))
#'    y[x<.25] <- 4
#'    y[x>.75] <- 10
#'    return(y)
#'  }, max=10.1)
#' intens[[3]] <- list(intens=function(x) return(8*(1/2-abs(x-1/2))), max=4.1)
#' intens[[4]] <- list(intens=function(x) return(100*x*exp(-8*x)), max=4.698493)
#' intens[[5]] <- list(intens=function(x) return(exp(3*x)*(sin(6*pi*x-pi/2)+1)/2), max=12.59369)
#' intens[[6]] <- list(intens=function(x) return(8.1*(exp(-6*abs(x-1/2))-.049)), max=7.8031)
#' 
#' generated_Q3_n20 <- generateDynppsbm(intens,Time=1,n=20,prop.groups=rep(1/3,3),directed=F)
#' }
#' 
#' @references
#' MATIAS, C., REBAFKA, T. & VILLERS, F. (2018).  A semiparametric extension of the stochastic block model for longitudinal networks. Biometrika. 105(3): 665-680. 
#' 
"generated_Q3_n20"



#' Output example of \link[ppsbm]{mainVEM}
#' 
#' Output of \link[ppsbm]{mainVEM} obtained on dataset \code{generated_Q3} with \code{hist} method and Qmin=1, Qmax=5.
#'
#' @format 
#' List of 5 components. 
#' Each one is the output of the algorithm with a different value of the number of clusters \eqn{Q} for \eqn{1\le Q \le 5} and given as a list of 8 components:  
#' \describe{
#'  \item{\code{tau}}{Matrix with size \eqn{Q\times n} containing the estimated probability in \eqn{(0,1)} that cluster \eqn{q} contains node \eqn{i}.}
#'  \item{\code{rho}}{Sparsity parameter - 1 in this case (non sparse method).}
#'  \item{\code{beta}}{Sparsity parameter - 1 in this case (non sparse method).}   
#'  \item{\code{logintensities.ql}}{Matrix with size \eqn{Q(Q+1)/2\times K}. Each row contains estimated values of the log of the intensity function \eqn{\log(\alpha^{(q,l)})} on a regular partition (in \eqn{K} parts) of the time interval [0,Time].}  
#'  \item{\code{best.d}}{Vector with length \eqn{Q(Q+1)/2} (undirected case) with estimated value for the exponent of the best partition to estimate intensity \eqn{\alpha^{(q,l)}}. The best number of parts is \eqn{K=2^d}.}
#'  \item{\code{J}}{Estimated value of the ELBO}
#'  \item{\code{run}}{Which run of the algorithm gave the best solution. A run relies on a specific initialization of the algorithm. A negative value maybe obtained in the decreasing phase (for Q) of the algorithm.}  
#'  \item{\code{converged}}{Boolean. If TRUE, the algorithm stopped at convergence. Otherwise it stopped at the maximal number of iterations.} 
#' }
#'            
#' @references
#' MATIAS, C., REBAFKA, T. & VILLERS, F. (2018).  A semiparametric extension of the stochastic block model for longitudinal networks. Biometrika. 105(3): 665-680. 
#' 
#' @details
#' This solution was (randomly) obtained using the following code
#' \preformatted{
#' Nijk <- statistics(generated_Q3$data,n=50,K=8,directed=FALSE)
#' generated_sol_hist <- mainVEM(list(Nijk=Nijk,Time=1),n=50,Qmin=1,Qmax=5,directed=FALSE,method='hist')
#' }
#' 
"generated_sol_hist"



#' Output example of \link[ppsbm]{mainVEM} 
#' 
#' Output of \link[ppsbm]{mainVEM} obtained on dataset \code{generated_Q3} with \code{kernel} method and Qmin=Qmax=5.
#'
#' @format Solution for Q=5 clusters, containing 5 components:
#' \describe{
#'  \item{\code{tau}}{Matrix with size \eqn{Q\times n} containing the estimated probability in \eqn{(0,1)} that cluster \eqn{q} contains node \eqn{i}.}
#'  \item{\code{J}}{Estimated value of the ELBO}
#'  \item{\code{run}}{Which run of the algorithm gave the best solution. A run relies on a specific initialization of the algorithm. A negative value maybe obtained in the decreasing phase (for Q) of the algorithm.}  
#'  \item{\code{converged}}{Boolean. If TRUE, the algorithm stopped at convergence. Otherwise it stopped at the maximal number of iterations.} 
#' }
#' 
#' @references
#' MATIAS, C., REBAFKA, T. & VILLERS, F. (2018).  A semiparametric extension of the stochastic block model for longitudinal networks. Biometrika. 105(3): 665-680. 
#' 
#' @details
#' This solution was (randomly) obtained using the following code
#' \preformatted{
#' # WARNING - This is very long
#' generated_sol_kernel <- mainVEM(generated_Q3$data,n=50,Qmin=5,directed=FALSE,method='kernel')[[1]]
#' }
#' 
#' 
"generated_sol_kernel"

