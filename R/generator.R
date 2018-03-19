####################################################
#### Functions to generate data under dynppsbm
####################################################



#' Poisson process
#'
#' Generate realizations of an inhomogeneous Poisson process with an intensity function
#'
#' @param intens Intensity function defined on [0,Time] (needs to be positive)
#' @param Time Final time
#' @param max.intens Upper bound of intensity on [0,Time]
#'
#' @return Vector of realizations of the PP
#'
#' @export
#'
#' @examples
#' # Generate a Poisson Process with intensity function
#' # intens= function(x) 100*x*exp(-8*x)
#' # and max.intens = 5
#'
#' intens <- function(x) 100*x*exp(-8*x)
#'
#' poissonProcess <- generatePP(intens, Time=30, max.intens=1)
#'
generatePP <- function(intens, Time, max.intens){
  area <- Time*max.intens
  M <- rpois(1,area)
  candid.points <- runif(M,0,Time)
  prob <- intens(candid.points)/max.intens
  thinning <- rbinom(rep(1,M),1,prob)
  points <- candid.points[thinning==1]
  return(points)
}



#' Data under dynppsbm
#'
#' Generate data under dynppsbm
#'
#' @param intens List containing intensity functions \eqn{\alpha^{(q,l)}} and upper bounds of intensities
#' @param Time Final time
#' @param n Total number of nodes
#' @param prop.groups Vector of group proportions (probability to belong to a group), should be of length \eqn{Q}
#' @param directed Boolean for directed (TRUE) or undirected (FALSE) case. If directed=TRUE then intens should be of length \eqn{Q^2} and if directed =FALSE then length \eqn{Q*(Q+1)/2}
#'
#' @return Simulated data, latent group variables and intensities \eqn{\alpha^{(q,l)}}
#'
#' @export
#'
#' @references
#'
#' ANDERSEN, P. K., BORGAN, Ø., GILL, R. D. & KEIDING, N. (1993). Statistical models based on counting processes. Springer Series in Statistics. Springer-Verlag, New York.
#'
#' DAUDIN, J.-J., PICARD, F. & ROBIN, S. (2008). A mixture model for random graphs. Statist. Comput. 18, 173–183.
#'
#' MATIAS, C., REBAFKA, T. & VILLERS, F. (2018).  A semiparametric extension of the stochastic block model for longitudinal networks. Biometrika.
#'
#' MATIAS, C. & ROBIN, S. (2014). Modeling heterogeneity in random graphs through latent space models: a selective review. Esaim Proc. & Surveys 47, 55–74.
#'
#' @examples
#' # Generate data from an undirected graph with n=10 individuals and Q=2 clusters
#'
#' # equal cluster proportions
#' prop.groups <- c(0.5,0.5)
#'
#' # 3 different intensity functions :
#' intens <- list(NULL)
#' intens[[1]] <- list(intens= function(x) 100*x*exp(-8*x),max=5)
#'     # (q,l) = (1,1)
#' intens[[2]] <- list(intens= function(x) exp(3*x)*(sin(6*pi*x-pi/2)+1)/2,max=13)
#'     # (q,l) = (1,2)
#' intens[[3]] <- list(intens= function(x) 8.1*(exp(-6*abs(x-1/2))-.049),max=8)
#'     # (q,l) = (2,2)
#'
#' # generate data :
#' obs <- generateDynppsbm(intens,Time=1,n=10,prop.groups,directed=FALSE)
#'
#' # latent variables (true clustering of the individuals)
#' obs$z
#'
#' # number of time events :
#' length(obs$data$time.seq)
#'
#' # number of interactions between each pair of individuals:
#' table(obs$data$type.seq)
#'
generateDynppsbm <- function(intens,Time,n,prop.groups,directed=TRUE){
  Q <- length(prop.groups)

  # test coherence
  if (directed){
    if (length(intens)!=Q^2){ stop("not a correct number of intensities")}}
  else{
    if (length(intens)!=Q*(Q+1)/2){stop("not a correct number of intensities")}}

  if (directed){
    N <- n*(n-1)}
  else{
    N <- n*(n-1)/2}

  # Draw hidden groups
  z <- rmultinom(n, 1, prop.groups)
    group <- colMaxs(z,FALSE) # apply(z,2,which.max)

  # Generate the process
  time.seq <- NULL
  type.seq <- NULL
  if (directed){ # directed
    for (i in 1:n){
      for (j in (1:n)[-i]){
        type.ql <- convertGroupPair(group[i],group[j],Q,directed)
        proc <- generatePP(intens[[type.ql]]$intens,Time,intens[[type.ql]]$max)
        if (length(proc)>0){
          time.seq <- c(time.seq,proc)  	# not ordered
          type.seq <- c(type.seq,rep(convertNodePair(i,j,n,directed),length(proc)))
        }
      }
    }
  }
  else{ # undirected
    for (i in 1:(n-1)){
      for (j in (i+1):n){
        type.ql <- convertGroupPair(group[i],group[j],Q,directed)
        proc <- generatePP(intens[[type.ql]]$intens,Time,intens[[type.ql]]$max)
        if (length(proc)>0){
          time.seq <- c(time.seq,proc)  	# not ordered
          type.seq <- c(type.seq,rep(convertNodePair(i,j,n,directed),length(proc)))
        }
      }
    }
  }
  ordre <- order(time.seq)
  time.seq <- time.seq[ordre]
  type.seq <- type.seq[ordre]
  data <-	list(time.seq=time.seq,type.seq=type.seq,Time=Time)

  return(list(data=data,z=z,intens=intens))
}

####################################################
####  PIECEWISE CONSTANT
####################################################



#' Poisson process  with piecewise constant intensities
#'
#' Generate realizations of a Poisson process with piecewise constant intensities
#'
#' @param intens Vector with the constants of the intensities (defined on a regular partition of interval [0,Time])
#' @param Time Time
#'
#' @export
#'
#' @examples
#' intens <- c(1,3,8)
#' constpp <- generatePPConst(intens, 10)
#'
generatePPConst <- function(intens,Time){
  K <- length(intens)
  long.subint <- Time/K
  area <- intens*long.subint
  nb.points.interv <- rpois(rep(1,K),area)
  points <- NULL
  for (k in 1:K)
    if (nb.points.interv[k]>0)
      points <- c(points,runif(nb.points.interv[k],(k-1)*long.subint,k*long.subint))
  return(points)
}


#' Data under dynppsbm with piecewise constant intensities
#'
#' Generate data under dynppsbm with piecewise constant intensities
#'
#' @param intens Matrix with piecewise constant intensities \eqn{\alpha^{(q,l)}} (each row gives the constants of the piecewise constant intensity for a group pair \eqn{(q,l)})
#' @param Time Time
#' @param n Total number of nodes
#' @param prop.groups Vector of group proportions, should be of length \eqn{Q}
#' @param directed Boolean for directed (TRUE) or undirected (FALSE) case
#'
#' If directed then intens should be of length \eqn{Q^2} and if undirected then length \eqn{Q*(Q+1)/2}
#'
#' @export
#'
#' @examples
#' intens1 <- c(1,3,8)
#' intens2 <- c(2,3,6)
#'
#' intens <- matrix(c(intens1,intens2,intens1,intens2),4,3)
#'
#' Time <- 10
#' n <- 20
#' prop.groups <- c(0.2,0.3)
#' dynppsbm <- generateDynppsbmConst(intens,Time,n,prop.groups,directed=TRUE)
#'
generateDynppsbmConst <- function(intens,Time,n,prop.groups,directed=TRUE){
  Q <- length(prop.groups)

  # test coherence
  N_Q <- if (directed)  Q^2 else Q*(Q+1)/2
  if (nrow(intens)!=N_Q) stop("not a correct number of intensities")

  # Draw hidden groups
  z <- rmultinom(n, 1, prop.groups)
  group <- colMaxs(z, FALSE) # apply(z, 2, which.max)

  # Generate the process
  time.seq <- NULL
  type.seq <- NULL
  vec_i <- if (directed)  1:n else 1:(n-1)
  for (i in vec_i){
    vec_j <- if (directed) (1:n)[-i] else (i+1):n
    for (j in vec_j){
      type.ql <- convertGroupPair(group[i],group[j],Q,directed)
      proc <- generatePPConst(intens[type.ql,],Time)
      if (length(proc)>0){
        time.seq <- c(time.seq,proc)  	# not ordered
        type.seq <- c(type.seq,rep(convertNodePair(i,j,n,directed),length(proc)))
      }
    }
  }
  ordre <- order(time.seq)
  time.seq <- time.seq[ordre]
  type.seq <- type.seq[ordre]
  data <-	list(time.seq=time.seq,type.seq=type.seq,Time=Time)

  return(list(data=data,z=z,intens=intens))
}
