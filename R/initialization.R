####################################
#### INITIALIZATION OF VEM
####################################


####################################
#### INITIAL TAU
####################################

#' List of initial values for \eqn{\tau}
#'
#' Same function whatever directed or undirected case
#'
#' The (maximal) total number of initializations is \eqn{d_{part}*(1+n_{perturb}) + n_{random}}
#'
#' @param data Data : only needs the \eqn{N_{ijk}} field of data
#' @param n Total number of nodes
#' @param Q Total number of groups
#' @param d_part Maximal level for finest partitions of time interval [0,T], used for kmeans initializations.
#'   \itemize{
#'     \item Algorithm takes partition up to depth \eqn{2^d} with \eqn{d=1,...,d_{part}}
#'     \item Explore partitions \eqn{[0,T], [0,T/2], [T/2,T], ... [0,T/2^d], ...[(2^d-1)T/2^d,T]}
#'     \item Total number of partitions \eqn{npart= 2^{(d_part +1)} - 1}
#'   }
#' @param n_perturb Number of different perturbations on k-means result
#' @param perc_perturb Percentage of labels that are to be perturbed (= randomly switched)
#' @param n_random Number of completely random initial points. If not zero there will be n_random taus uniformly sampled in the initialization.
#' @param directed Boolean for directed (TRUE) or undirected (FALSE) case
#'
#' @return List of matrixes of initial values for \eqn{\tau}
#'
#' @keywords internal
#'
# # Examples
# # Generate initial tau for generated_Q3 data
#
# n <- 50
# Dmax <- 2^3
# Q <- 3
# d_part <- 1 # less than 3 (owing to Dmax)
# n_perturb <- 2
# perc_perturb <- 0.2
# n_random <- 1
# directed <- FALSE
#
# data <- list(Nijk = statistics(generated_Q3$data, n, Dmax, directed = FALSE))
#
# tau <- tauInitial(data,n,Q,d_part,n_perturb,perc_perturb,n_random,directed)
#
tauInitial <- function(data,n,Q,d_part,n_perturb,perc_perturb,n_random,directed){
  nbstart <- (2^(d_part+1)-1)*(1+n_perturb)
  init.values <- list()
  if (Q==1)
    init.values <- list(matrix(1,1,n))
  else{
    if (nbstart>0){
      nbswitch <- min(c(max(c(round(n*perc_perturb),2)),n))
      Dmax <- ncol(data$Nijk)
      #      nn <- min(c(2^d_part,Dmax))
      #      n_levels <- floor(log2(nn))
      n_levels <- min(c(d_part,log2(Dmax)))
      for (d_iter in 0:n_levels){
        statistics <- data$Nijk
        while(ncol(statistics)>2^d_iter){
          statistics <- as.matrix(statistics[,seq(1,ncol(statistics)-1,by=2)]+statistics[,seq(2,ncol(statistics),by=2)])
        }
        for (k in 1:2^d_iter){
          if (sum(statistics[,k]>0)>=Q){ # condition for kmeans to cluster into Q groups
            tau <- tauKmeansSbm(statistics[,k],n,Q,directed)
            tau <- apply(tau,2,correctTau)
            init.values <- append(init.values,list(tau))
            # perturbation by label switching
            if (n_perturb>0){
              for (l in 1:n_perturb){
                tau.copy <- tau
                label2switch <- sample(1:n,nbswitch)
                tau.copy[,label2switch] <- tau[,sample(label2switch,nbswitch)]
                tau.copy <- apply(tau.copy,2,correctTau)
                init.values <- append(init.values,list(tau.copy))
              }
            }
          }
        }# end loop for k in 1:2^d_iter
      } # end loop for d_iter
    }
    if (n_random>0){
      for (k in (nbstart+1):(nbstart+n_random)){
        tau <- matrix(rdirichlet(n,rep(1,Q)),Q,n,byrow=T)
        tau <- apply(tau,2,correctTau)
        init.values <- append(init.values,list(tau))
      }
    }
  } # end else Q!=1
  return(init.values)
}



#' Sparse setup - \eqn{\rho} parameter
#'
#' @param tau \eqn{\tau}
#' @param data Data : only needs the \eqn{N_{ijk}} field of data
#' @param n Total number of nodes
#' @param Q Total number of groups
#' @param directed Boolean for directed (TRUE) or undirected (FALSE) case
#'
#' @return Both \eqn{\tau} and \eqn{\rho}.
#'
#' @keywords internal
#'
# # Examples
# # Generate first initial tau for generated_Q3 data
#
# n <- 50
# Dmax <- 2^3
# Q <- 3
# d_part <- 1 # less than 3 (owing to Dmax)
# n_perturb <- 2
# perc_perturb <- 0.2
# n_random <- 1
# directed <- FALSE
#
# data <- list(Nijk = statistics(generated_Q3$data, n, Dmax, directed = FALSE))
#
# tau <- tauInitial(data,n,Q,d_part,n_perturb,perc_perturb,n_random,directed)
#
# taurho <- taurhoInitial(tau[[1]],data,n,Q,directed=FALSE)
#
taurhoInitial <- function(tau,data,n,Q,directed=TRUE){
  if (directed){
    N <- n*(n-1)
    N_Q <- Q^2
  }
  else{ # undirected
    N <- n*(n-1)/2
    N_Q <- Q*(Q+1)/2
  }
  ind.all <- listNodePairs(n,directed)
  Nij <- rowSums(data$Nijk)
  rho <- 1:N_Q
  for (ind.ql in 1:N_Q){
    if (Q==1)
      tauql <- rep(1,N)
    else{
      ql <- find_ql(ind.ql,Q,directed)
      q <- ql[1]
      l <- ql[2]
      if (directed) # directed case
        tauql <- tau[q, ind.all[,1]]*tau[l, ind.all[,2]] # N-vector
      else{ # undirected
        if ((q==l))
          tauql <- tau[q, ind.all[,1]]*tau[l, ind.all[,2]]
        else
          tauql <- tau[q, ind.all[,1]]*tau[l, ind.all[,2]] + tau[q, ind.all[,2]]*tau[l, ind.all[,1]]
      }
    }
    nb.procql.connect <- sum((Nij>0)*tauql)
    Yql <- sum(tauql)
    beta <- nb.procql.connect/Yql
    AqlT <- sum(Nij*tauql)/nb.procql.connect
    rho[ind.ql] <- 1/((1/beta-1)*exp(AqlT)+1)
  }
  return(list(tau=tau,rho=rho))
}


#' Construct initial \eqn{\tau} from \eqn{Q-1}
#'
#' Construct initial \eqn{\tau} with \eqn{Q} groups from value obtained at \eqn{Q-1} groups
#'
#'
#' @param tau \eqn{\tau}
#' @param n_perturb  Number of different perturbations on k-means result
#'
#' @return List of matrixes of initial values for \eqn{\tau} for \eqn{Q} groups from value obtained at \eqn{Q-1}
#'
#' @keywords internal
#'
# # Examples
# # Generate first initial tau for generated_Q3 data
#
# n <- 50
# Dmax <- 2^3
# Q <- 3
# d_part <- 1 # less than 3 (owing to Dmax)
# n_perturb <- 2
# perc_perturb <- 0.2
# n_random <- 1
# directed <- FALSE
#
# data <- list(Nijk = statistics(generated_Q3$data, n, Dmax, directed = FALSE))
#
# tau <- tauInitial(data,n,Q,d_part,n_perturb,perc_perturb,n_random,directed)
#
# tau.list <- tauUp_Q(tau[[1]],1)
#
tauUp_Q <- function(tau,n_perturb=1){
  Q <- dim(tau)[1]
  n <- dim(tau)[2]
  tau.list <- list()

  if (Q==1){
    tau.list <- list(matrix(1/2,2,n))
  } else{ # Q>=2
    if (n_perturb>=Q){ # split all the groups once
      for (q in 1:Q){
        tau.new <- tau
        tau.new[q,] <- tau.new[q,]/2
        tau.new <- rbind(tau.new,tau.new[q,])
        tau.list <- append(tau.list,list(tau.new))
      }
    } else { # always split the largest group
      largest.gp <- which.max(rowSums(tau))
      tau.new <- tau
      tau.new[largest.gp,] <- tau.new[largest.gp,]/2
      tau.new <- rbind(tau.new,tau.new[largest.gp,])
      tau.list <- append(tau.list,list(tau.new))
      if (n_perturb >=2){
        # and split the component with largest entropy
        ent <- tau*log(tau)
        ent[is.na(ent)] <- 0
        ent.gp <- which.min(rowSums(ent))
        tau.new <- tau
        tau.new[ent.gp,] <- tau.new[ent.gp,]/2
        tau.new <- rbind(tau.new,tau.new[ent.gp,])
        tau.list <- append(tau.list,list(tau.new))
      }
      # And split a number n_perturb-2 of groups which are randomly chosen
      if (n_perturb >2){
        split.gp <- sample((1:Q)[-c(largest.gp,ent.gp)],n_perturb-2,replace=F)
        for (q in split.gp){
          tau.new <- tau
          tau.new[q,] <- tau.new[q,]/2
          tau.new <- rbind(tau.new,tau.new[q,])
          tau.list <- append(tau.list,list(tau.new))
        }
      }
    } # end if/else n_perturb>=Q
  } # end if/else Q==1

  return(tau.list)
}


#' Construct initial \eqn{\tau} from \eqn{Q+1}
#'
#' Construct initial \eqn{\tau} with \eqn{Q} groups from value obtained at \eqn{Q+1} groups
#'
#' @param tau \eqn{\tau}
#' @param n_perturb  Number of different perturbations on k-means result
#'
#' @return List of matrixes of initial values for \eqn{\tau} for \eqn{Q} groups from value obtained at \eqn{Q+1}
#'
#' @keywords internal
#'
# # Examples
# # Generate first initial tau for generated_Q3 data
#
# n <- 50
# Dmax <- 2^3
# Q <- 3
# d_part <- 1 # less than 3 (owing to Dmax)
# n_perturb <- 2
# perc_perturb <- 0.2
# n_random <- 1
# directed <- FALSE
#
# data <- list(Nijk = statistics(generated_Q3$data, n, Dmax, directed = FALSE))
#
# tau <- tauInitial(data,n,Q,d_part,n_perturb,perc_perturb,n_random,directed)
#
# tau.list <- tauDown_Q(tau[[1]],1)
#
tauDown_Q <- function(tau,n_perturb=1){
  Q <- dim(tau)[1]
  n <- dim(tau)[2]
  tau.list <- list()

  if (Q==2){
    tau.list <- list(matrix(1,1,n))
  } else {
    if (n_perturb>=Q*(Q-1)/2){ # merge all the possible pairs of groups (q,l) with q<l
      for (q in 1:(Q-1)){
        for (l in ((q+1):Q)){
          tau.new <- tau
          tau.new <- rbind(tau.new[-c(q,l),],tau.new[q,]+tau.new[l,])
          tau.list <- append(tau.list,list(tau.new))
        }
      }
    } else {
      # start by merging the 2 components with largest entropies
      ent <- tau*log(tau)
      ent[is.na(ent)] <- 0
      entropy.per.gp <- rowSums(ent)
      q <- order(entropy.per.gp)[1]
      l <- order(entropy.per.gp)[2]
      tau.new <- tau
      tau.new <- rbind(tau.new[-c(q,l),],tau.new[q,]+tau.new[l,])
      tau.list <- append(tau.list,list(tau.new))
      # merge a number n_perturb-1 of pairs of groups (q,l) q<l which are randomly chosen
      merge.gp <- sample(1:(Q*(Q-1)/2),n_perturb-1,replace=F)
      for (ind.ql in merge.gp){
        tau.new <- tau
        q <- find_ql_diff(ind.ql,Q)[1]  # this is defined in auxiliary.R
        l <- find_ql_diff(ind.ql,Q)[2]
        tau.new <- rbind(tau.new[-c(q,l),],tau.new[q,]+tau.new[l,])
        tau.list <- append(tau.list,list(tau.new))
      }
    } # end if/else n_perturb>=Q*(Q-1)/2
  } # end if/else Q==2
  return(tau.list)
}



######################
####   K-MEANS FOR SBM
######################


#' Function for k-means
#'
#' @param cl Label list of nodes
#'
#' @return x    : class indicator matrix
#'
#' @keywords internal
#'
classInd<-function (cl){
  ## INPUT
  ##    cl    : label list of nodes
  ## OUPUT
  ##    x    : class indicator matrix
  ## classInd builds a class indicator matrix
  ## from a list of labels
  ## nq  <- rmultinom(1, size=7, prob = c(1/2, 1/8, 1/8))
  ## cl  <- rep(1:3, nq)
  ## cl = [1] 1 1 1 1 2 2 3
  ## classInd(cl)
  ##      1 2 3
  ## [1,] 1 0 0
  ## [2,] 1 0 0
  ## [3,] 1 0 0
  ## [4,] 1 0 0
  ## [5,] 0 1 0
  ## [6,] 0 1 0
  ## [7,] 0 0 1

  n         <- length(cl)
  cl         <- as.factor(cl)
  x         <- matrix(0, n, length(levels(cl)))
  x[(1:n) + n * (unclass(cl) - 1)] <- 1
  dimnames(x) <- list(names(cl), levels(cl))
  x
}


####
## k-means for SBM
####

#' k-means for SBM
#'
#' @param statistics Statistics matrix \eqn{N_{ijk}}, counting the events for the nodes pair \eqn{(i,j)} during the subinterval \eqn{k}
#' @param n Total number of nodes \eqn{n}
#' @param Q Total number of groups \eqn{Q}
#' @param directed Boolean for directed (TRUE) or undirected (FALSE) case
#'
#' @return Initial values for \eqn{\tau}
#'
#' @keywords internal
#'
# #Examples
#
# n <- 50
# Q <- 3
#
# Dmax <- 2^3
#
# Nijk <- statistics(generated_Q3$data,n,Dmax,directed=FALSE)
#
# tau <- tauKmeansSbm(Nijk,n,Q,FALSE)
#
tauKmeansSbm <- function(statistics,n,Q,directed){
  X   <- matrix(0,n,n)
  if (!directed){
    k <- 1
    for (i in 1:(n-1)){
      X[i,(i+1):n] <- statistics[k:(k+n-i-1)]
      k <- k+n-i
    }
    X <- X + t(X)
  }
  else{
    for (i in 1:n){
      X[i,-i] <- statistics[((i-1)*(n-1)+1):(i*(n-1))]
    }
  }

  try.km <- try(kmeans(X,Q,nstart=50),silent=T)
  if (is.list(try.km)) {
    tau <- t(classInd(try.km$cluster))
  } else {
    tau <- matrix(rdirichlet(n,rep(1,Q)),Q,n,byrow=T)}
  tau <- apply(tau,2,correctTau)

  return(tau)
}

