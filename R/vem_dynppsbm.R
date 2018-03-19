############################################################
## VEM FOR ADAPTIVE NONPARAMETRIC ESTIMATION IN DYNPPSBM  ##
############################################################



###################################################
## ADAPTIVE VEM ALGORITHM
###################################################

#' Adaptative VEM algorithm
#'
#' Principal adaptative VEM algorithm for histogram with model selection or for kernel method.
#'
#' The sparse version works only for the histogram approach.
#'
#' @param data Data format depends on the estimation method used!!
#' \enumerate{
#'   \item Data with \strong{hist} method - list with 2 components:
#'     \describe{
#'       \item{data$Time}{[0,data$Time] is the total time interval of observation}
#'
#'       \item{data$Nijk}{Data matrix with counts per process \eqn{N_{ij}} and sub-intervals ; matrix of size \eqn{N*Dmax}  where \eqn{N = n(n-1)} or \eqn{n(n-1)/2} is the number of possible node pairs in the graph and \eqn{Dmax = 2^{dmax}} is the size of the finest partition in the  histrogram approach
#'
#'       Counts are pre-computed - Obtained through function 'statistics' (auxiliary.R) on data with second format}
#'     }
#'
#'   \item Data with \strong{kernel} method - list with 3 components:
#'     \describe{
#'       \item{data$time.seq}{Sequence of observed time points of the m-th event (M-vector)}
#'
#'       \item{data$type.seq}{Sequence of observed values convertNodePair(i,j,n,directed) (auxiliary.R) of process that produced the mth event (M-vector).}
#'
#'     \item{data$Time}{[0,data$Time] is the total time interval of observation}
#'     }
#' }
#' @param n Total number of nodes
#' @param Qmin Minimum number of groups
#' @param Qmax Maximum number of groups
#' @param directed Boolean for directed (TRUE) or undirected (FALSE) case
#' @param sparse Boolean for sparse (TRUE) or not sparse (FALSE) case
#' @param method List of string. Can be "hist" for histogram method or "kernel" for kernel method
#' @param init.tau List of initial values of \eqn{\tau} - all tau's are matrices with size \eqn{Q\times n} (might be with different values of Q)
#' @param cores Number of cores for parallel execution
#'
#' If set to 1 it does sequential execution
#'
#' Beware: parallelization with fork (multicore) : doesn't work on Windows!
#'
#' @param d_part Maximal level for finest partition of time interval [0,T] used for k-means initializations.
#'   \itemize{
#'     \item Algorithm takes partition up to depth \eqn{2^d} with \eqn{d=1,...,d_{part}}
#'     \item Explore partitions \eqn{[0,T], [0,T/2], [T/2,T], ... [0,T/2^d], ...[(2^d-1)T/2^d,T]}
#'     \item Total number of partitions \eqn{npart= 2^{(d_{part} +1)} -1}
#'   }
#'
#' @param n_perturb Number of different perturbations on k-means result
#'
#'       When \eqn{Qmin < Qmax}, number of perturbations on the result with \eqn{Q-1} or \eqn{Q+1} groups
#' @param perc_perturb Percentage of labels that are to be perturbed (= randomly switched)
#' @param n_random Number of completely random initial points. The total number of initializations for the VEM is \eqn{npart*(1+n_{perturb}) +n_{random}}
#' @param nb.iter Number of iterations of the VEM algorithm
#' @param fix.iter Maximum number of iterations of the fixed point into the VE step
#' @param epsilon Threshold for the stopping criterion of VEM and fixed point iterations
#' @param filename Name of the file where to save the results along the computation (increasing steps for \eqn{Q}, these are the longest).
#'
#'    The file will contain a list of 'best' results.
#'
#' @export
#'
#' @references
#'
#' DAUDIN, J.-J., PICARD, F. & ROBIN, S. (2008). A mixture model for random graphs. Statist. Comput. 18, 173–183.
#'
#' DEMPSTER, A. P., LAIRD, N. M. & RUBIN, D. B. (1977). Maximum likelihood from incomplete data via the EM algorithm. J. Roy. Statist. Soc. Ser. B 39, 1–38.
#'
#' JORDAN, M., GHAHRAMANI, Z., JAAKKOLA, T. & SAUL, L. (1999). An introduction to variational methods for graphical models. Mach. Learn. 37, 183–233.
#'
#' MATIAS, C., REBAFKA, T. & VILLERS, F. (2018).  A semiparametric extension of the stochastic block model for longitudinal networks. Biometrika.
#'
#' MATIAS, C. & ROBIN, S. (2014). Modeling heterogeneity in random graphs through latent space models: a selective review. Esaim Proc. & Surveys 47, 55–74.
#'
#' @examples
#' # load data of a synthetic graph with 50 individuals and 3 clusters
#' n <- 20
#' Q <- 3
#'
#' Time <- generated_Q3_n20$data$Time
#' data <- generated_Q3_n20$data
#' z <- generated_Q3_n20$z
#'
#' step <- .001
#' x0 <- seq(0,Time,by=step)
#' intens <-  generated_Q3_n20$intens
#'
#' # VEM-algo kernel
#' sol.kernel <- mainVEM(data,n,Q,directed=FALSE,method='kernel', d_part=0,
#'     n_perturb=0)[[1]]
#' # compute smooth intensity estimators
#' sol.kernel.intensities <- kernelIntensities(data,sol.kernel$tau,Q,n,directed=FALSE)
#' # eliminate label switching
#' intensities.kernel <- sortIntensities(sol.kernel.intensities,z,sol.kernel$tau,
#'     directed=FALSE)
#'
#' # VEM-algo hist
#' # compute data matrix with precision d_max=3
#' Dmax <- 2^3
#' Nijk <- statistics(data,n,Dmax,directed=FALSE)
#' sol.hist <- mainVEM(list(Nijk=Nijk,Time=Time),n,Q,directed=FALSE, method='hist',
#'     d_part=0,n_perturb=0,n_random=0)[[1]]
#' log.intensities.hist <- sortIntensities(sol.hist$logintensities.ql,z,sol.hist$tau,
#'      directed=FALSE)
#'
#' # plot estimators
#' par(mfrow=c(2,3))
#' ind.ql <- 0
#' for (q in 1:Q){
#'   for (l in q:Q){
#'     ind.ql <- ind.ql + 1
#'     true.val <- intens[[ind.ql]]$intens(x0)
#'     values <- c(intensities.kernel[ind.ql,],exp(log.intensities.hist[ind.ql,]),true.val)
#'     plot(x0,true.val,type='l',xlab=paste0("(q,l)=(",q,",",l,")"),ylab='',
#'         ylim=c(0,max(values)+.1))
#'     lines(seq(0,1,by=1/Dmax),c(exp(log.intensities.hist[ind.ql,]),
#'         exp(log.intensities.hist[ind.ql,Dmax])),type='s',col=2,lty=2)
#'     lines(seq(0,1,by=.001),intensities.kernel[ind.ql,],col=4,lty=3)
#'   }
#' }
#'
mainVEM <- function(data,
                     n,
                     Qmin,
                     Qmax=Qmin,
                     directed=TRUE,
                     sparse=FALSE,
                     method=c('hist','kernel'),
                     init.tau=NULL,
                     cores=1,
                     d_part=5,
                     n_perturb=10,
                     perc_perturb=.2,
                     n_random=0,
                     nb.iter=50,
                     fix.iter=10,
                     epsilon=1e-6,
                     filename=NULL){

  #### Sanity checks and default settings
  if (is.null(data$Time))
    stop("You didn't specify time interval")

  if (missing(method)) # default is 'hist' method
    method <- 'hist'

  if (method=='hist'){
    if (is.null(data$Nijk))
      stop("Compute 'statistics' before running vem with 'hist' (default) method")
  }

  if (method=='kernel'){
    if (is.null(data$time.seq))
      stop("Check input format for 'kernel' inference method")
    if (is.null(data$type.seq))
      stop("Check input format for 'kernel' inference method")
    if (sparse)
      stop("'kernel' inference method not implemented for sparse model")

    # Compute statistics to be used in initialisations
    data$Nijk <- statistics(data,n,2^d_part,directed)
  }

  if (missing(n))
    stop("Number of individuals is missing")

  if (missing(Qmin))
    stop("Number of groups is missing")
  #### end of sanity checks

  best <- list(NULL,Qmax-Qmin+1)
  for (Q in Qmin:Qmax){ # Start loop for increasing values of Q
    cat("----increasing Q----",Q,"\n")
    # get a list of initial points for tau
    init.values <- tauInitial(data,n,Q,d_part,n_perturb,perc_perturb,n_random,directed)
    if (!(is.null(init.tau))) {
      for (i in 1:length(init.tau)){
        if (dim(init.tau[[i]])[1]==Q)
        {
          init.values <- append(init.values,list(init.tau[[i]]))
        }
      }
    }
    if (Q>Qmin) {# use perturbed result with Q-1 groups as additional init points
      init.values <- append(init.values,tauUp_Q(best[[Q-Qmin]]$tau,n_perturb))
    }

    # start the different runs
    if (cores==1){
      best[[Q-Qmin+1]] <- list(J=-Inf)
      for (run in 1:length(init.values)){
        cat("----run----",run,"\n")
        # Initialisation of tau
        if (sparse){ # sparse
          VE <- taurhoInitial(init.values[[run]],data,n,Q,directed)
        } else { # not sparse
          VE <- list(tau=init.values[[run]],rho=NULL)}

        convergence <- list(converged=FALSE, nb.iter.achieved=FALSE)
        J.old <- -Inf
        it <- 0
        while (sum(convergence==T)==0){
          it <- it+1
          cat("iteration",it,'\n')
          if (method=='hist'){ # histogram
            mstep <- Mstep_hist(data,VE,directed,sparse)
          } else { # kernel
            mstep <- Mstep_kernel(data,VE,directed)}
          VE <-	 VEstep(VE,mstep,directed,sparse,method,epsilon,fix.iter,data)
          J <- JEvalMstep(VE,mstep,data,directed,sparse,method)$J
          convergence$nb.iter.achieved <- (it > nb.iter+1)
          convergence$converged <-  (abs((J-J.old)/J)< epsilon)
          J.old    <- J
        }

        if (method=='hist'){ cat("values of d:",mstep$best.d,"\n")}

        cat("value of J:",J,"\n")
        if (J > best[[Q-Qmin+1]]$J)
          if (method=='hist'){
            best[[Q-Qmin+1]] <- list(tau=VE$tau,rho=VE$rho,beta=mstep$beta,logintensities.ql=mstep$logintensities.ql,best.d=mstep$best.d,J=J,run=run,converged= convergence$converged)
          } else { # method=='kernel' - note that sparse inference not implemented
            best[[Q-Qmin+1]] <- list(tau=VE$tau,logintensities.ql.ij=mstep$logintensities.ql.ij,J=J,run=run,converged= convergence$converged)
          }

      } # end loop for run
    }
    else{ # parallelization
      sol <- mclapply(init.values, function(init.point){
        mainVEMPar(init.point,n,Q,data,directed,sparse,method,nb.iter,fix.iter,epsilon)
      },mc.cores=cores)

      sol.J <- lapply(sol,function(x)return(x$J))
      best[[Q-Qmin+1]] <- sol[[which.max(sol.J)]]
    }

    if (!is.null(filename))
      save(best,file=filename)
  }# end of increasing loop for Q


  if ((Qmin<Qmax)&(n_perturb>0)){ # do a loop for decreasing values of Q
    for (Q in (Qmax-1):Qmin){ # use perturbed result with Q+1 groups as init points
      cat("----decreasing Q----",Q,"\n")
      init.values <- tauDown_Q(best[[Q-Qmin+2]]$tau,n_perturb)
      # start the different runs
      if (cores==1){
        for (run in 1:length(init.values)){
          # encode run with negative integer to distinguish it from increasing phase
          cat("----run----",run,"-\n")
          # Initialisation of tau
          if (sparse){
            VE <- taurhoInitial(init.values[[run]],data,n,Q,directed)
          } else {VE <- list(tau=init.values[[run]],rho=NULL)}

          convergence <- list(converged=FALSE, nb.iter.achieved=FALSE)
          J.old <- -Inf
          it <- 0
          while (sum(convergence==T)==0){
            it <- it+1
            cat("iteration",it,'\n')
            if (method=='hist'){
              mstep <- Mstep_hist(data,VE,directed,sparse)
            } else {mstep <- Mstep_kernel(data,VE,directed)}
            VE <-	 VEstep(VE,mstep,directed,sparse,method,epsilon,fix.iter,data)
            J <- JEvalMstep(VE,mstep,data,directed,sparse,method)$J
            convergence$nb.iter.achieved <- (it > nb.iter+1)
            convergence$converged <-  (abs((J-J.old)/J)< epsilon)
            J.old    <- J
          }

          if (method=='hist'){ cat("values of d:",mstep$best.d,"\n")}

          cat("value of J:",J,"\n")
          if (J > best[[Q-Qmin+1]]$J)
            if (method=='hist'){
              # encode run with negative integer to distinguish from increasing phase
              best[[Q-Qmin+1]] <- list(tau=VE$tau,rho=VE$rho,beta=mstep$beta,logintensities.ql=mstep$logintensities.ql,best.d=mstep$best.d,J=J,run=-run,converged= convergence$converged)
            } else { # method=='kernel' - note that sparse inference not implemented
              best[[Q-Qmin+1]] <- list(tau=VE$tau,logintensities.ql.ij=mstep$logintensities.ql.ij,J=J,run=-run,converged= convergence$converged)
            }

        } # end loop for run
      }
      else{ # parallelization
        sol <- mclapply(init.values, function(init.point){
          mainVEMPar(init.point,n,Q,data,directed,sparse,method,nb.iter,fix.iter,epsilon)
        },mc.cores=cores)

        sol.J <- lapply(sol,function(x)return(x$J))
        if ( max(unlist(sol.J)) > best[[Q-Qmin+1]]$J){
          best[[Q-Qmin+1]] <- sol[[which.max(sol.J)]]
        }
      }

      if (!is.null(filename))
        save(best,file=filename)
    }# end decreasing loop
  }

  return(best)
}

###################################################
## MAIN VEM - STEPS
###################################################

# init.point <- init.values[[run]]

#' VEM step for parallel version
#'
#' @param init.point Initial point
#' @param n Total number of nodes
#' @param Q Total number of groups
#' @param data Data same of \link[ppsbm]{mainVEM}
#' @param directed  Boolean for directed (TRUE) or undirected (FALSE) case
#' @param sparse Boolean for sparse (TRUE) or not sparse (FALSE) case
#' @param method List of string. Can be "hist" for histogram method or "kernel" for kernel method
#' @param nb.iter Number of iterations
#' @param fix.iter Maximum number of iterations of the fixed point
#' @param epsilon Threshold for the stopping criterion of VEM and fixed point iterations
#'
#' @export
#'
mainVEMPar <- function(init.point,n,Q,data,directed,sparse,method,nb.iter,fix.iter,epsilon){
  # Initialisation of tau
  if (sparse){
    VE <- taurhoInitial(init.point,data,n,Q,directed)
  } else {VE <- list(tau=init.point,rho=NULL)}

  convergence <- list(converged=FALSE, nb.iter.achieved=FALSE)
  J.old <- -Inf
  it <- 0
  while (sum(convergence==T)==0){
    it <- it+1
    if (method=='hist'){
      mstep <- Mstep_hist(data,VE,directed,sparse)
    } else {mstep <- Mstep_kernel(data,VE,directed)}
    VE <-	 VEstep(VE,mstep,directed,sparse,method,epsilon,fix.iter,data)
    J <- JEvalMstep(VE, mstep,data,directed,sparse,method)$J
    convergence$nb.iter.achieved <- (it > nb.iter+1)
    convergence$converged <-  (abs((J-J.old)/J)< epsilon)
    J.old    <- J
  }

  if (method=='hist'){
    sol.run <- list(tau=VE$tau,rho=VE$rho,beta=mstep$beta,logintensities.ql=mstep$logintensities.ql,best.d=mstep$best.d,J=J,converged= convergence$converged)
  } else { # method='kernel' - note that sparse inference not implemented
    sol.run <- list(tau=VE$tau,logintensities.ql.ij=mstep$logintensities.ql.ij,J=J,converged= convergence$converged)
  }

  return(sol.run)
}



###################################################
## VE-STEP
###################################################

#' VE step
#'
#' @param VE Results of the previous VE step for iterative computation
#' @param mstep Results of the previous mstep for iterative computation
#' @param directed Boolean for directed (TRUE) or undirected (FALSE) case
#' @param sparse Boolean for sparse (TRUE) or not sparse (FALSE) case
#' @param method List of string. Can be "hist" for histogram method or "kernel" for kernel method
#' @param epsilon Threshold for the stopping criterion of VEM and fixed point iterations
#' @param fix.iter Maximum number of iterations of the fixed point
#' @param data Data same of \link[ppsbm]{mainVEM}
#'
#'
VEstep <- function(VE,mstep,directed,sparse,method,epsilon,fix.iter,data){
  # variables
  Q <- dim(VE$tau)[[1]]

  # update rho in sparse case
  rho <- if (sparse) 1/((1/mstep$beta-1)*exp(mstep$AqlT)+1) else 0

  # compute pi
  pi <- if (Q!=1) apply(VE$tau,1,mean) else 1

  # solve the fixed point equation
  tau <- VE$tau
  converged <- FALSE
  it <- 0
  while ((!converged)&(it<fix.iter)){
    it <- it+1
    tau.old <- tau
    tau <- tauUpdate(tau,pi,mstep,data,directed,sparse,method,rho)
    converged <- ((max(abs(tau-tau.old))<epsilon))
  }

  return(list(tau=tau,rho=rho))
}

###############
## one update of tau by the fixed point equation
###############

#' Update \eqn{\tau}
#'
#' One update of \eqn{\tau} by the fixed point equation
#'
#' @param tau Old \eqn{\tau}
#' @param pi Estimator of group probabilities \eqn{\pi}
#' @param mstep Results of the previous mstep for iterative computation
#' @param data Data same of \link[ppsbm]{mainVEM}
#' @param directed  Boolean for directed (TRUE) or undirected (FALSE) case
#' @param sparse Boolean for sparse (TRUE) or not sparse (FALSE) case
#' @param method List of string. Can be "hist" for histogram method or "kernel" for kernel method
#' @param rho Old \eqn{\rho} (only for sparse model, set to 0 otherwise)
#'
#'
tauUpdate <- function (tau,pi,mstep,data,directed,sparse,method,rho){
  Q <- dim(tau)[1]
  if ((method=='hist') & (sparse)){
    no.obs <- (rowSums(data$Nijk)==0)
  }

  if (Q==1)
    tau.new <- tau
  else{
    n <- dim(tau)[2]
    logtau <- log(tau)
    tau.new <- tau

    if (sparse){
      psi_rho <- rho*log(rho)+(1-rho)*log(1-rho)# N_Q - vector
      psi_rho[is.na(psi_rho)] <- 0
    }

    if (directed){ # directed
      for (i in 1:n){
        j.in.ij <- convertNodePair(rep(i,1,n-1),(1:n)[-i],n,directed)	 # n-1
        j.in.ji <- convertNodePair((1:n)[-i],rep(i,1,n-1),n,directed)	 # n-1
        sum_tauj_Nijk <- tau[,-i]%*%as.matrix(data$Nijk[j.in.ij,]) # QxDmax
        sum_tauj_Njik <- tau[,-i]%*%as.matrix(data$Nijk[j.in.ji,]) # QxDmax

        for (q in 1:Q){
          ind.ql <- convertGroupPair(rep(q,1,Q),1:Q,Q,directed)  # indices of all groups (q,l) for l=1,...,Q
          ind.lq <- convertGroupPair(1:Q,rep(q,1,Q),Q,directed)  # indices of all groups (l,q) for l=1,...,Q

          if (method=='hist'){
            term1 <- sum(mstep$logintensities.ql[ind.ql,]*sum_tauj_Nijk)
            term1 <- term1 + sum(mstep$logintensities.ql[ind.lq,]*sum_tauj_Njik)
          } else{ # kernel
            term1 <- sum(mstep$logintensities.ql.ij[ind.ql,j.in.ij])
            term1 <- term1 + sum(mstep$logintensities.ql.ij[ind.lq,j.in.ji])
          }


          ###
          AqlT <- mstep$AqlT[ind.ql] # Q
          AlqT <- mstep$AqlT[ind.lq]
          nonzero <- (mstep$sum_rhotau[ind.ql]!=0)&(mstep$sum_rhotau[ind.lq]!=0) # Q

          if (sparse){
            rho.mat1 <- matrix(1,Q,n-1)
            rho.mat1[,no.obs[j.in.ij]] <- rho[ind.ql]
            sum_taujl_rhoijql <- rowSums(tau[,-i]*rho.mat1) # Q
            rho.mat2 <- matrix(1,Q,n-1)  # Q x (n-1)
            rho.mat2[,no.obs[j.in.ji]] <- rho[ind.lq]
            sum_taujl_rhojilq <- rowSums(tau[,-i]*rho.mat2) # Q  (sum over j)
            term2 <- sum((sum_taujl_rhoijql*AqlT+sum_taujl_rhojilq*AlqT)[nonzero])

            term3 <- sum(log(mstep$beta[ind.ql]) * sum_taujl_rhoijql)
            term3 <- term3 + sum(log(1-mstep$beta[ind.ql]) * rowSums(tau[,-i]*(1-rho.mat1)))

            sum_tau_jl <- rowSums(tau[,-i]*matrix(no.obs[j.in.ij],nrow=Q,ncol=n-1,byrow=T))
            term4 <- sum(psi_rho[ind.ql]*sum_tau_jl)
            sum_tau_jl <- rowSums(tau[,-i]*matrix(no.obs[j.in.ji],nrow=Q,ncol=n-1,byrow=T))
            term4 <- term4 + sum(psi_rho[ind.lq]*sum_tau_jl)

          } else{ # non sparse
            sum_taujl <- rowSums(tau[,-i])
            term2 <- sum((sum_taujl*(AqlT+AlqT))[nonzero])
            term3 <- 0
            term4 <- 0
          }

          logtau[q,i] <- log(pi[q]) + term1 - term2 + term3 - term4
        }
      }
      ## Normalizing in the log space to avoid numerical problems
      ## and going back to exponential with the same normalization
      logtau[,i]<-logtau[,i]-max(logtau[,i])
      tau.new[,i] <- exp(logtau[,i])
      tau.new[,i] <- correctTau(tau.new[,i])

    }

    else{ # undirected
      for (i in 1:n){
        j.in.ij <- convertNodePair(rep(i,1,n-1),(1:n)[-i],n,directed)	 # n-1
        sum_tauj_Nijk <- tau[,-i]%*%as.matrix(data$Nijk[j.in.ij,]) # QxDmax

        for (q in 1:Q){
          ind.ql <- convertGroupPair(rep(q,1,Q),1:Q,Q,directed)  # indices of all groups (q,l) for l=1,...,Q
          if (method=='hist'){
            term1 <- sum(mstep$logintensities.ql[ind.ql,]*sum_tauj_Nijk)
          } else{ # kernel
            term1 <- sum(mstep$logintensities.ql.ij[ind.ql,j.in.ij])
          }


          ###
          AqlT <- mstep$AqlT[ind.ql] # Q
          nonzero <- (mstep$sum_rhotau[ind.ql]!=0) # Q

          if (sparse){
            rho.mat1 <- matrix(1,Q,n-1) # Q x (n-1)
            rho.mat1[,no.obs[j.in.ij]] <- rho[ind.ql]
            sum_taujl_rhoijql <- rowSums(tau[,-i]*rho.mat1) # Q
            term2 <- sum((sum_taujl_rhoijql*AqlT)[nonzero])

            term3 <- sum(log(mstep$beta[ind.ql]) * sum_taujl_rhoijql)

            sum_tau_jl <- rowSums(tau[,-i]*matrix(no.obs[j.in.ij],nrow=Q,ncol=n-1,byrow=T))
            term4 <- sum(psi_rho[ind.ql]*sum_tau_jl)

          } else{ # non sparse
            sum_taujl <- rowSums(tau[,-i])
            term2 <- sum((sum_taujl*AqlT)[nonzero])
            term3 <- 0
            term4 <- 0
          }

          logtau[q,i] <- log(pi[q]) + term1 - term2 + term3 - term4
        }
      }
      ## Normalizing in the log space to avoid numerical problems
      ## and going back to exponential with the same normalization
      logtau[,i]<-logtau[,i]-max(logtau[,i])
      tau.new[,i] <- exp(logtau[,i])
      tau.new[,i] <- correctTau(tau.new[,i])
    }

  }
  return(tau.new)
}

###################################################
##   M-STEP for histograms
###################################################


#' M step for histograms
#'
#' M step for histograms estimator
#'
#' @param data Data same of \link[ppsbm]{mainVEM}
#' @param VE Results of the previous VE for iterative computation
#' @param directed Boolean for directed (TRUE) or undirected (FALSE) case
#' @param sparse  Boolean for sparse (TRUE) or not sparse (FALSE) case
#'
#' @references
#'
#' BARAUD, Y. & BIRGÉ, L. (2009). Estimating the intensity of a random measure by histogram type estimators. Probab. Theory Related Fields 143, 239–284.
#'
#' MATIAS, C., REBAFKA, T. & VILLERS, F. (2018).  A semiparametric extension of the stochastic block model for longitudinal networks. Biometrika.
#'
#' REYNAUD -BOURET, P. (2006). Penalized projection estimators of the Aalen multiplicative intensity. Bernoulli 12, 633–661.
#'
Mstep_hist <- function(data,VE,directed,sparse){
  Q <- dim(VE$tau)[1]
  n <- dim(VE$tau)[2]

  N <- if (directed) n*(n-1) else n*(n-1)/2
  N_Q <- if (directed) Q^2 else Q*(Q+1)/2

  Dmax <- ncol(data$Nijk)
  D <- log2(Dmax)
  logintensities.ql <- matrix(0,N_Q,Dmax)
  sum_rhotau <- rep(0,N_Q)
  sum_rhotau_obs <- rep(0,N_Q)
  AqlT <- rep(0,N_Q)
  beta <- rep(0,N_Q) # return 0 values in non sparse case

  if (sparse){
    no.obs <- (rowSums(data$Nijk)==0)
  }

  crit_pen <- rep(Inf,N_Q)
  pen_cste <- rep(NA,N_Q)
  Nijk.d <- data$Nijk
  best.d <- rep(NA,N_Q)

  if (directed){ # directed case
    ind.all <- listNodePairs(n)
    for (d in D:0){
      ind.ql <- 0
      for (q in 1:Q){
        for (l in 1:Q){
          ind.ql <- ind.ql + 1
          if (Q==1)
            tauql <- rep(1,N)
          else{
            tauql <- VE$tau[q, ind.all[,1]]*VE$tau[l, ind.all[,2]] # N-vector
          }
          if (sparse){ # sparse
            rho.vec <- rep(1,N)
            rho.vec[no.obs] <- VE$rho[ind.ql]
            tauql_rho <- tauql*rho.vec
            wN_k <- colSums(Nijk.d*matrix(tauql_rho,N,2^d))  # N-vector
          } else{ # not sparse
            wN_k <- colSums(Nijk.d*matrix(tauql,N,2^d)) # N-vector
          }

          if (d==D) {
            pen_cste[ind.ql] <- 2*Dmax*max(wN_k)

            if (sparse){ # sparse
              sum_tauql <- sum(tauql)
              sum_rhotau_obs[ind.ql] <- sum(wN_k)
              sum_rhotau[ind.ql] <- sum(tauql_rho)
              if (sum_tauql>0)
                beta[ind.ql] <- sum_rhotau[ind.ql]/sum_tauql
            } else{ # non sparse
              sum_rhotau[ind.ql] <-  sum(tauql)
              sum_rhotau_obs[ind.ql] <- sum(wN_k)}
          }

          if (sum_rhotau[ind.ql]>0)
            AqlT[ind.ql] <- sum_rhotau_obs[ind.ql]/sum_rhotau[ind.ql]

          crit_pen_d = 2^d*(-sum(wN_k^2) + pen_cste[ind.ql])
          if (crit_pen_d<crit_pen[ind.ql]){
            best.d[ind.ql] <- 2^d
            crit_pen[ind.ql] <- crit_pen_d
            if (sum_rhotau[ind.ql]>0){
              loglam <- log(2^d/sum_rhotau[ind.ql]*wN_k/data$Time)
              loglam[loglam==-Inf] <- 0
              logintensities.ql[ind.ql,] <- rep(loglam,each=2^(D-d))
            }

          } # end case for criterion

        } # end for l
      }# end for q

      if (d>0)
        Nijk.d <- Nijk.d[,seq(1,2^d-1,by=2)]+Nijk.d[,seq(2,2^d,by=2)]

    } # end for d
  } else{ # undirected
    ind.all <- listNodePairs(n,directed)
    for (d in D:0){
      ind.ql <- 0
      for (q in 1:Q){
        for (l in q:Q){
          ind.ql <- ind.ql + 1
          if (Q==1)
            tauql <- rep(1,N)
          else{ # undirected
            if ((q==l))
              tauql <- VE$tau[q, ind.all[,1]]*VE$tau[l, ind.all[,2]] # N-vector
            else
              tauql <- VE$tau[q, ind.all[,1]]*VE$tau[l, ind.all[,2]] + VE$tau[q, ind.all[,2]]*VE$tau[l, ind.all[,1]]
          }
          if (sparse){ # sparse
            rho.vec <- rep(1,N)
            rho.vec[no.obs] <- VE$rho[ind.ql]
            tauql_rho <- tauql*rho.vec
            wN_k <- colSums(Nijk.d*matrix(tauql_rho,N,2^d))  # N-vector
          } else{ # not sparse
            wN_k <- colSums(Nijk.d*matrix(tauql,N,2^d)) # N-vector
          }

          if (d==D) {
            pen_cste[ind.ql] <- 2*Dmax*max(wN_k)

            if (sparse){ # sparse
              sum_tauql <- sum(tauql)
              sum_rhotau_obs[ind.ql] <- sum(wN_k)
              sum_rhotau[ind.ql] <- sum(tauql_rho)
              if (sum_tauql>0)
                beta[ind.ql] <- sum_rhotau[ind.ql]/sum_tauql
            } else{ # non sparse
              sum_rhotau[ind.ql] <-  sum(tauql)
              sum_rhotau_obs[ind.ql] <- sum(wN_k)}
          }

          if (sum_rhotau[ind.ql]>0)
            AqlT[ind.ql] <- sum_rhotau_obs[ind.ql]/sum_rhotau[ind.ql]

          crit_pen_d = 2^d*(-sum(wN_k^2) + pen_cste[ind.ql])
          if (crit_pen_d<crit_pen[ind.ql]){
            best.d[ind.ql] <- 2^d
            crit_pen[ind.ql] <- crit_pen_d
            if (sum_rhotau[ind.ql]>0){
              loglam <- log(2^d/sum_rhotau[ind.ql]*wN_k/data$Time)
              loglam[loglam==-Inf] <- 0
              logintensities.ql[ind.ql,] <- rep(loglam,each=2^(D-d))
            }

          } # end case for criterion

        } # end for l
      }# end for q

      if (d>0)
        Nijk.d <- Nijk.d[,seq(1,2^d-1,by=2)]+Nijk.d[,seq(2,2^d,by=2)]
    } # end for d
  } # end for undirected case

  if (sparse){ # sparse parameters trimming - avoid numercial problems
    eps <- 1e-10
    beta[beta>(1-eps)] <- 1-eps
    beta[beta<eps] <- eps
  }

  return(list(beta=beta, sum_rhotau_obs=sum_rhotau_obs,
              sum_rhotau=sum_rhotau, AqlT = AqlT,
              logintensities.ql=logintensities.ql,best.d=best.d))
}


###################################################
## M-STEP FOR KERNEL ESTIMATOR
###################################################
# sparse method not implemented !!!


#' M step for kernel
#'
#' M step for kernel estimator
#'
#' @param data Data same of \link[ppsbm]{mainVEM}
#' @param VE Results of the previous VE for iterative computation
#' @param directed Boolean for directed (TRUE) or undirected (FALSE) case
#'
#' @references
#'
#' GRÉGOIRE , G. (1993). Least squares cross-validation for counting process intensities. Scand. J. Statist. 20, pp. 343–360.
#'
#' MATIAS, C., REBAFKA, T. & VILLERS, F. (2018).  A semiparametric extension of the stochastic block model for longitudinal networks. Biometrika.
#'
#' RAMLAU-HANSEN, H. (1983). Smoothing counting process intensities by means of kernel functions. Ann. Statist. 11, pp. 453–466.
#'
Mstep_kernel <- function(data,VE,directed){
  Q <- dim(VE$tau)[1]
  n <- dim(VE$tau)[2]

  N <- if (directed) n*(n-1) else n*(n-1)/2
  N_Q <- if (directed) Q^2 else Q*(Q+1)/2

  ### Note the different output
  ### logintensities.ql.ij = matrix(0,N_Q,N)  vs logintensities.ql = matrix(0,N_Q,Dmax)
  output <- list(sum_rhotau=rep(0,N_Q),sum_rhotau_obs=rep(0,N_Q),logintensities.ql.ij=matrix(0,N_Q,N),AqlT=rep(0,N_Q),beta=rep(0,N_Q))

  ind.all <- listNodePairs(n,directed)
  if (directed){ # directed
    ind <- 0
    for (q in 1:Q){
      for (l in 1:Q){
        ind <- ind + 1
        if (Q==1)
          tauql <- rep(1,N)
        else{
          tauql <- VE$tau[q, ind.all[,1]]*VE$tau[l, ind.all[,2]] # N-vector
        }
        sum_tauql <- sum(tauql)
        if (sum_tauql>0){
          output$sum_rhotau[ind] <- sum_tauql
          if (Q==1){
            poids <-rep(1/length(data$type.seq),length(data$type.seq))
            output$sum_rhotau_obs[ind] <- length(data$type.seq)
          }
          else {
            poids <- tauql[data$type.seq]
            output$sum_rhotau_obs[ind] <- sum(poids)
            if (is.finite(sum(poids/output$sum_rhotau_obs[ind]))){
              poids <- poids/output$sum_rhotau_obs[ind]
            }
            else poids <- rep(0,length(poids))
          }
          dens.estim <- density(data$time.seq,weights=poids,kernel="epanechnikov")
          for (ind.ij in 1:N){
            time.seq.ij <- data$time.seq[data$type.seq==ind.ij]
            if (length(time.seq.ij)>0){
              logintensities.obs <- log(output$sum_rhotau_obs[ind] /sum_tauql*approxfun(dens.estim$x,dens.estim$y)(time.seq.ij))  # values of log(intensities.ql) at the time points t_m
              output$logintensities.ql.ij[ind,ind.ij] <- sum(logintensities.obs[logintensities.obs>-Inf])
            }
          }
        }
      } # end loop for l
    } # end loop for q
  } else{ # undirected
    ind <- 0
    for (q in 1:Q){
      for (l in q:Q){
        ind <- ind + 1
        if (Q==1)
          tauql <- rep(1,N)
        if ((q==l) & (Q!=1))
          tauql <- VE$tau[q, ind.all[,1]]*VE$tau[l, ind.all[,2]]
        if ((q!=l) & (Q!=1))
          tauql <- VE$tau[q, ind.all[,1]]*VE$tau[l, ind.all[,2]] + VE$tau[q, ind.all[,2]]*VE$tau[l, ind.all[,1]]

        sum_tauql <- sum(tauql)
        if (sum_tauql>0){
          output$sum_rhotau[ind] <- sum_tauql
          if (Q==1){
            poids <-rep(1/length(data$type.seq),length(data$type.seq))
            output$sum_rhotau_obs[ind] <- length(data$type.seq)
          }
          else {
            poids <- tauql[data$type.seq]
            output$sum_rhotau_obs[ind] <- sum(poids)
            if (is.finite(sum(poids/output$sum_rhotau_obs[ind]))){
              poids <- poids/output$sum_rhotau_obs[ind]
            }
            else poids <- rep(0,length(poids))
          }
          dens.estim <- density(data$time.seq,weights=poids,kernel="epanechnikov")
          for (ind.ij in 1:N){
            time.seq.ij <- data$time.seq[data$type.seq==ind.ij]
            if (length(time.seq.ij)>0){
              logintensities.obs <- log(output$sum_rhotau_obs[ind] /sum_tauql*approxfun(dens.estim$x,dens.estim$y)(time.seq.ij))  # values of log(intensities.ql) at the time points t_m
              output$logintensities.ql.ij[ind,ind.ij] <- sum(logintensities.obs[logintensities.obs>-Inf])
            }
          }
        }
      } # end loop for l
    } # end loop for q
  } # end undirected case

  return(output)
}





###################################################
## EVALUATION OF CRITERION J
###################################################


#' Evaluation of criterion J
#'
#' Evaluation of the criterion J to verify the convergence of the VEM algorithm
#'
#' @param VE Results of the previous VE for iterative computation
#' @param mstep Results of the previous mstep for iterative computation
#'   \itemize{
#'     \item mstep$sum_rhotau : N_Q vector (not needed in the function)
#'     \item mstep$sum_rhotau_obs : N_Q vector
#'     \item mstep$logintensities.ql :  N_Q x Dmax matrix
#'     \item m.step$beta : N_Q vector
#'   }
#' @param data Data same of \link[ppsbm]{mainVEM}
#' @param directed Boolean for directed (TRUE) or undirected (FALSE) case
#' @param sparse Boolean for sparse (TRUE) or not sparse (FALSE) case
#' @param method List of string. Can be "hist" for histogram method or "kernel" for kernel method
#'
JEvalMstep  <- function(VE, mstep,data,directed,sparse,method='hist'){
  Q <- dim(VE$tau)[1]
  n <- dim(VE$tau)[2]

  if (method=='hist'){
    Dmax <- dim(data$Nijk)[2]
  }

  N <- if (directed) n*(n-1) else n*(n-1)/2

  J <- -sum(mstep$sum_rhotau_obs)

  if (Q==1){
    if (method=='hist'){
      logintensities.ql.ij <- rowSums(matrix(mstep$logintensities.ql,N,Dmax,byrow=TRUE)*data$Nijk)  # N - vector
    } else {logintensities.ql.ij <- mstep$logintensities.ql.ij[1,]}
    if (sparse) {
      nb.no.obs <- sum((rowSums(data$Nijk)==0))
      nb.obs <- N
    }
    J <- J + sum(logintensities.ql.ij)
  }
  else{
    # Q!=1
    ind.all <- listNodePairs(n,directed)
    ind.ql <- 0
    N_Q <- if (directed) Q^2 else Q*(Q+1)/2
    nb.no.obs <- rep(NA,N_Q)
    nb.obs <- rep(NA,N_Q)

    for (q in 1:Q){
      ind.l <- if (directed) 1:Q else q:Q
      for (l in ind.l){
        ind.ql <- ind.ql+1
        if (directed) { tauql <- VE$tau[q, ind.all[,1]]*VE$tau[l, ind.all[,2]] }
        else {
          if ((q==l))
            tauql <- VE$tau[q, ind.all[,1]]*VE$tau[l, ind.all[,2]] # N-vector
          else
            tauql <- VE$tau[q, ind.all[,1]]*VE$tau[l, ind.all[,2]] + VE$tau[q, ind.all[,2]]*VE$tau[l, ind.all[,1]]
        }# end else undirected

        if (method=='hist'){
          logintensities.ql.ij <- rowSums(matrix(mstep$logintensities.ql[ind.ql,],N,Dmax,byrow=TRUE)*data$Nijk)  # N - vector
        } else {logintensities.ql.ij <- mstep$logintensities.ql.ij[ind.ql,]}
        J <- J + sum(tauql*logintensities.ql.ij)
        # Sparse case
        if (sparse) {
          nb.no.obs[ind.ql] <- sum(tauql * (rowSums(data$Nijk)==0))
          nb.obs[ind.ql] <- sum(tauql)
        }
      } # end loop for l
    } # end loop for q
  }# end of Q!=1

  log.tau <- log(VE$tau)
  log.tau[VE$tau==0] <- 0
  pi <- rowMeans(VE$tau)
  log.pi  <- log(pi)
  log.pi[pi==0]   <- 0

  if (sparse){ # sparse case
    log.rho <- log(VE$rho)
    log.rho[VE$rho==0] <- 0
    log.1mrho <- log(1-VE$rho)
    log.1mrho[VE$rho==1] <- 0
    log.beta <- log(mstep$beta)
    log.beta[mstep$beta==0] <- 0
    log.1mbeta <- log(1-mstep$beta)
    log.1mbeta[mstep$beta==1] <- 0
    term.beta.sparse <- sum ( log.beta *(nb.no.obs*(VE$rho -1) + nb.obs) +log.1mbeta * (nb.no.obs *(1-VE$rho)))
    psi <- VE$rho*log.rho + (1-VE$rho*log.1mrho)
    term.sparse <- sum(nb.no.obs* psi)
  } else {
    term.beta.sparse <- 0
    term.sparse <- 0
  } # no entropy for sparse parameter in non sparse case

  compl.log.lik <- J + sum ( log.pi%*% VE$tau ) +term.beta.sparse
  J <- J - sum( VE$tau * log.tau ) + sum ( log.pi%*% VE$tau ) +term.beta.sparse - term.sparse

  return(list(J=J,compl.log.lik=compl.log.lik))
}


