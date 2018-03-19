###################################################
## MODEL SELECTION BY Integrated Classification Likelihood criterion  &
## ANALYSIS OF RESULTS
###################################################



#' Selects the number of groups with ICL
#'
#' Selects the number of groups with Integrated Classification Likelihood Criterion
#'
#' @param data List with 2 components:
#'   \itemize{
#'     \item $Time - [0,data$Time] is the total time interval of observation
#'     \item $Nijk - data matrix with the statistics per process \eqn{N_{ij}} and sub-intervals \eqn{k}
#'   }
#' @param n Total number of nodes \eqn{n}
#' @param Qmin Minimum number of groups
#' @param Qmax Maximum number of groups
#' @param directed Boolean for directed (TRUE) or undirected (FALSE) case
#' @param sparse Boolean for sparse (TRUE) or not sparse (FALSE) case
#' @param sol.hist.sauv List of size Qmax-Qmin+1 obtained from running mainVEM(data,n,Qmin,Qmax,method='hist')
#'
#' @export
#'
#' @references
#'
#' BIERNACKI, C., CELEUX, G. & GOVAERT, G. (2000). Assessing a mixture model for clustering with the integrated completed likelihood. IEEE Trans. Pattern Anal. Machine Intel. 22, 719–725.
#'
#' CORNELI, M., LATOUCHE, P. & ROSSI, F. (2016). Exact ICL maximization in a non-stationary temporal extension of the stochastic block model for dynamic networks. Neurocomputing 192, 81 – 91.
#'
#' DAUDIN, J.-J., PICARD, F. & ROBIN, S. (2008). A mixture model for random graphs. Statist. Comput. 18, 173–183.
#'
#' MATIAS, C., REBAFKA, T. & VILLERS, F. (2018).  A semiparametric extension of the stochastic block model for longitudinal networks. Biometrika.
#'
#' @examples
#' # load data of a synthetic graph with 50 individuals and 3 clusters
#' n <- 50
#'
#' # compute data matrix with precision d_max=3
#' Dmax <- 2^3
#' data <- list(Nijk=statistics(generated_Q3$data,n,Dmax,directed=FALSE),
#'     Time=generated_Q3$data$Time)
#'
#' # ICL-model selection
#' sol.selec_Q <- modelSelection_Q(data,n,Qmin=1,Qmax=4,directed=FALSE,
#'     sparse=FALSE,generated_sol_hist)
#'
#' # best number Q of clusters:
#' sol.selec_Q$Qbest
#'
modelSelection_Q <- function(data,
                             n,
                             Qmin=1,
                             Qmax,
                             directed=TRUE,
                             sparse=FALSE,
                             sol.hist.sauv){
  N <- if (directed) n*(n-1) else n*(n-1)/2

  ICLbest <- -Inf
  sol.Qbest <- NA
  all.ICL <- rep(NA,Qmax-Qmin+1)
  all.compl.log.likelihood <- rep(NA,Qmax-Qmin+1)
  all.J <- rep(NA,Qmax-Qmin+1)
  all.pen  <- rep(NA,Qmax-Qmin+1)
  for (Qcand in Qmin:Qmax){
    cat("----Qcand----",Qcand,"\n")
    sol.hist <-sol.hist.sauv[[Qcand]]
    all.J[Qcand-Qmin+1]<-sol.hist$J

    # arguments for complete.log.lik
    VE <- list(tau=sol.hist$tau,rho=sol.hist$rho)
    no.obs <- (rowSums(data$Nijk)==0)
    Nij <- rowSums(data$Nijk)
    # compute sum_rho_obs
    N_Q <- if (directed) Qcand^2 else Qcand*(Qcand+1)/2
    sum_rhotau_obs <- rep(0,N_Q)
    if (sparse){
      ind.all <- listNodePairs(n,directed)
      ind.ql <- 0
      for (q in 1:Qcand){
        ind.l <- if (directed) 1:Qcand else q:Qcand
        for (l in ind.l){
          ind.ql <- ind.ql + 1
          if (Qcand==1)	{
            tauql <- rep(1,N)
          }
          else {
            if (directed){ # directed case
              tauql <- sol.hist$tau[q, ind.all[,1]]*sol.hist$tau[l, ind.all[,2]] # N-vector
            } else { # undirected case
              if ((q==l))
                tauql <- sol.hist$tau[q, ind.all[,1]]*sol.hist$tau[l, ind.all[,2]] # N-vector
              else
                tauql <- sol.hist$tau[q, ind.all[,1]]*sol.hist$tau[l, ind.all[,2]] + sol.hist$tau[q, ind.all[,2]]*sol.hist$tau[l, ind.all[,1]]
            }
          } # end else Q!=1
          rhoql.vec <- rep(1,N)
          rhoql.vec[no.obs] <- sol.hist$rho[ind.ql]
          sum_rhotau_obs[ind.ql] <- sum(tauql * rhoql.vec * Nij)
        } # end for l
      } # end for q
    } # end if sparse

    mstep <- list(logintensities.ql=sol.hist$logintensities.ql,sum_rhotau_obs=sum_rhotau_obs,beta=sol.hist$beta)
    all.compl.log.likelihood[Qcand-Qmin+1] <- JEvalMstep(VE,mstep,data,directed,sparse)$compl.log.lik

    # penalty term
    sum.d.ql <- sum(sol.hist$best.d)
    all.pen[Qcand-Qmin+1] <- (Qcand-1)/2*log(n) + sum.d.ql/2*log(N)
    if (directed*sparse){ # penalizing beta values
      all.pen[Qcand-Qmin+1] <- all.pen[Qcand-Qmin+1]  + Qcand^2/2*log(N)
    }
    if ((!directed)*sparse){ # penalizing beta values
      all.pen[Qcand-Qmin+1] <- all.pen[Qcand-Qmin+1]  + Qcand*(Qcand+1)/4*log(N)}

    ICL <- all.compl.log.likelihood[Qcand-Qmin+1] - all.pen[Qcand-Qmin+1]
    all.ICL[Qcand-Qmin+1] <- ICL
    if (ICL>ICLbest){
      ICLbest <- ICL
      Qbest <- Qcand
      sol.Qbest <- sol.hist
    }
  } # end loop for Qcand
  return(list(Qbest=Qbest, sol.Qbest=sol.Qbest, Qmin=Qmin,all.J=all.J,all.ICL=all.ICL,all.compl.log.likelihood=all.compl.log.likelihood,all.pen=all.pen))
}




#' Plots for model selection
#'
#' @param model.selec_Q Output from modelSelection_Q()
#'
#' @export
#'
#' @examples
#' # load data of a synthetic graph with 50 individuals and 3 clusters
#' n <- 50
#'
#' # compute data matrix with precision d_max=3
#' Dmax <- 2^3
#' data <- list(Nijk=statistics(generated_Q3$data,n,Dmax,directed=FALSE),
#'     Time=generated_Q3$data$Time)
#'
#' # ICL-model selection
#' sol.selec_Q <- modelSelection_Q(data,n,Qmin=1,Qmax=4,directed=FALSE,
#'     sparse=FALSE,generated_sol_hist)
#'
#' # plot ICL
#' modelSelec_QPlot(sol.selec_Q)
#'
modelSelec_QPlot <- function(model.selec_Q){
  Qmin <- model.selec_Q$Qmin
  Qmax <- length(model.selec_Q$all.ICL)+Qmin-1
  par(mfrow=c(1,2))
  ymin <- min(model.selec_Q$all.compl.log.likelihood,model.selec_Q$all.ICL)
  ymax <- max(model.selec_Q$all.compl.log.likelihood,model.selec_Q$all.ICL)
  plot(Qmin:Qmax,model.selec_Q$all.compl.log.likelihood,type='l',col='blue',xlab=paste('Qbest =',model.selec_Q$Qbest),ylab='ICL (red), complete likelihood (blue)',ylim=c(ymin,ymax))
  lines(Qmin:Qmax,model.selec_Q$all.ICL,col='red')
  plot(Qmin:Qmax,model.selec_Q$all.J,ylab='J criterion',type='l',xlab='Q')
  return()
}


#######################################
##  Bootstrap and Confidence Interval
#######################################

#' Bootstrap and Confidence Interval
#'
#' Not for sparse models and only for histograms
#'
#' @param sol sol
#' @param Time time
#' @param R Number of bootstrap samples
#' @param alpha Level of confidence : \eqn{1- \alpha}
#' @param nbcores Number of cores for parallel execution
#'
#' If set to 1 it does sequential execution
#'
#' Beware: parallelization with fork (multicore) : doesn't work on Windows!
#' @param d_part Maximal level for finest partitions of time interval [0,T], used for kmeans initializations.
#'   \itemize{
#'     \item Algorithm takes partition up to depth \eqn{2^d} with \eqn{d=1,...,d_{part}}
#'     \item Explore partitions \eqn{[0,T], [0,T/2], [T/2,T], ... [0,T/2^d], ...[(2^d-1)T/2^d,T]}
#'     \item Total number of partitions \eqn{npart= 2^{(d_{part} +1)} - 1}
#'   }
#' @param n_perturb Number of different perturbations on k-means result
#' @param perc_perturb Percentage of labels that are to be perturbed (= randomly switched)
#' @param directed Boolean for directed (TRUE) or undirected (FALSE) case
#' @param filename filename
#'
#' @export
#'
#' @examples
#'
#' # data of a synthetic graph with 50 individuals and 3 clusters
#'
#' n <- 50
#' Q <- 3
#'
#' Time <- generated_Q3$data$Time
#' data <- generated_Q3$data
#' z <- generated_Q3$z
#'
#' Dmax <- 2^3
#'
#' # VEM-algo hist
#' sol.hist <- mainVEM(list(Nijk=statistics(data,n,Dmax,directed=FALSE),Time=Time),
#'      n,Qmin=3,directed=FALSE,method='hist',d_part=1,n_perturb=0)[[1]]
#'
#' # compute bootstrap confidence bands
#' boot <- bootstrap_and_CI(sol.hist,Time,R=10,alpha=0.1,nbcores=1,d_part=1,n_perturb=0,
#'      directed=FALSE)
#'
#' # plot confidence bands
#' alpha.hat <- exp(sol.hist$logintensities.ql)
#' vec.x <- (0:Dmax)*Time/Dmax
#' ind.ql <- 0
#' par(mfrow=c(2,3))
#' for (q in 1:Q){
#'   for (l in q:Q){
#'     ind.ql <- ind.ql+1
#'     ymax <- max(c(boot$CI.limits[ind.ql,2,],alpha.hat[ind.ql,]))
#'     plot(vec.x,c(alpha.hat[ind.ql,],alpha.hat[ind.ql,Dmax]),type='s',col='black',
#'         ylab='Intensity',xaxt='n',xlab= paste('(',q,',',l,')',sep=""),
#'         cex.axis=1.5,cex.lab=1.5,ylim=c(0,ymax),main='Confidence bands')
#'     lines(vec.x,c(boot$CI.limits[ind.ql,1,],boot$CI.limits[ind.ql,1,Dmax]),col='blue',
#'         type='s',lty=3)
#'     lines(vec.x,c(boot$CI.limits[ind.ql,2,],boot$CI.limits[ind.ql,2,Dmax]),col='blue',
#'         type='s',lty=3)
#'   }
#' }
#'
bootstrap_and_CI <- function(sol,Time,R,alpha=0.05,nbcores=1,d_part=5,n_perturb=10,perc_perturb=.2,directed,filename=NULL){
  Dmax <- ncol(sol$logintensities.ql)
  Q <- nrow(sol$tau)
  n <- ncol(sol$tau)
  N_Q <- if (directed)  Q^2 else Q*(Q+1)/2
  boot.sol <- array(NA,dim=c(R,N_Q,Dmax)) # storage of logintensities.ql values
  all.intens <- exp(sol$logintensities.ql)
  all.intens[sol$logintensities.ql==0] <- 0
  prop.groups=rowMeans(sol$tau)
  for (i in 1:R){
    cat('R:',i,'\n')
    data <- generateDynppsbmConst(all.intens,Time,n,prop.groups,directed)
    Nijk <- statistics(data$data,n,Dmax,directed)
    obs <- list(Nijk=Nijk,Time=Time)
    sol <- mainVEM(obs,n,Qmin=Q,Qmax=Q,directed,method='hist',cores=nbcores,d_part=d_part,n_perturb=n_perturb,perc_perturb=perc_perturb)[[1]]

    # label switching
    boot.sol[i,,] <- sortIntensities(sol$logintensities.ql,data$z,sol$tau,directed)
    if (!is.null(filename))
      save(boot.sol,file=filename)
  }

  CIlim <- confidenceInterval(boot.sol,alpha)
  return(list(CI.limits=CIlim, boot.sol=boot.sol))
}




#' Confidence Interval
#'
#' Compute confidence bands for all pair of groups \eqn{(q,l)}
#'
#' @param boot.sol Bootstrap list of estimators
#' @param alpha Level of confidence : 1 - \eqn{\alpha}
#'
#' @export
#'
#' @examples
#'
#' # data of a synthetic graph with 50 individuals and 3 clusters
#'
#' n <- 50
#' Q <- 3
#'
#' Time <- generated_Q3$data$Time
#' data <- generated_Q3$data
#' z <- generated_Q3$z
#'
#' Dmax <- 2^3
#'
#' # VEM-algo hist
#' sol.hist <- mainVEM(list(Nijk=statistics(data,n,Dmax,directed=FALSE),Time=Time),
#'      n,Qmin=3,directed=FALSE,method='hist',d_part=1,n_perturb=0)[[1]]
#'
#' # compute bootstrap confidence bands
#' boot <- bootstrap_and_CI(sol.hist,Time,R=5,alpha=0.1,nbcores=1,d_part=1,n_perturb=0,
#'      directed=FALSE)
#'
#' boot.sol <- boot$boot.sol
#'
#' confidenceInterval(boot.sol)
#'
confidenceInterval <- function(boot.sol,alpha=0.05){
  dim.sol <- dim(boot.sol)
  R <- dim.sol[1]
  N_Q <- dim.sol[2]
  Dmax <- dim.sol[3]
  CIlim <- array(NA,dim=c(N_Q,3,Dmax))
  for (ql in 1:N_Q)
    CIlim[ql,,] <- apply(boot.sol[,ql,],2,function(x){quantile(x,probs=c(alpha/2,1-alpha/2,.5))})
  CIlim <- exp(CIlim)
  return(CIlim)
}

############################################################
############################################################


#' Direct kernel estimator intensities
#'
#' Compute smooth intensities with direct kernel estimation of intensities relying on a classification tau.
#' This can be used with the values \eqn{\tau} obtained on a dataset with mainVEM function run with 'hist' method.
#'
#' Warning : sparse case not implemented !!!
#'
#' @param data List with 3 components:
#'   \itemize{
#'     \item data$time.seq : sequence of observed time points of the m-th event (M-vector)
#'     \item data$type.seq : sequence of observed values convertNodePair(i,j,n,directed) (auxiliary.R)
#'         of process that produced the mth event (M-vector)
#'     \item $Time - [0,data$Time] is the total time interval of observation
#'   }
#' @param tau \eqn{\tau}
#' @param Q Total number of groups
#' @param n Total number of nodes
#' @param directed Boolean for directed (TRUE) or undirected (FALSE) case
#' @param rho \eqn{\rho}
#' @param sparse Boolean for sparse (TRUE) or not sparse (FALSE) case
#' @param nb.points Number of points
#'
#' @export
#'
#' @examples
#'
#' # The generated_sol_kernel was generated calling mainVEM with kernel method on the generated_Q3 data
#' # (50 individuals and 3 clusters)
#'
#' data <- generated_Q3$data
#'
#' n <- 50
#' Q <- 3
#'
#'
#' # compute smooth intensity estimators
#' sol.kernel.intensities <- kernelIntensities(data,generated_sol_kernel$tau,Q,n,directed=FALSE)
#'
kernelIntensities <- function(data,tau,Q,n,directed,rho=1,sparse=FALSE,nb.points=1000){
  N <- if (directed) n*(n-1) else n*(n-1)/2
  N_Q <- if (directed) Q^2 else Q*(Q+1)/2

  if (sparse){#
    stop("Sorry sparse method not implemented yet")
  }

  K = 1
  step <- data$Time/nb.points
  partition <- seq(0,data$Time, by=step)
  hat.intensities.ql <- matrix(0,N_Q,length(partition))
  Nijk <- statistics(data, n, K, directed)
  no.obs <- (Nijk==0)

  if (directed){ # directed case
    ind.all <- listNodePairs(n,directed)
    ind <- 0
    for (q in 1:Q){
      for (l in 1:Q){
        ind <- ind + 1
        if (Q==1){
          rhotauql <- rep(1,N)
          if (sparse){rhotauql[no.obs] <- rho}
        } else {
          rhotauql <- tau[q, ind.all[,1]]*tau[l, ind.all[,2]]  # N vector
          if (sparse){rhotauql[no.obs] <- rho[convertGroupPair(q,l,Q,directed)]*rhotauql[no.obs]}
        }
        Yql <- sum(rhotauql)
        if (Yql>0){
          if (Q==1){
            poids <-rep(1/length(data$type.seq),length(data$type.seq))
            sum_tauql_obs <- length(data$type.seq)
          }
          else {
            rhotau_obs <- rhotauql[data$type.seq]
            rhotau_obs_sum <- sum(rhotau_obs)
            #
            if (is.finite(sum(rhotau_obs/rhotau_obs_sum))){
              poids <- rhotau_obs/rhotau_obs_sum
            }
            else poids <- rep(0,length(rhotau_obs))
            dens.estim <- density(data$time.seq,weights=poids,kernel="epanechnikov")
            hat.intensities.ql[ind,] <- rhotau_obs_sum /Yql*(approxfun(dens.estim$x,dens.estim$y)(partition))
          }
        }
      } # end loop for l
    } # end loop for q
  } else { # undirected case
    ind.all <- listNodePairs(n,directed)
    ind <- 0
    for (q in 1:Q){
      for (l in q:Q){
        ind <- ind + 1
        if (Q==1){
          rhotauql <- rep(1,N)
          rhotauql[no.obs] <- rho
        } else {
          if ((q==l)){
            rhotauql <- tau[q, ind.all[,1]]*tau[l, ind.all[,2]]  # N vector
            if (sparse){rhotauql[no.obs] <- rho[convertGroupPair(q,l,Q,directed)]*rhotauql[no.obs]}
          } else {
            rhotauql <- tau[q, ind.all[,1]]*tau[l, ind.all[,2]] + tau[q, ind.all[,2]]*tau[l, ind.all[,1]]
            if (sparse){rhotauql[no.obs] <- rho[convertGroupPair(q,l,Q,directed)]*rhotauql[no.obs]}
          }
        }
        Yql <- sum(rhotauql)
        if (Yql>0){
          if (Q==1){
            poids <-rep(1/length(data$type.seq),length(data$type.seq))
            sum_tauql_obs <- length(data$type.seq)
          }
          else {
            rhotau_obs <- rhotauql[data$type.seq]
            rhotau_obs_sum <- sum(rhotau_obs)
            if (is.finite(sum(rhotau_obs/rhotau_obs_sum))){
              poids <- rhotau_obs/rhotau_obs_sum
            }
            else poids <- rep(0,length(rhotau_obs))
            dens.estim <- density(data$time.seq,weights=poids,kernel="epanechnikov")
            hat.intensities.ql[ind,] <- rhotau_obs_sum /Yql*(approxfun(dens.estim$x,dens.estim$y)(partition))
          }
        }
      } # end loop for l
    } # end loop for q
  }
  return(hat.intensities.ql)
}



#' Adjusted Rand Index (ARI)
#'
#' Compute the Adjusted Rand Index (ARI) between the true latent variables and the estimated latent variables
#'
#' @param z Matrix of size  \eqn{Q \times n} with entries = 0 or 1 : 'true' latent variables
#' @param hat.z Matrix of \eqn{Q \times n}  with 0<entries<1 : estimated latent variables
#'
#' @export
#'
#' @examples
#'
#' z <- matrix(c(1,1,0,0,0,0, 0,0,1,1,0,0, 0,0,0,0,1,1), nrow = 3, byrow = TRUE)
#' hat.z <- matrix(c(0,0,1,1,0,0, 1,1,0,0,0,0, 0,0,0,0,1,1), nrow = 3, byrow = TRUE)
#'
#' ARI(z, hat.z)
#'
ARI <- function(z, hat.z){
  clustering.z <- colMaxs(z, FALSE) # apply(z,2,which.max)
  clustering.hat.z <- colMaxs(hat.z, FALSE) # apply(hat.z,2,which.max)

  # first, get crosstabs
  ctab <- table(clustering.z,clustering.hat.z);

  # now calculate 4 intermediary sums
  cellsum <- sum(ctab*(ctab-1)/2)
  totsum <- sum(ctab)*(sum(ctab)-1)/2

  # use matrix multiplication to get row and column marginal sums
  rows <- ctab %*% rep(1,ncol(ctab))
  rowsum <- sum(rows*(rows-1)/2)
  cols <- rep(1,nrow(ctab)) %*% ctab
  colsum <- sum(cols*(cols-1)/2)

  # now put them together
  adj.rand <- (cellsum - (rowsum*colsum/totsum))/(.5*(rowsum +colsum)-(rowsum*colsum/totsum))
  return (adj.rand)
}

###################################################
##   COMPARISON OF THE Z ESTIMATORS WITH TRUE VALUES
###################################################

#' Optimal matching between 2 clusterings
#'
#' Compute the permutation of the rows of hat.z that has to be applied to obtain the "same order" as z.
#' Compute optimal matching between 2 clusterings using Hungarian algorithm
#'
#'
#' @param z Matrice of size  \eqn{Q \times n}
#' @param hat.z Matrice of size  \eqn{Q \times n}
#'
#' @export
#'
#' @references
#'
#' HUBERT, L. & ARABIE, P. (1985). Comparing partitions. J. Classif. 2, 193–218.
#'
#' MATIAS, C., REBAFKA, T. & VILLERS, F. (2018).  A semiparametric extension of the stochastic block model for longitudinal networks. Biometrika.
#'
#' @examples
#'
#' z <- matrix(c(1,1,0,0,0,0, 0,0,1,1,0,0, 0,0,0,0,1,1), nrow = 3, byrow = TRUE)
#' hat.z <- matrix(c(0,0,1,1,0,0, 1,1,0,0,0,0, 0,0,0,0,1,1), nrow = 3, byrow = TRUE)
#'
#' perm <- permuteZEst(z,hat.z)
#'
permuteZEst <- function(z,hat.z){
  map.z <- colMaxs(z, FALSE) # apply(z,2,which.max)
  map.hatz <- colMaxs(hat.z, FALSE) # apply(hat.z,2,which.max)
  tab <- table(map.z,map.hatz)
  if (ncol(tab) < nrow(tab) ){
    res <- solve_LSAP(t(tab),maximum = TRUE)
  } else {res <- solve_LSAP(tab,maximum = TRUE)}
  return(res)
}


#' Sort intensities
#'
#' Sort intensities associated with hat.z "in the same way" as the original intensities associated with z by permutation of rows
#'
#' @param intensities Intensities \eqn{\alpha}
#' @param z Matrice of size  \eqn{Q \times n}
#' @param hat.z Matrice of size  \eqn{Q \times n}
#' @param directed Boolean for directed (TRUE) or undirected (FALSE) case
#'
#' @export
#'
#' @references
#'
#' HUBERT, L. & ARABIE, P. (1985). Comparing partitions. J. Classif. 2, 193–218.
#'
#' MATIAS, C., REBAFKA, T. & VILLERS, F. (2018).  A semiparametric extension of the stochastic block model for longitudinal networks. Biometrika.
#'
#' @examples
#'
#' z <- matrix(c(1,1,0,0,0,0, 0,0,1,1,0,0, 0,0,0,0,1,1), nrow = 3, byrow = TRUE)
#' hat.z <- matrix(c(0,0,1,1,0,0, 1,1,0,0,0,0, 0,0,0,0,1,1), nrow = 3, byrow = TRUE)
#'
#' intens <- matrix(c(1,1,1,2,2,2,3,3,3),9)
#'
#' sortIntensities(intens,z,hat.z, TRUE)
#'
sortIntensities <- function(intensities,z,hat.z,directed){
  Q <- nrow(z)
  perm <- permuteZEst(z,hat.z)

  # permutation of rows of intensities
  sol <- intensities
  for (q in 1:Q){
    vec.l <- if (directed) 1:Q else q:Q
    for (l in vec.l)
      if (!is.na(perm[q]*perm[l]))
        sol[convertGroupPair(q,l,Q,directed),]<- intensities[convertGroupPair(perm[q],perm[l],Q,directed),]
  }
  return(sol)
}



