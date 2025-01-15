###################################################
## MODEL SELECTION BY Integrated Classification Likelihood criterion  &
## ANALYSIS OF RESULTS
###################################################



#' Selects the number of groups with ICL criterion
#'
#' Selects the number of groups with Integrated Classification Likelihood (ICL) criterion.
#'
#' @param data List with 2 components:
#'   \itemize{
#'     \item \code{Time} - Positive real number. [0,Time] is the total time interval of observation.
#'     \item \code{Nijk} - Data matrix with the statistics per process \eqn{N_{ij}} and sub-intervals \eqn{1\le k\le K}.
#'   }
#' @param n Total number of nodes,  \eqn{1\le i \le n}.
#' @param Qmin Minimum number of groups.
#' @param Qmax Maximum number of groups.
#' @param directed Boolean for directed (TRUE) or undirected (FALSE) case.
#' @param sparse Boolean for sparse (TRUE) or not sparse (FALSE) case.
#' @param sol.hist.sauv List of size Qmax-Qmin+1 obtained from running \link[ppsbm]{mainVEM} on the data with method='hist'.
#'
#' @return The function outputs a list of 7 components:
#' \itemize{
#'  \item {\code{Qbest}} Selected value of the number of groups in [Qmin, Qmax].
#'  \item {\code{sol.Qbest}} Solution of the \link[ppsbm]{mainVEM} function for the number of groups Qbest.
#'  \item {\code{Qmin}} Minimum number of groups used.
#'  \item {\code{all.J}} Vector of length Qmax-Qmin+1. Each value is the estimated ELBO function \eqn{J} for estimation with \eqn{Q} groups, \eqn{Qmin \le Q \le Qmax}.
#'  \item {\code{all.ICL}} Vector of length Qmax-Qmin+1. Each value is the ICL value for estimation with \eqn{Q} groups, \eqn{Qmin \le Q \le Qmax}.
#'  \item {\code{all.compl.log.likelihood}} Vector of length Qmax-Qmin+1. Each value is the estimated complete log-likelihood value for estimation with \eqn{Q} groups, \eqn{Qmin \le Q \le Qmax}.
#'  \item {\code{all.pen}} Vector of length Qmax-Qmin+1. Each value is the penalty term in ICL for estimation with \eqn{Q} groups, \eqn{Qmin \le Q \le Qmax}.
#' }
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
#' MATIAS, C., REBAFKA, T. & VILLERS, F. (2018).  A semiparametric extension of the stochastic block model for longitudinal networks. Biometrika. 105(3): 665-680.
#'
#' @examples
#' # load data of a synthetic graph with 50 individuals and 3 clusters
#' n <- 50
#'
#' # compute data matrix of counts per subinterval with precision d_max=3
#' # (ie nb of parts K=2^{d_max}=8).
#' K <- 2^3
#' data <- list(Nijk=statistics(generated_Q3$data,n,K,directed=FALSE),
#'     Time=generated_Q3$data$Time)
#'
#' # ICL-model selection with groups ranging from 1 to 4
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
    sol.hist <-sol.hist.sauv[[Qcand-Qmin+1]]
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
#' Plots the Integrated Classification Likelihood (ICL) criterion, the Complete Log-Likelihood (CLL) and the ELBO (J criterion).
#'
#' @param model.selec_Q Output from \link[ppsbm]{modelSelection_Q}.
#'
#' @export
#'
#' @examples
#' # load data of a synthetic graph with 50 individuals and 3 clusters
#' n <- 50
#'
#' # compute data matrix of counts per subinterval with precision d_max=3 :
#' # (ie nb of parts K=2^{d_max}=8).
#' K <- 2^3
#' data <- list(Nijk=statistics(generated_Q3$data,n,K,directed=FALSE),
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
  plot(Qmin:Qmax,model.selec_Q$all.compl.log.likelihood,type='l',col='blue',xlab=paste('Qbest =',model.selec_Q$Qbest),ylab='ICL (red), CLL (blue)',ylim=c(ymin,ymax))
  lines(Qmin:Qmax,model.selec_Q$all.ICL,col='red')
  plot(Qmin:Qmax,model.selec_Q$all.J,ylab='ELBO (J criterion)',type='l',xlab='Q')
  return()
}


#######################################
##  Bootstrap and Confidence Interval
#######################################

#' Bootstrap and Confidence Bands
#'
#' Plots confidence bands for estimated intensities between pairs of groups obtained by bootstrap.
#'
#' Not for sparse models and only for histogram method.
#'
#' @param sol One list (for one value of \eqn{Q}) output by \link[ppsbm]{mainVEM} with \code{hist} method.
#' @param Time Positive real number. [0,Time] is the total time interval of observation.
#' @param R Number of bootstrap samples.
#' @param alpha Level of confidence: \eqn{1- \alpha}.
#' @param nbcores Number of cores for parallel execution.
#'
#' If set to 1 it does sequential execution.
#'
#' Beware: parallelization with fork (multicore): doesn't work on Windows!
#' @param d_part Maximal level for finest partitions of time interval [0,T], used for kmeans initializations on the bootstrap samples
#'   \itemize{
#'     \item Algorithm takes partition up to depth \eqn{2^d} with \eqn{d=1,...,d_{part}}
#'     \item Explore partitions \eqn{[0,T], [0,T/2], [T/2,T], ... [0,T/2^d], ...[(2^d-1)T/2^d,T]}
#'     \item Total number of partitions \eqn{npart= 2^{(d_{part} +1)} - 1}
#'   }
#' @param n_perturb Number of different perturbations on k-means result on the bootstrap samples.
#' @param perc_perturb Percentage of labels that are to be perturbed (= randomly switched)  on the bootstrap samples.
#' @param directed Boolean for directed (TRUE) or undirected (FALSE) case.
#' @param filename Name of the file where to save the results.
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
#' K <- 2^3
#'
#' # VEM-algo hist:
#' sol.hist <- mainVEM(list(Nijk=statistics(data,n,K,directed=FALSE),Time=Time),
#' n,Qmin=3,directed=FALSE,method='hist',d_part=1,n_perturb=0)[[1]]
#'
#' # compute bootstrap confidence bands
#' boot <- bootstrap_and_CI(sol.hist,Time,R=10,alpha=0.1,nbcores=1,d_part=1,n_perturb=0,
#'      directed=FALSE)
#'
#' # plot confidence bands
#' alpha.hat <- exp(sol.hist$logintensities.ql)
#' vec.x <- (0:K)*Time/K
#' ind.ql <- 0
#' par(mfrow=c(2,3))
#' for (q in 1:Q){
#'   for (l in q:Q){
#'     ind.ql <- ind.ql+1
#'     ymax <- max(c(boot$CI.limits[ind.ql,2,],alpha.hat[ind.ql,]))
#'     plot(vec.x,c(alpha.hat[ind.ql,],alpha.hat[ind.ql,K]),type='s',col='black',
#'         ylab='Intensity',xaxt='n',xlab= paste('(',q,',',l,')',sep=""),
#'         cex.axis=1.5,cex.lab=1.5,ylim=c(0,ymax),main='Confidence bands')
#'     lines(vec.x,c(boot$CI.limits[ind.ql,1,],boot$CI.limits[ind.ql,1,K]),col='blue',
#'         type='s',lty=3)
#'     lines(vec.x,c(boot$CI.limits[ind.ql,2,],boot$CI.limits[ind.ql,2,K]),col='blue',
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
#' @keywords internal
#'
#' @return Array with 3 dimensions and size \eqn{Q(Q+1)/2\times 3\times K} (if undirected) or \eqn{Q^2\times 3\times K} (when undirected)
#' containing for each pair of groups \eqn{(q,l)} (first dimension) and each \eqn{k}-th subinterval (third dimension) the 3 quantiles at levels \eqn{(\alpha/2,1-\alpha/2,.5)} (second dimension).
#'
# # Example
# # data of a synthetic graph with 50 individuals and 3 clusters
#
# n <- 50
# Q <- 3
#
# Time <- generated_Q3$data$Time
# data <- generated_Q3$data
# z <- generated_Q3$z
#
# K <- 2^3
#
# # VEM-algo hist
# sol.hist <- mainVEM(list(Nijk=statistics(data,n,K,directed=FALSE),Time=Time),
#     n,Qmin=3,directed=FALSE,method='hist',d_part=1,n_perturb=0)[[1]]
#
# # compute bootstrap confidence bands
# boot <- bootstrap_and_CI(sol.hist,Time,R=5,alpha=0.1,nbcores=1,d_part=1,n_perturb=0,
#      directed=FALSE)
#
# boot.sol <- boot$boot.sol
#
# confidenceInterval(boot.sol)
#
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
#' Compute smooth intensities with direct kernel estimation of intensities relying on a classification \eqn{\tau}.
#' This can be used with the values \eqn{\tau} obtained on a dataset with \link[ppsbm]{mainVEM} function.
#'
#' Warning: sparse case not implemented !!!
#'
#' @param data List with 3 components:
#'   \itemize{
#'     \item \code{time.seq} - Vector of observed time points of the events (length \eqn{M}).
#'     \item \code{type.seq} - Vector of observed types of node pairs (as encoded through \link[ppsbm]{convertNodePair} of the events (length \eqn{M})).
#'     \item \code{Time} - [0,Time] is the total time interval of observation.
#'   }
#' @param tau Matrix with size \eqn{Q\times n} and values in \eqn{(0,1)}, containing the (estimated) probability that cluster \eqn{q} contains node \eqn{i}.
#' @param Q Total number of groups.
#' @param n Total number of nodes,  \eqn{1\le i \le n}.
#' @param directed Boolean for directed (TRUE) or undirected (FALSE) case
#' @param rho Either 1 (non sparse case) or vector with length \eqn{Q(Q+1)/2} (undirected case) or \eqn{Q^2} (directed case) with (estimated) values for the sparsity parameters \eqn{\rho^{(q,l)}}. See Section S6 in the supplementary material paper of Matias et al. (Biometrika, 2018) for more details.
#' @param sparse Boolean for sparse (TRUE) or not sparse (FALSE) case.
#' @param nb.points Number of points for the kernel estimation.
#'
#' @export
#'
#' @examples
#'
#' # The generated_sol_kernel solution was generated calling mainVEM
#' # with kernel method on the generated_Q3$data dataset.
#' # (50 individuals and 3 clusters)
#'
#' data <- generated_Q3$data
#'
#' n <- 50
#' Q <- 3
#'
#'
#' # Compute smooth intensity estimators
#' sol.kernel.intensities <- kernelIntensities(data,generated_sol_kernel$tau,Q,n,directed=FALSE)
#'
kernelIntensities <- function(data,tau,Q,n,directed,rho=1,sparse=FALSE,nb.points=1000*data$Time){
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
#' @param z Matrix of size  \eqn{Q \times n} with entries = 0 or 1: 'true' latent variables
#' @param hat.z Matrix of \eqn{Q \times n}  with 0<entries<1: estimated latent variables
#'
#' @references
#' HUBERT, L. & ARABIE, P. (1985). Comparing partitions. J. Classif. 2, 193–218.
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
#' Compute optimal matching between 2 clusterings using Hungarian algorithm.
#'
#'
#' @param z Matrix of size  \eqn{Q \times n}
#' @param hat.z Matrix of size  \eqn{Q \times n}
#'
#' @keywords internal
#'
#' @references
#'
#' HUBERT, L. & ARABIE, P. (1985). Comparing partitions. J. Classif. 2, 193–218.
#'
#' MATIAS, C., REBAFKA, T. & VILLERS, F. (2018).  A semiparametric extension of the stochastic block model for longitudinal networks. Biometrika. 105(3): 665-680.
#'
# # Example
# z <- matrix(c(1,1,0,0,0,0, 0,0,1,1,0,0, 0,0,0,0,1,1), nrow = 3, byrow = TRUE)
# hat.z <- matrix(c(0,0,1,1,0,0, 1,1,0,0,0,0, 0,0,0,0,1,1), nrow = 3, byrow = TRUE)
#
# perm <- permuteZEst(z,hat.z)
#
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
#' Sort intensities associated with the estimated clustering \eqn{\hat z} "in the same way" as the original intensities associated with true clustering \eqn{z} by permutation of rows.
#'
#' @param intensities Matrix whose rows contain piecewise constant intensities \eqn{\alpha^{(q,l)}}, given as a vector of values on a regular partition of the time interval.
#' @param z Matrix of size  \eqn{Q \times n}.
#' @param hat.z Matrix of size  \eqn{Q \times n}.
#' @param directed Boolean for directed (TRUE) or undirected (FALSE) case.
#'
#' @export
#'
#' @references
#'
#' HUBERT, L. & ARABIE, P. (1985). Comparing partitions. J. Classif. 2, 193–218.
#'
#' MATIAS, C., REBAFKA, T. & VILLERS, F. (2018).  A semiparametric extension of the stochastic block model for longitudinal networks. Biometrika. 105(3): 665-680.
#'
#' @examples
#' # True and estimated clusters for n=6 nodes clustered into Q=3 groups
#' z <- matrix(c(1,1,0,0,0,0, 0,0,1,1,0,0, 0,0,0,0,1,1), nrow = 3, byrow = TRUE)
#' hat.z <- matrix(c(0,0,1,1,0,0, 1,1,0,0,0,0, 0,0,0,0,1,1), nrow = 3, byrow = TRUE)
#'
#' # Set constant intensities for each directed pair (q,l)
#' intens <- matrix(c(1,1,1,2,2,2,3,3,3),9)
#'
#' # Permute the rows according to the permutation that "matches" hat.z with z
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
