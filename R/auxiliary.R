###################################################
## AUXILIARY FUNCTIONS
###################################################


###################################################
##   HANDLING OF NODES (i,j)
###################################################


#' Convert node pair \eqn{(i,j)}
#'
#' Convert node pair \eqn{(i,j)} into an index in \eqn{\{1,\dots,N\}} where \eqn{N=n(n-1)} (directed case) or \eqn{N=n(n-1)/2} (undirected case).
#' 
#' Interacting individuals \eqn{(i,j)} must be encoded into integer values in  \eqn{\{1,\dots,N\}} in describing the data. 
#'  
#' \describe{
#'   \item{\strong{Directed case :}}{
#'     \itemize{
#'       \item The node pair \eqn{(i,j)} with \eqn{(i\neq j)} is converted into the index \eqn{(i-1)*(n-1)+j-(i<j)}
#'     }
#'   }
#'   \item{\strong{Undirected case :}}{
#'     \itemize{
#'       \item The node pair \eqn{(i,j)} with \eqn{(i\neq j)} is converted into the index \eqn{(2*n-i)*(i-1)/2 +j-i}
#'     }
#'   }
#' }
#'
#' The number of possible node pairs is
#'     \itemize{
#'       \item \eqn{N = n*(n-1)} for the directed case
#'       \item \eqn{N = n*(n-1)/2}  for the undirected case
#'     }
#' which corresponds to the range of values for \code{data$type.seq}
#'
#' @param i Node \eqn{i} : \eqn{i\in \{1, \ldots, n\} }
#' @param j Node \eqn{j} : \eqn{j\in \{1, \ldots, n\} }
#' @param n Total number of nodes,  \eqn{1\le i \le n}
#' @param directed Boolean for directed (TRUE) or undirected (FALSE) case
#'
#' @return Index corresponding to the node pair
#'
#' @export
#'
#' @examples
#' # Convert the node pair (3,7) into an index, where the total number of nodes is 10,
#' # for directed and undirected interactions
#'
#' i <- 3
#' j <- 7
#' n <- 10
#'
#' directedIndex <- convertNodePair(i,j,n,TRUE)
#' undirectedIndex <- convertNodePair(i,j,n,FALSE)
#'
convertNodePair <- function(i,j,n,directed){
  if (sum((i>n) | (j>n))>0){
    stop("Your index is out of range")
  }
  if (directed){#directed case
    dyads = (i-1)*(n-1)+j-(i<j)
  } else {#undirected case
    dyads <- c(0,cumsum((n-1):1))[pmin(i,j)] + abs(j-i)
  }
  return(dyads)
}


#' List node pairs
#'
#' Create the list of all node pairs
#'
#' @param n Total number of nodes,  \eqn{1 \le i \le n}
#' @param directed Boolean for directed (TRUE) or undirected (FALSE) case
#'
#' @return Matrix with two columns which lists all the possible node pairs. Each row is a node pair.
#'
#' @export
#'
#' @examples
#' # List all the node pairs with 10 nodes, for directed and undirected graphs
#'
#' n <- 10
#' listNodePairs(n, TRUE)
#' listNodePairs(n, FALSE)
#'
listNodePairs <- function(n,directed=TRUE){
  N <- if (directed) n*(n-1) else n*(n-1)/2
  index <- matrix(0,N,2)
  if (directed){ # directed
    index[,1] <- rep(1:n,each=n-1)
    k <- (1:n^2)[-seq(1,n^2,by=n+1)]
    index[,2] <- rep(1:n,n)[k]
  }else { # undirected
    index[,1] <- rep(1:(n-1),times=(n-1):1)
    toto <- c()
    for (k in 1:(n-2)){
      toto <- c(toto,k*(n-1)+ 1:k)
    }
    index[,2] <- rep(2:n,n-1)[-toto]
  }
  return(index)
}



###################################################
##   HANDLING OF GROUP INDICES (q,l)
###################################################


#' Convert group pair \eqn{(q,l)}
#'
#' Gives the index in \eqn{1, \ldots, Q^2} (directed) or \eqn{1, \ldots, Q(Q+1)/2} (undirected) that corresponds to group pair \eqn{(q,l)}. Works also for vectors of indices \eqn{q} and \eqn{l}.
#'
#' Relations between groups \eqn{(q,l)} are stored in vectors, whose indexes depend on whether the graph is directed or undirected.
#' \describe{
#'   \item{\strong{Directed case :}}{
#'     \itemize{
#'       \item The \eqn{(q,l)} group pair is converted into the index \eqn{(q-1)Q+l}
#'     }
#'   }
#'   \item{\strong{Undirected case :}}{
#'     \itemize{
#'       \item The \eqn{(q,l)} group pair with \eqn{q\leq l} is converted into the index \eqn{(2Q-q+2)*(q-1)/2 +l-q+1}
#'     }
#'   }
#' }
#' @param q Group index \eqn{q}
#' @param l Group index \eqn{l}
#' @param Q Total number of groups \eqn{Q}
#' @param directed Boolean for directed (TRUE) or undirected (FALSE) case
#'
#' @return Index corresponding to the group pair \eqn{(q,l)}
#'
#' @export
#'
#' @examples
#' # Convert the group pair (3,2) into an index, where the total number of groups is 3,
#' # for directed and undirected graph
#'
#' q <- 3
#' l <- 2
#' Q <- 3
#'
#' directedIndex <- convertGroupPair(q,l,Q)
#' undirectedIndex <- convertGroupPair(q,l,Q, FALSE)
#'
convertGroupPair <- function(q,l,Q,directed=TRUE){
  if (sum((q>Q) | (l>Q))>0){
    stop("Your index is out of range")
  }
  if (directed){
    index <- (q-1)*Q+l
  } else { # undirected
    qp <- pmin(q,l)
    lp <- pmax(q,l)
    index <- (2*Q-qp+2)*(qp-1)/2 +lp-qp+1
  }
  return(index)
}



#' Convert index into group pair
#'
#' This function is the inverse of the conversion \eqn{\{(q,l), 1 \le q,l\le Q \} } into \eqn{\{1,...,Q^2\}} for the directed case and of \eqn{\{(q,l), 1 \le q \le  l \le Q \}} into \eqn{\{1,...,Q(Q+1)/2\}} for the undirected case.
#' It takes the integer index corresponding to \eqn{(q,l)} and returns the pair \eqn{(q,l)}.
#'
#' @param ind_ql Converted \eqn{(q,l)} index
#' @param Q Total number of groups \eqn{Q}
#' @param directed Boolean for directed (TRUE) or undirected (FALSE) case
#'
#' @return Group pair \eqn{(q,l)} corresponding to the given index
#'
#' @export
#'
#' @examples
#' # Convert the index 5 into a group pair for undirected graph
#' # and the index 8 into a group pair for directed graph
#' # where the total number of groups is 3
#'
#' ind_ql_dir <- 8
#' ind_ql_undir <- 5
#'
#' Q <- 3
#'
#' directedIndex <- find_ql(ind_ql_dir,Q)
#' undirectedIndex <- find_ql(ind_ql_undir,Q, FALSE)
#'
find_ql <- function(ind_ql, Q,directed=TRUE){
  # test coherence
  nb.groups.ql <- if (directed)  Q^2 else Q*(Q+1)/2
  if (ind_ql > nb.groups.ql) stop("Your index is out of range")

  if (directed){ # directed
    q <- ceiling(ind_ql/Q)
    l <- ind_ql - Q*(q-1)
  }else{ # undirected
    w <- cumsum(Q:1)
    q <- which.max(ind_ql<=w)
    w <- c(0,w)
    l <- ind_ql - w[q] + q - 1
  }
  return(c(q,l))
}



#' Convert index into group pair in tauDown_Q
#'
#' This function is the inverse of the conversion \eqn{{(q,l), q<l}} into \eqn{{1,...,Q*(Q-1)/2}}. Used only in tauDown_Q.
#'
#' @param ind_ql Converted \eqn{(q,l)} index
#' @param Q Total number of groups \eqn{Q}
#'
#' @return Group pair \eqn{(q,l)} corresponding to the given index
#'
#' @keywords internal
#'
find_ql_diff <- function(ind_ql,Q){
  if (ind_ql > Q*(Q-1)/2){
    stop("Your index is out of range")
  }
  w <- cumsum((Q-1):1)
  q <- which.max(ind_ql<=w)
  w <- c(0,w)
  l <- ind_ql - w[q] + q
  return(c(q,l))
}

###################################################
##   HANDLING OF VALUES OF TAU
###################################################


#' Handling of values of \eqn{\tau}
#'
#' Avoid values of \eqn{\tau} to be exactly 0 and exactly 1.
#'
#' @param tau \eqn{\tau}
#'
#' @keywords internal
#'
correctTau <- function(tau){
  tau <- pmin(tau,.Machine$double.xmax)
  tau <- pmax(tau,.Machine$double.xmin)
  tau <- tau/sum(tau)
  tau <- pmin(tau,1-1e-7)
  tau <- pmax(tau,1e-7)
  tau <- tau/sum(tau)

  return(tau)
}


###################################################
##   COMPUTE STATISTICS FROM DATA
###################################################


#' Compute statistics
#'
#' Convert the initial data into the statistics matrix \eqn{N_{(ij),k}}, by counting the number of events for the nodes pair types \eqn{(i,j)} during the \eqn{k}-th subinterval of a regular partition (in \eqn{K} parts) of the time interval.
#'
#'
#' @param data List with 3 components \code{$type.seq, $time.seq, $Time} (see \link[ppsbm]{mainVEM} for more details). 
#' @param n Total number of nodes,  \eqn{1\le i \le n}.
#' @param K Size of the regular partition, i.e. number of subintervals of the time interval. When used as input in the VEM algorithm (with \code{hist} method), \eqn{K} must be a power of 2.
#' @param directed Boolean for directed (TRUE) or undirected (FALSE) case.
#'
#' @return N[(i,j),k] = matrix with \eqn{K} columns, each row contains the number of events for the node pair \eqn{(i,j)} during the k-th subinterval
#'
#' @export
#'
#' @examples
#' # Convert the generated data into the statistics matrix N_ijk with 8 subintervals
#'
#' n <- 50
#' K <- 2^3
#'
#' obs <- statistics(generated_Q3$data,n,K,directed=FALSE)
#'
statistics <- function(data,n,K,directed=TRUE){
  partition <- data$Time*seq(1/K,1,by=1/K)
  N <- if (directed)  n*(n-1) else n*(n-1)/2
  Nijk <- matrix(0,N,K)
  for (ind in 1:N){
    events.ij <- data$time.seq[data$type.seq==ind]
    nb.ij <- length(events.ij)
    if (nb.ij>0){
      counts <- rowSums(matrix(events.ij,nrow=K,ncol=nb.ij,byrow=TRUE)<matrix(partition,nrow=K,ncol=nb.ij))
      Nijk[ind,] <- if (K>1) counts-c(0,counts[1:(K-1)]) else counts
    }
  }
  return(Nijk)
}

