#' Generated graph with 50 individuals and 3 clusters
#'
#' @format A data frame
#' \describe{
#'   \item{data}{List of 3}
#'   \item{z}{Latent variables}
#'   \item{intens}{Intensities}
#' }
"generated_Q3"

#' Generated graph with 20 individuals and 3 clusters
#'
#' @format A data frame
#' \describe{
#'   \item{data}{List of 3}
#'   \item{z}{Latent variables}
#'   \item{intens}{Intensities}
#' }
"generated_Q3_n20"

#' Generated solution with histogram method
#'
#' @format List of 5 iterations of the algorithm, each one containing
#' \describe{
#'   \item{List of 8}{tau, rho, beta, logintensities.ql, best.d, J, run, converged}
#' }
"generated_sol_hist"

#' Generated solution with kernel method
#'
#' @format Solution containing
#' \describe{
#'   \item{List of 8}{tau, logintensities.ql.ij, J, run, converged}
#' }
"generated_sol_kernel"

