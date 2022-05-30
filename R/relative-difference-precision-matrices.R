#' The Relative Difference between two Precision Matrices
#' 
#' Returns the relative \eqn{L1} difference between precision matrix \eqn{\Theta(k+1)}
#' (parameter \code{Theta_new}) and \eqn{\Theta(k)} (parameter \code{Theta_old}).
#' 
#' This is used for checking whether the stopping condition has been met.
#' 
#' @param Theta_new A list with matrices with the updated values of \eqn{\Theta}
#' @param Theta_old A list with matrices with the old values of \eqn{\Theta}
#' @param n_cores Number of cores used (Default: 1)
#' 
#' @return The relative difference between \eqn{\Theta(k+1)} and \eqn{\Theta(k)}
#'        
#' @export
relative_difference_precision_matrices <- function(Theta_new, Theta_old, n_cores = 1) { 
  sum( mcmapply(function(theta_new, theta_old) {norm(theta_new - theta_old)}, 
             Theta_new, Theta_old, mc.cores = n_cores)) / sum( mcmapply(norm, Theta_old, mc.cores = n_cores) )
}