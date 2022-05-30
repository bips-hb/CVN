#' The \eqn{Z}-update step for the ADMM
#' 
#' Returns the updated value of \eqn{Z} for the ADMM given the previously 
#' updated values of \eqn{\Theta} and \eqn{Y}
#' 
#' @param m Number of graphs
#' @param Theta A list with matrices with the current values of \eqn{\Theta}
#' @param Z A list with matrices with the current values of \eqn{Z}
#' @param Y A list with matrices with the current values of \eqn{Y} 
#' @param D A matrix used for solving the Generalized LASSO 
#' @param n_obs A \eqn{m}-dimensional vector with the number of observations per graph 
#' @param rho The \eqn{\rho} ADMM's penalty parameter (Default: 1) 
#' @param n_cores Number of cores used (Default: 1)
#' 
#' @return A list with matrices with the new values of \eqn{Z}
#'
#' @seealso \code{\link{create_matrix_D}}        
#' @export
updateZ <- function(m, Theta, Z, Y, D, n_obs, rho = 1, n_cores = 1) { 
  # Y_new = Y_old + Theta_new - Z_new
  #mapply(function(y, theta, z) {y + theta - z}, 
  #       Y_old, Theta_new, Z_new, SIMPLIFY = FALSE)
}