#' The \eqn{\Theta}-update step for the ADMM
#' 
#' \code{updateTheta} returns the updated value of \eqn{\Theta} for the 
#' ADMM given the previously updated values of \eqn{Z} and \eqn{Y}
#' 
#' @param m Number of graphs
#' @param Z A list with matrices with the current values of \eqn{Z}
#' @param Y A list with matrices with the current values of \eqn{Y}
#' @param Sigma A list with empirical covariance matrices \eqn{{\Sigma}}
#' @param n_obs A \eqn{m}-dimensional vector with the number of observations per graph 
#' @param rho The \eqn{\rho} ADMM's penalty parameter (Default: 1) 
#' 
#' @return A list with matrices with the new values of \eqn{\Theta}
#' @noRd
#'        
#' @keywords internal
updateTheta <- function(m, Z, Y, Sigma, n_obs, rho = 1) { 
  
  # we use the fact that we can separate the optimization problem for each 
  # individual graph i = 1, 2,..., m 
  lapply(1:m, function(i) { 
    # perform an eigendecomposition 
    eigen_decomposition <- eigen(Sigma[[i]] - rho * Z[[i]]/ n_obs[i] + rho * Y[[i]] / n_obs[i])
    
    # obtain matrices Q and Lambda (Lambda is a diagonal matrix, so first only the diagonal is stored)
    Q <- eigen_decomposition$vectors 
    Lambda <- eigen_decomposition$values
    
    # Lambda is updated 
    Lambda <- n_obs[i]/(2 * rho) * (-Lambda + sqrt(Lambda^2 + 4*rho/n_obs[i]))
    
    # The new Theta matrix
    Q %*% diag(Lambda) %*% t(Q)
  })
}