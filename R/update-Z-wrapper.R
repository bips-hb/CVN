#' Wrapper for the \eqn{Z}-update Step for the ADMM
#' 
#' A wrapper for the \code{C} function that returns the 
#' updated value of \eqn{Z} for the ADMM given the previously 
#' updated values of \eqn{\Theta} and \eqn{Y}
#' 
#' @param m Number of graphs
#' @param p Number of variables
#' @param nrow_D Number of rows of the \eqn{D}-matrix
#' @param Theta A list with matrices with the current values of \eqn{\Theta}
#' @param Y A list with matrices with the current values of \eqn{Y} 
#' @param D A matrix used for solving the Generalized LASSO 
#' 
#' @return A list with matrices with the new values of \eqn{Z}
#'
#' @seealso \code{\link{create_matrix_D}} 
#' 
#' @export
updateZ_wrapper <- function(m, p, nrow_D, 
                            Theta, Y, W, eta1, eta2, a, 
                            rho_genlasso, maxiter_genlasso, eps_genlasso, 
                            truncate_genlasso) { 
  
  updateZRcpp(m, p, nrow_D, 
               Theta, Y, W, eta1, eta2, a, 
               rho_genlasso, maxiter_genlasso, eps_genlasso, 
               truncate_genlasso)  
  
}
  