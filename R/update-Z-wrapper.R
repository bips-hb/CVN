#' Wrapper for the \eqn{Z}-update Step for the ADMM
#'
#' A wrapper for the \code{C} function that returns the
#' updated value of \eqn{Z} for the ADMM given the previously
#' updated values of \eqn{\Theta} and \eqn{Y}
#'
#' @param m Number of graphs
#' @param p Number of variables
#' @param Y A list with matrices with the current values of \eqn{Y} 
#' @param W Weight matrix 
#' @param eta1  lambda1 / rho, with \eqn{\rho} is penalty parameter for the global ADMM algorithm (Default: \code{1})
#' @param eta2 
#' @param rho_genlasso The \eqn{\rho} penalty parameter for the ADMM algorithm 
#' @param eps_genlasso If the relative difference between two update steps is 
#'                smaller than \eqn{\epsilon}, the algorithm stops
#' @param maxiter_genlasso Maximum number of iterations for solving
#'                the generalized LASSO problem
#' @param truncate_genlasso All values of the final \eqn{\hat{\beta}} below
#'                 \code{truncate_genlasso} will be set to \code{0}.
#'
#' @return A list with matrices with the new values of \eqn{Z}
#' @noRd
#' 
#' @keywords internal
updateZ_wrapper <- function(m, p,  
                            Theta, Y, W, eta1, eta2, a, 
                            rho_genlasso, maxiter_genlasso, eps_genlasso, 
                            truncate_genlasso=1e-4) { 
  
  updateZRcpp(m, p, 
               Theta, Y, W, eta1, eta2, a, 
               rho_genlasso, maxiter_genlasso, eps_genlasso, 
               truncate_genlasso)  
  
}
