#' Wrapper for \code{genlassoRcpp}
#'
#' See for details \code{\link{genlassoRcpp}}
#' 
#' @param y A numeric vector 
#' @param W Weight matrix 
#' @param m Number of graphs
#' @param c constant
#' @param eta1 lambda1 / rho, with \eqn{\rho} is penalty parameter for the global ADMM algorithm (Default: \code{1})
#' @param eta2 lambda2 / rho
#' @param a  constant
#' @param rho the \eqn{\rho} penalty parameter for the ADMM algorithm 
#' @param max_iter Maximum number of iterations (Default: \code{100})
#' @param eps If the relative difference between two update steps is 
#'                smaller than \eqn{\epsilon}, the algorithm stops. 
#'                See \code{\link[CVN]{relative_difference_precision_matrices}}
#'                (Default: \code{1e-4})
#' @param truncate All values of the final \eqn{\hat{\Theta}_i}'s below \code{truncate} will be 
#'                 set to \code{0}. (Default: \code{1e-4})
#'
#' @seealso \code{\link[CVN]{relative_difference_precision_matrices}}
#' @keywords internal
#'
#' @seealso \code{\link{genlassoRcpp}}
#' @export
# FIXME: This passes an addiitonal argument 'c' that is not listed in the definition
# of genlassoRcpp in CVN.cpp
genlasso_wrapper <- function(y, W, m, c, eta1, eta2, a, rho, max_iter, eps, truncate) {
  genlassoRcpp(y, W, m, c, eta1, eta2, a, rho, max_iter, eps, truncate)
}
