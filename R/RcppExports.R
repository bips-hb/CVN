# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Generalized LASSO
#' 
#' Solving Generalized LASSO with fixed \eqn{\lambda = 1}
#' Solves efficiently the generalized LASSO problem of the form
#' \deqn{
#'   \hat{\beta} = \text{argmin } \frac{1}{2} || y - \beta ||_2^2 + ||D\beta||_1
#' }
#' where \eqn{\beta} and \eqn{y} are \eqn{m}-dimensional vectors and
#' \eqn{D} is a \eqn{(c \times m)}-matrix where \eqn{c \geq m}.
#' We solve this optimization problem using an adaption of the ADMM
#' algorithm presented in Zhu (2017).
#'
#' @param Y The \eqn{y} vector of length \eqn{m}
#' @param W The weight matrix \eqn{W} of dimensions \eqn{m x m}
#' @param m The number of graphs
#' @param eta1 Equals \eqn{\lambda_1 / rho}
#' @param eta2 Equals \eqn{\lambda_2 / rho}
#' @param a Value added to the diagonal of \eqn{-D'D} so that
#'          the matrix is positive definite, see
#'          \code{\link{matrix_A_inner_ADMM}}
#' @param rho The ADMM's parameter
#' @param max_iter Maximum number of iterations
#' @param eps Stopping criterion. If differences
#'            are smaller than \eqn{\epsilon}, algorithm
#'            is halted
#' @param truncate Values below \code{truncate} are
#'                 set to \code{0}
#' @return The estimated vector \eqn{\hat{\beta}}
#' @author Louis Dijkstra
#' @references
#' Zhu, Y. (2017). An Augmented ADMM Algorithm With Application to the
#' Generalized Lasso Problem. Journal of Computational and Graphical Statistics,
#' 26(1), 195–204. https://doi.org/10.1080/10618600.2015.1114491
#'
#' @seealso \code{\link{genlasso_wrapper}}
genlassoRcpp <- function(Y, W, m, eta1, eta2, a, rho, max_iter, eps, truncate) {
    .Call('_CVN_genlassoRcpp', PACKAGE = 'CVN', Y, W, m, eta1, eta2, a, rho, max_iter, eps, truncate)
}

#' The \eqn{Z}-update Step
#' 
#' A \code{C} implementation of the \eqn{Z}-update step. We
#' solve a generalized LASSO problem repeatedly for each of the
#' individual edges
#'
#' @param m The number of graphs
#' @param p The number of variables
#' @param Theta A list of matrices with the \eqn{\Theta}-matrices
#' @param Y A list of matrices with the \eqn{Y}-matrices
#' @param W The weight matrix \eqn{W} of dimensions \eqn{m x m}
#' @param eta1 Equals \eqn{\lambda_1 / rho}
#' @param eta2 Equals \eqn{\lambda_2 / rho}
#' @param a Value added to the diagonal of \eqn{-D'D} so that
#'          the matrix is positive definite, see
#'          \code{\link{matrix_A_inner_ADMM}}
#' @param rho The ADMM's parameter
#' @param max_iter Maximum number of iterations
#' @param eps Stopping criterion. If differences
#'            are smaller than \eqn{\epsilon}, algorithm
#'            is halted
#' @param truncate Values below \code{truncate} are
#'                 set to \code{0}
#' @author Louis Dijkstra
#' @return The estimated vector \eqn{\hat{\beta}}
#'
updateZRcpp <- function(m, p, Theta, Y, W, eta1, eta2, a, rho, max_iter, eps, truncate) {
    .Call('_CVN_updateZRcpp', PACKAGE = 'CVN', m, p, Theta, Y, W, eta1, eta2, a, rho, max_iter, eps, truncate)
}

