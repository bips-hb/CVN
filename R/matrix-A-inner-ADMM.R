#' Determine matrix \eqn{A} for inner-ADMM for the \eqn{Z}-update step
#' 
#' The \eqn{Z}-update step requires 
#' us to solve a special Generalized LASSO problem of the form 
#' \deqn{
#'   \hat{\beta} = \text{argmin } \frac{1}{2} || y - \beta ||_2^2 + ||D\beta||_1 
#' }
#' where \eqn{\beta} and \eqn{y} are \eqn{m}-dimensional vectors and 
#' \eqn{D} is a \eqn{(c \times m)}-matrix where \eqn{c = (m^2 + m) / 2}. 
#' We solve this optimization problem using an adaption of the ADMM
#' algorithm presented in Zhu (2017). 
#' This algorithm requires the choice of a matrix \eqn{A} such that 
#' \eqn{A - D'D} is positive semidefinite. In order to optimize the ADMM, 
#' we choose the matrix \eqn{A} to be diagonal with a fixed value \eqn{a}. 
#' This function determines the smallest value of \eqn{a} such that 
#' \eqn{A - D'D} is indeed positive semidefinite. We do this be determining 
#' the largest eigenvalue
#' 
#' @param W Weight matrix \eqn{W}
#' @param eta1,eta2 The values \eqn{\eta_1 = \lambda_1 / \rho} and 
#'                   \eqn{\eta_2 = \lambda_2 / \rho}
#' 
#' @return Value of \eqn{a}
#' @references 
#' Zhu, Y. (2017). An Augmented ADMM Algorithm With Application to the 
#' Generalized Lasso Problem. Journal of Computational and Graphical Statistics, 
#' 26(1), 195–204.\cr \doi{10.1080/10618600.2015.1114491}
#' @author Louis Dijkstra
#' 
#' @export
matrix_A_inner_ADMM <- function(W, eta1, eta2) { 
  
  # determine number of graphs
  m <- ncol(W) 
  
  # determine D'D 
  DtD <- -eta2^2 * (W*W) + eta1^2 * diag(m) + eta2^2*diag(rowSums(W*W))
  
  # determine the eigenvalues
  eigenvectors <- eigen(DtD)
  
  # return the value of a
  return(max(eigenvectors$values))
}