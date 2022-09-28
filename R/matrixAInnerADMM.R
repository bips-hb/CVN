#' Determine matrix \eqn{A} for inner-ADMM for the \eqn{Z}-update step
#' 
#' The \eqn{Z}-update step, see \code{\link{updateZ}}, requires 
#' us to solve a special Generalized LASSO problem of the form 
#' \deqn{
#'   \hat{\beta} = \text{argmin } \frac{1}{2} || y - \beta ||_2^2 + ||D\beta||_1 
#' }
#' where \eqn{\beta} and \eqn{y} are \eqn{m}-dimensional vectors and 
#' \eqn{D} is a \eqn{(c \times m)}-matrix where \eqn{c \geq m}. 
#' We solve this optimization problem using an adaption of the ADMM
#' algorithm presented in Zhu (2017). 
#' This algorithm requires the choice of a matrix \eqn{A} such that 
#' \eqn{A - D'D} is positive semidefinite. In order to optimize the ADMM, 
#' we choose the matrix \eqn{A} to be diagonal with a fixed value \eqn{a}. 
#' This function determines the smallest value of \eqn{a} such that 
#' \eqn{A - D'D} is indeed positive semidefinite. It does so by solving 
#' a semidefinite program very much related to the student educational problem. 
#' 
#' @param m Number of graphs
#' @param D Matrix D, see \code{\link{create_matrix_D}}
#' 
#' @return Value of \eqn{a}
#' @references 
#' Zhu, Y. (2017). An Augmented ADMM Algorithm With Application to the 
#' Generalized Lasso Problem. Journal of Computational and Graphical Statistics, 
#' 26(1), 195–204. https://doi.org/10.1080/10618600.2015.1114491
#' 
#' https://math.stackexchange.com/questions/665026/adding-elements-to-diagonal-of-symmetric-matrix-to-ensure-positive-definiteness
#' @export
matrix_A_inner_ADMM <- function(m, D) { 
  
  # create the variables
  A <- CVXR::Variable(m, m) 
  a <- CVXR::Variable(1) 
  
  # objective function
  objective <- CVXR::Minimize(abs(sum(A)) / m)
  
  # constraints:
  R <- t(D) %*% D
  
  #Q <- diag(a,m)
  #print(Q)
  
  #constraint1 <- A == Q # A must be a diagonal matrix with fixed a
  #constraint2 <- A %>>% R       # A - D'D must be positive semidefinite 
  
  #constraints <- lapply(1:m, function(i) A[i,i] == a)
  #constraints <- append(constraints, constraint2)
  
  # define the problem using CVXR:
  
  problem <- CVXR::Problem(objective, constraints = list(A == diag(a,m), A %>>% R))
  #problem <- CVXR::Problem(objective, constraints = list(constraint1, constraint2))
  #problem <- CVXR::Problem(objective, constraints = constraints)
  
  
  # solve
  solution <- CVXR::solve(problem, solver = "SCS")
  
  # return the value of a
  return(solution$getValue(a))
}