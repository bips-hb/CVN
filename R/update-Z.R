#' The \eqn{Z}-update step for the ADMM
#' 
#' Returns the updated value of \eqn{Z} for the ADMM given the previously 
#' updated values of \eqn{\Theta} and \eqn{Y}
#' 
#' @param m Number of graphs
#' @param Theta A list with matrices with the current values of \eqn{\Theta}
#' @param Y A list with matrices with the current values of \eqn{Y} 
#' @param D A matrix used for solving the Generalized LASSO 
#' @param n_cores Number of cores used (Default: 1)
#' 
#' @return A list with matrices with the new values of \eqn{Z}
#'
#' @seealso \code{\link{create_matrix_D}}        
#' @export
updateZ <- function(m, Theta, Y, D, n_cores = 1) { 
  
  # Initialize the Z matrices and immediately update the diagonal 
  # entries given Theta and Y
  Z <- mapply(function(theta, y) { diag( diag(theta) + diag(y) ) }, 
              Theta, Y, SIMPLIFY = FALSE)
  
  # Here, we use the fact that we can update the matrices Z for each 
  # individual entry (s,t)  by solving a generalized LASSO problem. 
  # In addition, we use the fact that these matrices are symmetric and 
  # we, therefore, need not solve the entire matrix, but only the upper 
  # diagonal. 
  
  B <- mapply('+', Theta, Y, SIMPLIFY = FALSE)
  
  # go over all unique pairs 
  combinations <- combn(1:m, 2, simplify = FALSE)
  
  mclapply(combinations, function(combination) {
      # obtain the vector y for the generalized LASSO
      y <- sapply(B, function(M) M[combination[1], combination[2]])
      
      # apply the generalized LASSO 
      out <- genlasso(y, diag(m, 1), D, minlam = 1)
      beta <- coef(out, lambda = 1)$beta
      
      # update the matrix Z (use that it is symmetric)
      for (i in 1:m) { 
        Z[[i]][s,t] <<- values[i] 
        Z[[i]][t,s] <<- values[i] 
      }
    }, mc.cores = n_cores)
  
  return(Z)
}