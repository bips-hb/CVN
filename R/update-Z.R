#' The \eqn{Z}-update step for the ADMM
#' 
#' Returns the updated value of \eqn{Z} for the ADMM given the previously 
#' updated values of \eqn{\Theta} and \eqn{Y}
#' 
#' @param m Number of graphs
#' @param p Number of variables
#' @param Theta A list with matrices with the current values of \eqn{\Theta}
#' @param Y A list with matrices with the current values of \eqn{Y} 
#' @param D A matrix used for solving the Generalized LASSO 
#' @param n_cores Number of cores used (Default: 1)
#' 
#' @return A list with matrices with the new values of \eqn{Z}
#'
#' @seealso \code{\link{create_matrix_D}}        
#' @export
updateZ <- function(m, p, nrow_D, 
                    Theta, Y, W, eta1, eta2, a, 
                    rho_genlasso, maxiter_genlasso, eps_genlasso, 
                    use_genlasso_package) { 
  
  # Initialize the Z matrices and immediately update the diagonal 
  # entries given Theta and Y
  Z <- mapply(function(theta, y) { diag( diag(theta) + diag(y) ) }, 
              Theta, Y, SIMPLIFY = FALSE)
  
  # Here, we use the fact that we can update the matrices Z for each 
  # individual entry (s,t)  by solving a generalized LASSO problem. 
  # In addition, we use the fact that these matrices are symmetric and 
  # we, therefore, need not solve the entire matrix, but only the upper 
  # diagonal
  B <- lapply(1:m, function(i) { Theta[[i]] + Y[[i]] })
  
  # go over all unique pairs 
  combinations <- combn(1:p, 2, simplify = FALSE)
  
  n_combinations <- length(combinations)
  
  # go over all combinations
  for (k in 1:n_combinations) { 
  #parLapply(cl, 1:n_combinations, function(k) { 
    
    # get the individual indices (s,t)
    s <- combinations[[k]][1]
    t <- combinations[[k]][2]
    
    # get the y-vector for the generalized LASSO
    y <- c()
    for (i in 1:m) { 
       y[i] <- B[[i]][s, t]
    }
    
    # apply the generalized LASSO 
    if (use_genlasso_package) { 
      D <- CVN::create_matrix_D(W, eta1, eta2) 
      out <- genlasso::genlasso(y, diag(1, m), D, minlam = 1)
      beta <- coef(out, lambda = 1)$beta
    } else {
      beta <- aug_genlassoRcpp(y, W, m, nrow_D, eta1, eta2, a,
                               rho_genlasso, maxiter_genlasso,
                               eps_genlasso)
    }
    
    # update the matrix Z (use that it is symmetric)
    for (i in 1:m) { 
      Z[[i]][s, t] <- beta[i] 
      Z[[i]][t, s] <- beta[i] 
    }
  }
  
  return(Z)
}