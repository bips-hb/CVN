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
updateZ <- function(m, p, Theta, Y, D, n_cores = 1) { 
  
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
  combinations <- combn(1:p, 2, simplify = FALSE)
  
  n_combinations <- length(combinations)
  
  # go over all combinations
  for (k in 1:n_combinations) { 
    
    # get the individual indices (s,t)
    s <- combinations[[k]][1]
    t <- combinations[[k]][2]
    
    # get the y-vector for the generalized LASSO
    y <- c()
    for (i in 1:m) { 
       y[i] <- B[[i]][s, t]
    }
    
    # apply the generalized LASSO 
    out <- genlasso(y, diag(1, m), D, minlam = 1)
    beta <- coef(out, lambda = 1)$beta

    #beta[abs(beta) <= 1e-10] <- 0
    # 
    #print(beta)
    # 
    # est = optim(par = rep(0,r+m), fn, method = "L-BFGS-B", lower = rep(-1,r+m), upper = rep(1,r+m), y = y, D = D)
    # b = y - t(D) %*% est$par
    # #
    # b[abs(b) <= 1e-10] <- 0
    # #print("value b")
    # #print(b)
    # beta <- b

    # update the matrix Z (use that it is symmetric)
    for (i in 1:m) { 
      Z[[i]][s, t] <- beta[i] 
      Z[[i]][t, s] <- beta[i] 
    }
  }
  
  return(Z)
}