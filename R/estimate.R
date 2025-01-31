#' Estimate a CVN
#' 
#' A function for estimating a CVN for a single \eqn{(\lambda_1, \lambda_2)}-value. 
#' See for more details \code{\link{CVN}}
#' @keywords internal
estimate <- function(m, p, W, Theta0, Z0, Y0, a, eta1, eta2, Sigma, n_obs, 
                     rho, rho_genlasso, eps, eps_genlasso,
                     maxiter, maxiter_genlasso, truncate, 
                     truncate_genlasso,
                     verbose = FALSE) {
  
  # keep track whether the algorithm finished, either by 
  # whether the stopping condition was met, or the maximum
  # number of iterations was reached
  converged  <- FALSE 
  iter       <- 1 
  difference <- 1 # stores the discrepancy 
  
  # Initialize a Temp variable for Theta
  Temp  <- Theta0
  Theta <- Theta0
  Z     <- Z0
  Y     <- Y0
  
  repeat{
    
    # Update Z -------------------------------------
    Z <- updateZRcpp(m, p, 
                     as.matrix(Theta), as.matrix(Y), W, eta1, eta2, a, 
                     rho_genlasso, maxiter_genlasso, eps_genlasso, 
                     truncate_genlasso)  
    Z <- as.list(Z[,1])
    
    # Update Y -------------------------------------
    Y <- updateY(Theta, Z, Y) 
    
    # Update Theta ---------------------------------
    Temp <- updateTheta(m, Z, Y, Sigma, n_obs, rho)
    Theta_previous <- Theta 
    Theta <- Temp
    
    # Check whether the algorithm is ready ----------
    difference <- relative_difference_precision_matrices(Theta, Theta_previous)
    
    if (verbose) { 
      if (((iter - 1) %% 10) == 0) { 
        cat("-------------------------\n") 
      }
      cat(sprintf("iteration %d  |  %f\n", iter, difference)) 
    }
    
    if (difference < eps) { 
      converged <- TRUE
      iter <- iter + 1
      break() 
    }
    
    if (iter >= maxiter) { 
      warning("Maximum number of iterations reached. Stopping criterion not met")
      break()
    }
    
    iter <- iter + 1
  }
  
  Z <- lapply(Z, function(z) { 
    ifelse(abs(z) <= truncate, 0, z)
  })
  
  adj_matrices <- lapply(Z, function(X) { 
    diag(X) <- 0 
    Matrix( as.numeric( abs(X) >= 2*.Machine$double.eps), ncol = ncol(X) , sparse = TRUE)
  })
  
  list(
    Theta = Theta,
    Z = Z, 
    Y = Y, 
    adj_matrices = adj_matrices, 
    converged = converged, 
    value = difference, 
    n_iterations = iter
  )
}