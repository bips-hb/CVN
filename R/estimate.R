

#' @export 
estimate <- function(Theta0, Z0, Y0, D, m, p, Sigma, n_obs, rho, epsilon, maxiter, truncate, verbose = FALSE) { 
  
  # keep track whether the algorithm finished, either by 
  # whether the stopping condition was met, or the maximum
  # number of iterations was reached
  converged <- FALSE 
  iter <- 1 
  difference <- 1 # stores the discrepancy 
  
  # Initialize a Temp variable for Theta
  Temp <- Theta0
  Theta <- Theta0
  Z <- Z0
  Y <- Y0
  
  repeat{
    
    # Update Z -------------------------------------
    Z <- CVN::updateZ(m, p, Theta, Y, D) 
    
    # Update Y -------------------------------------
    Y <- CVN::updateY(Theta, Z, Y) 
    
    # Update Theta ---------------------------------
    Temp <- CVN::updateTheta(m, Z, Y, Sigma, n_obs, rho)
    Theta_previous <- Theta 
    Theta <- Temp
    
    # Check whether the algorithm is ready ----------
    difference <- CVN::relative_difference_precision_matrices(Theta, Theta_previous)
    
    if (verbose) { 
      if (((iter - 1) %% 10) == 0) { 
        cat("-------------------------\n") 
      }
      cat(sprintf("iteration %d  |  %f\n", iter, difference)) 
    }
    
    if (difference < epsilon) { 
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