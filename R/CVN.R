#' Create a CVN Data Object
#' 
#' Creates a Covariate-varying network (CVN) object for estimating the 
#' graphs. It only requires the \code{raw_data}: a list of \code{m} different
#' matrices, each containing the observed data. Each column is a variable. Therefore,
#' each matrix should have the same number of columns. The number of observations
#' can differ between the datasets. 
#' 
#' @param raw_data A list with matrices. The number of columns should be the 
#'                 same for each matrix
#' @param W The \eqn{(m x m)}-dimensional uppertriangular weight matrix
#' @param normalized Data is normalized if TRUE. Otherwise the data is only
#'                   centered (Default: TRUE)
#' @param Verbose Verbose (Default: FALSE) 
#' 
#' @return A CVN data object; a list with entries
#'   \item{\code{data}}{The \code{raw_data}, but then normalized or centered}
#'   \item{\code{W}}{The \eqn{m x m} weight matrix}
#'   \item{\code{m}}{Number of graphs}
#'   \item{\code{p}}{Number of variables}
#'   \item{\code{n_obs}}{Number of observations for each graph}
#'        
#' @export
# TODO: change lambda1 and lambda2 to a grid!
CVN <- function(data, W, lambda1 = 1, lambda2 = 1, 
                rho = 1,
                epsilon = 10^(-5),
                maxiter = 1000, 
                n_cores = 1, 
                normalized = FALSE, 
                verbose = FALSE) { 
  
  # Check correctness input -------------------------------
  check_correctness_input(data, W, lambda1, lambda2, rho)
  
  # Extract variables -------------------------------------
  m <- length(data)       # total number of graphs  
  p <- ncol(data[[1]])    # total number of variables
  n_obs <- sapply(data, function(X) nrow(X))    # no. of observations per graph 
  
  # Center or normalize data -------------------------------
  data <- lapply(data, function(X) scale(X, scale = normalized))
  
  # Compute the empirical covariance matrices --------------
  Sigma <- lapply(data, cov) 
  
  # Initialize variables for the algorithm -----------------
  # Generate matrix D for the generalized LASSO 
  D <- CVN::create_matrix_D(W, lambda1, lambda2, rho)
  
  # Generate matrices for ADMM iterations
  Theta_old <- rep(list(diag(p)), m) # m (p x p)-dimensional identity matrices
  Theta_new <- rep(list(diag(p)), m)
  
  #Z_old <- rep(list(matrix(0, nrow = p, ncol = p)), m) # m (p x p)-dimensional zero matrices
  #Z_new <- rep(list(matrix(0, nrow = p, ncol = p)), m)
  Z <- rep(list(matrix(0, nrow = p, ncol = p)), m) # m (p x p)-dimensional zero matrices
  
  Y <- rep(list(matrix(0, nrow = p, ncol = p)), m) # m (p x p)-dimensional zero matrices
  
  #Y_old <- rep(list(matrix(0, nrow = p, ncol = p)), m) # m (p x p)-dimensional zero matrices
  #Y_new <- rep(list(matrix(0, nrow = p, ncol = p)), m)
  
  # Initialize a Temp variable for these matrices 
  #Temp <- rep(list(matrix(0, nrow = p, ncol = p)), m)
  
  # keep track whether the algorithm finished, either by 
  # whether the stopping condition was met, or the maximum
  # number of iterations was reached
  converged <- FALSE 
  iter <- 1 
  
  repeat{
    
    # Update Theta ---------------------------------
    # TDOD!!!!!! Make sure that Z_old and Y_old are updated!
    Temp <- updateTheta(m, Z, Y, Sigma, n_obs, rho, n_cores = n_cores)
    Theta_old <- Theta_new 
    Theta_new <- Temp 
    
    # Update Z -------------------------------------
    Z <- updateZ(m, Theta_new, Y, D, n_cores = n_cores) 
    #Z_old <- Z_new
    #Z_new <- Temp 
    
    # Update Y -------------------------------------
    Y <- updateY(Theta_new, Z, Y) 
    #Y_new <- Y_old
    #Y_new <- Temp
    
    # Check whether the algorithm is ready ----------
    difference <- relative_difference_precision_matrices(Theta_new, Theta_old, n_cores = 1)
    print(difference)
    if (difference < epsilon) { 
      converged <- TRUE
      iter <- iter + 1
      break() 
    }
    
    if (iter >= maxiter) { 
      warning("Maximum number of iterations reached. Stopping criterion not met")
      break()
    }
  }
  
  
  res <- list(
    data = data, 
    W = W, 
    lambda1 = lambda1,
    lambda2 = lambda2,
    rho = rho, 
    converged = converged,
    n_iterations = iter, 
    m = m, 
    p = p, 
    n_obs = n_obs
  )
  
  class(res) <- "CVN"
  
  return(res)
}
