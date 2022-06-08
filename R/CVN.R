#' Inferring a Covariate-Varying Network (CVN)
#' 
#' Estimates a covariate-varying network model (CVN), i.e., \eqn{m} 
#' Gaussian graphical models that change with (multiple) external covariate(s). 
#' The smoothing between the graphs is specified by the \eqn{(m x m)}-dimensional
#' weight matrix \eqn{W}. The function returns the estimated precision matrices 
#' for each graph. 
#' 
#' @param data A list with matrices. The number of columns should be the 
#'                 same for each matrix. Number of observations can differ
#' @param W The \code{(m x m)}-dimensional uppertriangular weight matrix
#' @param lambda1 The \eqn{\lambda1} LASSO penalty term 
#' @param lambda2 The \eqn{\lambda2} global smoothing parameter 
#' @param rho The \eqn{\rho} ADMM's penalty parameter (Default = 1)
#' @param epsilon If the relative difference between two update steps is 
#'                smaller than \eqn{\epsilon}, the algorithm stops
#'                (Default: \code{1e-5})
#' @param maxiter Maximum number of iterations (Default = 100)
#' @param n_cores Number of cores used (Default = 1)
#' @param truncate All values of the final Theta below \code{truncate} will be 
#'                 set to \code{0}. (Default = \code{1e-5})
#' @param normalized Data is normalized if \code{TRUE}. Otherwise the data is only
#'                   centered (Default: \code{FALSE})
#' @param warmstart If \code{TRUE}, use the \code{huge} package for estimating
#'                  the individual graphs first (Default: \code{FALSE})
#' @param verbose Verbose (Default: \code{FALSE}) 
#' 
#' @return A CVN object; a list with entries
#'    \item{\code{Theta}}{The estimated precision matrices}
#'    \item{\code{adj_matrices}}{The estimated adjacency matrices; 1 if there is an edge, 0 otherwise}
#'    \item{\code{Sigma}}{Empirical covariance matrices}
#'    \item{\code{m}}{Number of graphs}
#'    \item{\code{p}}{Number of variables}
#'    \item{\code{n_obs}}{Vector of length \eqn{m} with number of observations for each graph}
#'   \item{\code{data}}{The \code{data}, but then normalized or centered}
#'   \item{\code{normalized}}{If \code{TRUE}, \code{data} was normalized. Otherwise \code{data} was only centered}
#'   \item{\code{W}}{The \code{m x m} weight matrix}
#'   \item{\code{D}}{Matrix \eqn{D} used for the Generalized LASSO}
#'   \item{\code{lambda1}}{The \eqn{\lambda1} LASSO penalty term}
#'   \item{\code{lambda2}}{The \eqn{\lambda2} global smoothing parameter} 
#'   \item{\code{rho}}{The \eqn{\rho} ADMM's penalty parameter} 
#'   \item{\code{epsilon}}{The stopping criterion \eqn{\epsilon}} 
#'   \item{\code{converged}}{If \code{TRUE}, stopping condition has been met}
#'   \item{\code{n_iterations}}{Number of iterations}
#'   \item{\code{truncate}}{Truncation value for \eqn{\Theta}}
#'   
#' @export
# TODO: change lambda1 and lambda2 to a grid!
CVN <- function(data, W, lambda1 = 1, lambda2 = 1, 
                rho = 1,
                epsilon = 10^(-5),
                maxiter = 100, 
                n_cores = 1, 
                truncate = 1e-5, 
                normalized = FALSE, 
                warmstart = FALSE, 
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
  Sigma <- lapply(1:m, function(i) cov(data[[i]])*(n_obs[i] - 1) / n_obs[i]) 
 
  # Initialize variables for the algorithm -----------------
  # Generate matrix D for the generalized LASSO 
  D <- CVN::create_matrix_D(W, lambda1, lambda2, rho)
  
  # Generate matrices for ADMM iterations
  Theta_old <- rep(list(diag(p)), m) # m (p x p)-dimensional identity matrices
  Theta_new <- lapply(Sigma, function(S) diag(1/diag(S)))
  Z <- rep(list(matrix(0, nrow = p, ncol = p)), m) # m (p x p)-dimensional zero matrices
  Y <- rep(list(matrix(0, nrow = p, ncol = p)), m)
  
  # if warmstart, the individual graphs are estimated
  # separately with huge
  if (warmstart) { 
    Theta_old <- Theta_new
    
    # estimate the graph with the GLASSO. The GLASSO is the only
    # option here when we want to use huge.select for 
    # selecting the optimal lambda value
    Theta_new <- lapply(data, function(X) {
        est <- huge(X, method = "glasso", verbose = FALSE)
        huge.select(est, verbose = FALSE)$opt.icov
      })
  } 
  
  # Initialize a Temp variable for Theta
  Temp <- rep(list(matrix(0, nrow = p, ncol = p)), m)
  
  # keep track whether the algorithm finished, either by 
  # whether the stopping condition was met, or the maximum
  # number of iterations was reached
  converged <- FALSE 
  iter <- 1 
  
  repeat{
    
    # Update Theta ---------------------------------
    Temp <- updateTheta(m, Z, Y, Sigma, n_obs, rho, n_cores = n_cores)
    Theta_old <- Theta_new 
    Theta_new <- Temp 
    
    # Update Z -------------------------------------
    Z <- updateZ(m, p, Theta_new, Y, D, n_cores = n_cores) 
    
    # Update Y -------------------------------------
    Y <- updateY(Theta_new, Z, Y) 
    
    # Check whether the algorithm is ready ----------
    difference <- relative_difference_precision_matrices(Theta_new, Theta_old, n_cores = 1)
    
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
  
  Z <- lapply(res$Theta, function(z) { 
      ifelse(abs(z) <= 1e-5, 0, z)
    })
  

  adj_matrices <- lapply(Z, function(X) { 
    diag(X) <- 0 
    Matrix( as.numeric( abs(X) >= 2*.Machine$double.eps), ncol = ncol(X) , sparse = TRUE)
  })
                       
  res <- list(
    Theta = Z,
    adj_matrices = adj_matrices, 
    Sigma = Sigma,
    m = m, 
    p = p, 
    n_obs = n_obs, 
    data = data, 
    normalized = normalized,
    W = W, 
    D = D, 
    lambda1 = lambda1,
    lambda2 = lambda2,
    rho = rho, 
    epsilon = epsilon,
    converged = converged,
    n_iterations = iter, 
    truncate = truncate
  )
  
  class(res) <- "CVN"
  
  return(res)
}
