#' Estimating a Covariate-Varying Network (CVN)
#' 
#' Estimates a covariate-varying network model (CVN), i.e., \eqn{m} 
#' Gaussian graphical models that change with (multiple) external covariate(s). 
#' The smoothing between the graphs is specified by the \eqn{(m \times m)}-dimensional
#' weight matrix \eqn{W}. The function returns the estimated precision matrices 
#' for each graph. 
#' @section Reusing Estimates: When estimating the graph for different values of 
#' \eqn{\lambda_1} and \eqn{\lambda_2}, we use the graph estimated (if available) 
#' for other \eqn{\lambda_1} and \eqn{\lambda_2} values closest to them. 
#' 
#' @param data A list with matrices, each entry associated with a single graph. 
#'             The number of columns should be the same for each matrix. 
#'             Number of observations can differ
#' @param W The \eqn{(m \times m)}-dimensional symmetric 
#'          weight matrix \eqn{W}
#' @param lambda1 Vector with different \eqn{\lambda_1} LASSO penalty terms 
#'                (Default: \code{1:10})
#' @param lambda2 Vector with different \eqn{\lambda_2} global smoothing parameter values 
#'                (Default: \code{1:10})
#' @param rho The \eqn{\rho} penalty parameter for the global ADMM algorithm (Default: \code{1})
#' @param eps If the relative difference between two update steps is 
#'                smaller than \eqn{\epsilon}, the algorithm stops. 
#'                See \code{\link{relative_difference_precision_matrices}}
#'                (Default: \code{1e-5})
#' @param maxiter Maximum number of iterations (Default: \code{100})
#' @param truncate All values of the final \eqn{\hat{\Theta}_i}'s below \code{truncate} will be 
#'                 set to \code{0}. (Default: \code{1e-5})
#' @param rho_genlasso The \eqn{\rho} penalty parameter for the ADMM algorithm 
#'                used to solve the generalized LASSO (Default: \code{1})
#' @param eps_genlasso If the relative difference between two update steps is 
#'                smaller than \eqn{\epsilon}, the algorithm stops. 
#'                (Default: \code{1e-10})
#' @param maxiter_genlasso Maximum number of iterations for solving 
#'                the generalized LASSO problem (Default: \code{100})
#' @param truncate_genlasso All values of the final \eqn{\hat{\beta}} below 
#'                 \code{truncate_genlasso} will be set to \code{0}. 
#'                 (Default: \code{1e-4})
#' @param n_cores Number of cores used (Default: \code{1})
#' @param normalized Data is normalized if \code{TRUE}. Otherwise the data is only
#'                   centered (Default: \code{FALSE})
#' @param warmstart If \code{TRUE}, use the \code{\link[huge]{huge}} package for estimating
#'                  the individual graphs first (Default: \code{FALSE})
#' @param use_previous_estimate If \code{TRUE}, the estimated graph found for the previous 
#'                  values of \eqn{\lambda_1} and \eqn{\lambda_2} is used as starting point
#'                  for the next estimate (Default: \code{TRUE})
#' @param use_genlasso If \code{TRUE}, use the \code{genlasso} package in 
#'                  the \eqn{Z}-update step, rather then the ADMM (Default: \code{FALSE})
#' @param verbose Verbose (Default: \code{FALSE}) 
#' 
#' @return A \code{CVN} object containing the estimates for all the graphs 
#'    for each different value of \eqn{(\lambda_1, \lambda_2)}. General results for 
#'    the different values of \eqn{(\lambda_1, \lambda_2)} can be found in the data frame
#'    \code{results}. It consists of multiple columns, namely: 
#'    \item{\code{lambda1}}{\eqn{\lambda_1} value}
#'    \item{\code{lambda2}}{\eqn{\lambda_1} value}
#'    \item{\code{converged}}{whether algorithm converged or not}
#'    \item{\code{value}}{value of the negative log-likelihood function}
#'    \item{\code{n_iterations}}{number of iterations of the ADMM}
#'    \item{\code{aic}}{Aikake information criteration}
#'    \item{\code{id}}{The id. This corresponds to the indices of the lists}
#'    In addition, the results of for the different values of \eqn{(\lambda_1, \lambda_2)}
#'    can be found in the following lists: 
#'    \item{\code{Theta}}{The estimated precision matrices \eqn{\{ \hat{\Theta}_i \}_{i = 1}^m}}
#'    \item{\code{adj_matrices}}{The estimated adjacency matrices; 
#'                         \code{1} if there is an edge, \code{0} otherwise. 
#'                         Matrices are sparse using package \code{\link[Matrix]{Matrix}}}
#'    \item{\code{Sigma}}{Empirical covariance matrices \eqn{\{\hat{\Sigma}_i\}_{i = 1}^m}}
#'    \item{\code{m}}{Number of graphs}
#'    \item{\code{p}}{Number of variables}
#'    \item{\code{n_obs}}{Vector of length \eqn{m} with number of observations for each graph}
#'   \item{\code{data}}{The \code{data}, but then normalized or centered}
#'   \item{\code{normalized}}{If \code{TRUE}, \code{data} was normalized. Otherwise \code{data} was only centered}
#'   \item{\code{W}}{The \eqn{(m \times m)}-dimensional weight matrix \eqn{W}}
#'   \item{\code{D}}{Matrix \eqn{D} used for the Generalized LASSO}
#'   \item{\code{lambda1}}{The \eqn{\lambda_1} LASSO penalty term}
#'   \item{\code{lambda2}}{The \eqn{\lambda_2} global smoothing parameter} 
#'   \item{\code{rho}}{The \eqn{\rho} ADMM's penalty parameter} 
#'   \item{\code{epsilon}}{The stopping criterion \eqn{\epsilon}} 
#'   \item{\code{converged}}{If \code{TRUE}, stopping condition has been met; \code{FALSE} otherwise}
#'   \item{\code{value}}{The final relative difference}
#'   \item{\code{n_iterations}}{Number of iterations}
#'   \item{\code{truncate}}{Truncation value for \eqn{\{ \hat{\Theta}_i \}_{i = 1}^m}}
#'   
#' @examples 
#' data(grid)
#' m <- 9 # must be 9 for this example
#' 
#' #' Choice of the weight matrix W. 
#' #' (uniform random) 
#' W <- matrix(runif(m*m), ncol = m)
#' W <- W %*% t(W)
#' W <- W / max(W)
#' diag(W) <- 0
#' 
#' # lambdas:
#' lambda1 = 1:2
#' lambda2 = 1:2
#' 
#' (cvn <- CVN::CVN(grid, W, lambda1 = lambda1, lambda2 = lambda2, eps = 1e-3, maxiter = 1000, verbose = TRUE))
#' @export
CVN <- function(data, W, lambda1 = 1:10, lambda2 = 1:10, 
                rho = 1,
                rho_genlasso = 1,
                eps = 1e-5,
                maxiter = 100, 
                truncate = 1e-5, 
                eps_genlasso = 1e-10, 
                maxiter_genlasso = 100, 
                truncate_genlasso = 1e-4, 
                n_cores = 1, 
                normalized = FALSE, 
                warmstart = FALSE, 
                use_previous_estimate = FALSE,
                use_genlasso_package = FALSE, 
                verbose = FALSE) { 
  
  # Check correctness input -------------------------------
  CVN::check_correctness_input(data, W, lambda1, lambda2, rho)
  
  if (use_previous_estimate) { 
    warning("Using the previous estimate might not always work. We are working on a test") 
  }
  
  # Extract variables -------------------------------------
  m <- length(data)       # total number of graphs  
  p <- ncol(data[[1]])    # total number of variables
  n_obs <- sapply(data, function(X) nrow(X))    # no. of observations per graph 
  
  # Center or normalize data -------------------------------
  data <- lapply(data, function(X) scale(X, scale = normalized))
  
  # Compute the empirical covariance matrices --------------
  Sigma <- lapply(1:m, function(i) cov(data[[i]])*(n_obs[i] - 1) / n_obs[i])
  
  # Create initial values for the ADMM algorithm ----------
  Z <- rep(list(matrix(0, nrow = p, ncol = p)), m) # m (p x p)-dimensional zero matrices
  Y <- rep(list(matrix(0, nrow = p, ncol = p)), m)
  
  # if warmstart, the individual graphs are first estimated using the GLASSO.
  if (warmstart) { 
    
    if (verbose) { 
      cat("Warm start...\n") 
    }
    
    # We use the first lambda1 value in the vector
    Theta <- lapply(Sigma, function(S) { 
       est <- glasso::glasso(s = S, rho = lambda1[1])
       est$w
    })
  }  else { 
    Theta <- lapply(Sigma, function(S) diag(1/diag(S)))
  }
  
  # initialize results list ------------
  global_res <- list(
    Theta    = list(),          
    adj_matrices = list(),   
    Sigma    = Sigma,
    m        = m, 
    p        = p, 
    n_obs    = n_obs,
    data     = data,
    normalized = normalized,
    W        = W, 
    rho      = rho, 
    eps      = eps,
    truncate = truncate
  )
  
  # data frame with the results for each unique (lambda1,lambda2) pair
  res <- data.frame(expand.grid(lambda1 = lambda1, 
                                lambda2 = lambda2, 
                                converged = FALSE, 
                                value = NA, 
                                n_iterations = NA, 
                                aic = NA))
  res$id <- 1:nrow(res)
  
  for (i in 1:nrow(res)) { 
    
    if (verbose) { 
      cat(sprintf("\n(%.0f%%) lambda1: %g / lambda2: %g\n", (i-1)/nrow(res)*100, res$lambda1[i], res$lambda2[i])) 
    }
    
    # Initialize variables for the algorithm -----------------
    # Generate matrix D for the generalized LASSO 
    D <- CVN::create_matrix_D(W, res$lambda1[i], res$lambda2[i], rho)
    
    a <- CVN::matrix_A_inner_ADMM(m, D) + 1
    
    # Estimate the graphs -------------------------------------
    eta1 <- res$lambda1[i] / rho 
    eta2 <- res$lambda2[i] / rho 
    est <- CVN::estimate(m, p, nrow(D), Theta, Z, Y, a, eta1, eta2, Sigma, n_obs, 
                         rho, rho_genlasso, 
                         eps, eps_genlasso, 
                         maxiter, maxiter_genlasso, truncate = truncate, 
                         truncate_genlasso = truncate_genlasso, 
                         use_genlasso_package = use_genlasso_package, verbose = verbose)  
    
    # Process results -----------------------------------------
    global_res$Theta[[i]] <- est$Z
    global_res$adj_matrices[[i]] <- est$adj_matrices 
    
    res$converged[i]    <- est$converged
    res$value[i]        <- est$value
    res$n_iterations[i] <- est$n_iterations
    
    # determine the AIC
    res$aic[i] <- CVN::AIC(Theta = est$Z, 
                           adj_matrices = est$adj_matrices, 
                           Sigma = Sigma, 
                           n_obs = n_obs) 
    
    if (!use_previous_estimate) { 
      Theta <- lapply(Sigma, function(S) diag(1/diag(S)))
      Z     <- rep(list(matrix(0, nrow = p, ncol = p)), m) # m (p x p)-dimensional zero matrices
      Y     <- rep(list(matrix(0, nrow = p, ncol = p)), m)
    } else { 
      Theta <- est$Theta
      Z     <- est$Z
      Y     <- est$Y
    }
  }
  
  if (verbose) { 
    cat(sprintf("\n(100%%)\n")) 
  }
     
  # Collect all the results & input ---------------------------
  global_res$results  <- res                  
  
  #class(global_res) <- "CVN" # TODO
  
  return(global_res)
}

#' Print Function for the CVN Object Class
#'
#' @export
print.CVN <- function(cvn, ...) { 
  cat(sprintf("Covariate-varying Network (CVN)\n\n"))
  
  if (all(cvn$converged)) {
    cat(green(sprintf("\u2713 CONVERGED\n\n")))
  } else { 
    cat(red(sprintf("\u2717 DID NOT CONVERGE\n\n"))) 
  }
  
  cat(sprintf("   -- No. of graphs (m):    %d\n", cvn$m))
  cat(sprintf("   -- No. of variables (p): %d\n\n", cvn$p))
  print(cvn2$results)
}

#' Plot Function for CVN Object Class
#' 
#' @export
plot.CVN <- function(cvn, ...) { 
  cat("TO DO, see https://kateto.net/network-visualization\n")
}
