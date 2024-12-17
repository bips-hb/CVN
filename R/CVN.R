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
#' @param lambda1 Vector with different \eqn{\lambda_1}. LASSO penalty terms 
#'                (Default: \code{1:2})
#' @param lambda2 Vector with different \eqn{\lambda_2}. The global smoothing parameter values 
#'                (Default: \code{1:2})
#' @param gamma1 A vector of \eqn{\gamma_1}'s LASSO penalty terms, where
#'              \eqn{\gamma_1 = \frac{2 \lambda_1}{m p (1 - p)}}. If \code{gamma1}
#'              is set, the value of \code{lambda1} is ignored. (Default: \code{NULL}).
#' @param gamma2 A vector of \eqn{\gamma_2}'s global smoothing parameters, where
#'               that \eqn{\gamma_2 = \frac{4 \lambda_2}{m(m-1)p(p-1)}}. If \code{gamma2}
#'               is set, the value of \code{lambda2} is ignored.(Default: \code{NULL}).
#' @param rho The \eqn{\rho} penalty parameter for the global ADMM algorithm (Default: \code{1})
#' @param eps If the relative difference between two update steps is
#'                smaller than \eqn{\epsilon}, the algorithm stops.
#'                See \code{\link{relative_difference_precision_matrices}}
#'                (Default: \code{1e-4})
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
#' @param n_cores Number of cores used (Default: max. number of cores - 1, or
#'                the total number penalty term pairs if that is less)
#' @param normalized Data is normalized if \code{TRUE}. Otherwise the data is only
#'                   centered (Default: \code{FALSE})
#' @param warmstart If \code{TRUE}, use the \code{\link[huge]{huge}} package for estimating
#'                  the individual graphs first (Default: \code{TRUE})
#' @param minimal If \code{TRUE}, the returned \code{cvn} is minimal in terms of
#'                  memory, i.e., \code{Theta}, \code{data} and \code{Sigma} are not
#'                  returned (Default: \code{FALSE})
#' @param gamma_ebic Gamma value for the eBIC (Default: 0.5)                  
#' @param verbose Verbose (Default: \code{TRUE})
#'
#' @return A \code{CVN} object containing the estimates for all the graphs
#'    for each different value of \eqn{(\lambda_1, \lambda_2)}. General results for
#'    the different values of \eqn{(\lambda_1, \lambda_2)} can be found in the data frame
#'    \code{results}. It consists of multiple columns, namely:
#'    \item{\code{id}}{The id. This corresponds to the indices of the lists}
#'    \item{\code{lambda1}}{\eqn{\lambda_1} value}
#'    \item{\code{lambda2}}{\eqn{\lambda_2} value}
#'    \item{\code{gamma1}}{\eqn{\gamma_1} value}
#'    \item{\code{gamma2}}{\eqn{\gamma_2} value}    
#'    \item{\code{converged}}{whether algorithm converged or not}
#'    \item{\code{value}}{value of the negative log-likelihood function}
#'    \item{\code{n_iterations}}{number of iterations of the ADMM}
#'    \item{\code{aic}}{Aikake information criterion}
#'    \item{\code{bic}}{Bayesian information criterion}
#'    \item{\code{ebic}}{Extended Bayesian information criterion}
#'    \item{\code{sparsity}}{Median sparsity of the m estimated graphs}
#'    The estimates of the precision matrices and the corresponding adjacency matrices
#'    for the different values of \eqn{(\lambda_1, \lambda_2)} can be found
#'    \item{\code{Theta}}{A list with the estimated precision matrices \eqn{\{ \hat{\Theta}_i(\lambda_1, \lambda_2) \}_{i = 1}^m},
#'                        (only if \code{minimal = FALSE})}
#'    \item{\code{adj_matrices}}{A list with the estimated adjacency matrices corresponding to the
#'                               estimated precision matrices in \code{Theta}. The entries
#'                               are \code{1} if there is an edge, \code{0} otherwise.
#'                               The matrices are sparse using package \code{\link[Matrix]{Matrix}}}
#'    In addition, the input given to the CVN function is stored in the object as well:
#'    \item{\code{Sigma}}{Empirical covariance matrices \eqn{\{\hat{\Sigma}_i\}_{i = 1}^m},
#'                              (only if \code{minimal = FALSE})}
#'    \item{\code{m}}{Number of graphs}
#'    \item{\code{p}}{Number of variables}
#'    \item{\code{n_obs}}{Vector of length \eqn{m} with number of observations for each graph}
#'   \item{\code{data}}{The \code{data}, but then normalized or centered (only if \code{minimal = FALSE})}
#'   \item{\code{W}}{The \eqn{(m \times m)}-dimensional weight matrix \eqn{W}}
#'   \item{\code{maxiter}}{Maximum number of iterations for the ADMM}
#'   \item{\code{rho}}{The \eqn{\rho} ADMM's penalty parameter}
#'   \item{\code{eps}}{The stopping criterion \eqn{\epsilon}}
#'   \item{\code{truncate}}{Truncation value for \eqn{\{ \hat{\Theta}_i \}_{i = 1}^m}}
#'   \item{\code{maxiter_genlasso}}{Maximum number of iterations for the generarlzed LASSO}
#'   \item{\code{rho_genlasso}}{The \eqn{\rho} generalized LASSO penalty parameter}
#'   \item{\code{eps_genlasso}}{The stopping criterion \eqn{\epsilon} for the generalized LASSO}
#'   \item{\code{truncate_genlasso}}{Truncation value for \eqn{\beta} of the generalized LASSO}
#'   \item{\code{n_lambda_values}}{Total number of \eqn{(\lambda_1, \lambda_2)} value combinations}
#'   \item{\code{normalized}}{If \code{TRUE}, \code{data} was normalized. Otherwise \code{data} was only centered}
#'   \item{\code{warmstart}}{If \code{TRUE}, warmstart was used}
#'   \item{\code{minimal}}{If \code{TRUE}, \code{data}, \code{Theta} and \code{Sigma} are not added}
#'   \item{\code{hits_border_aic}}{If \code{TRUE}, the optimal model based on the AIC hits the border of \eqn{(\lambda_1, \lambda_2)}}
#'   \item{\code{hits_border_bic}}{If \code{TRUE}, the optimal model based on the BIC hits the border of \eqn{(\lambda_1, \lambda_2)}}
#'   
#' @name CVN function
#' @aliases CVN
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
#' lambda1 = 1  # can also be lambda1 = 1:2 
#' lambda2 = 1
#' 
#' (cvn <- CVN::CVN(grid, W, lambda1 = lambda1, lambda2 = lambda2, 
#'                  eps = 1e-3, maxiter = 1000, verbose = TRUE))
#' @export
CVN <- function(data, W, lambda1 = 1:2, lambda2 = 1:2,
                gamma1 = NULL, gamma2 = NULL,
                rho = 1,
                eps = 1e-4,
                maxiter = 100,
                truncate = 1e-5,
                rho_genlasso = 1,
                eps_genlasso = 1e-10,
                maxiter_genlasso = 100,
                truncate_genlasso = 1e-4,
                n_cores = min(length(lambda1)*length(lambda2), parallel::detectCores() - 1),
                normalized = FALSE,
                warmstart = TRUE,
                minimal = FALSE,
                gamma_ebic = 0.5,
                verbose = TRUE) {

  # Check correctness input -------------------------------
  check_correctness_input(data, W, lambda1, lambda2, gamma1, gamma2, rho)

  # Extract variables -------------------------------------
  m <- length(data)       # total number of graphs
  p <- ncol(data[[1]])    # total number of variables
  n_obs <- sapply(data, function(X) nrow(X))    # no. of observations per graph

  # convert the lambda values to gamma values or the other way around
  if (is.null(gamma1) && is.null(gamma2)) {
    gamma1 <- 2*lambda1 / (m*p*(p-1))
    gamma2 <- 4*lambda2 / (m*(m-1)*p*(p-1))
  } else {
    lambda1 <- (gamma1*m*p*(p-1)) / 2
    lambda2 <- (gamma2*m*(m-1)*p*(p-1)) / 4
  }

  # When the weight matrix is completely zero, there is no smoothing
  # between graphs. Therefore, the value of lambda2 is irrelevant. 
  # In this case, we fix lambda2 to 0 and inform the user of this fact
  if (sum(W) == 0) { 
    warning("Since weight matrix W is zero, there is no smoothing. lambda2 is fixed to 0")
    lambda2 <- 0
  }

  # Set-up cluster ---------------------------
  if (n_cores > 1) {
    cl <- makeCluster(n_cores)
    #registerDoSNOW(cl)              # 16Dec
    registerDoParallel(cl)           # 16Dec
    
    # Export CVN library to workers  # 16Dec
    clusterEvalQ(cl, library(CVN))
    
    opts <- list(progress = function(i) {i}) # empty place holder for foreach
  }

  if (verbose) {
    cat(sprintf("Estimating a CVN with %d graphs...\n\n", m))
    cat(sprintf("Number of cores: %d\n", n_cores))
    if (warmstart) {
      cat("Uses a warmstart...\n\n")
    } else {
      cat("No warmstart...\n\n")
    }

    # Set-up a progress bar ---------------------------------
    if (n_cores > 1) {
      pb <- progress::progress_bar$new(
        format = "estimating CVN [:bar] :percent eta: :eta",
        total = length(lambda1)*length(lambda2) + 1, clear = FALSE, width= 80, show_after = 0)
      pb$tick()

      progress <- function(i) pb$tick()
      opts     <- list(progress = progress) # used by doSNOW
    }
  }

  # Center or normalize data -------------------------------
  data <- lapply(data, function(X) scale(X, scale = normalized))

  # Compute the empirical covariance matrices --------------
  Sigma <- lapply(1:m, function(i) cov(data[[i]])*(n_obs[i] - 1) / n_obs[i])


  # initialize results list ------------
  global_res <- list(
    Theta    = list(),
    adj_matrices = list(),
    Sigma    = Sigma,
    m        = m,
    p        = p,
    n_obs    = n_obs,
    data     = data,
    W        = W,
    maxiter  = maxiter,
    rho      = rho,
    eps      = eps,
    truncate = truncate,
    rho_genlasso      = rho_genlasso,
    eps_genlasso      = eps_genlasso,
    maxiter_genlasso  = maxiter_genlasso,
    truncate_genlasso = truncate_genlasso,
    n_lambda_values   = length(lambda1) * length(lambda2),
    normalized = normalized,
    warmstart  = warmstart,
    minimal = minimal
  )

  # data frame with the results for each unique (lambda1,lambda2) pair
  res <- data.frame(id = NA,
                    expand.grid(lambda1 = lambda1,
                                lambda2 = lambda2,
                                gamma1 = NA,
                                gamma2 = NA,
                                converged = FALSE,
                                value = NA,
                                n_iterations = NA,
                                aic = NA,
                                bic = NA,
                                ebic = NA,
                                sparsity = NA,
                                sparsity_iqr = NA))
  res$gamma1 <- 2*res$lambda1 / (m*p*(p-1))
  res$gamma2 <- 4*res$lambda2 / (m*(m-1)*p*(p-1))
  res$id <- 1:nrow(res)

  # estimate the graphs for the different values of (lambda1, lambda2) --------
  estimate_lambda_values <- function(i) {
    # Create initial values for the ADMM algorithm ----------
    Z <- rep(list(matrix(0, nrow = p, ncol = p)), m) # m (p x p)-dimensional zero matrices
    Y <- rep(list(matrix(0, nrow = p, ncol = p)), m)

    # if warmstart, the individual graphs are first estimated using the GLASSO.
    if (warmstart) {
      # We use the lambda1 value
      Theta <- lapply(Sigma, function(S) {
        est <- glasso::glasso(s = S, rho = res$lambda1[i])
        est$w
      })
    }  else {
      Theta <- lapply(Sigma, function(S) diag(1/diag(S)))
    }

    # Initialize variables for the algorithm -----------------
    if (sum(W) == 0) { # if the weight matrix is empty, the value a = eta_1^2 + 1
      a <- (res$lambda1[i] / rho)^2 + 1
    } else {
      # Determine the value of the diagonal matrix A such that A - D'D > 0 (positive definite)
      a <- matrix_A_inner_ADMM(W, res$lambda1[i] / rho, res$lambda2[i] / rho) + 1
    }

    # Estimate the graphs -------------------------------------
    eta1 <- res$lambda1[i] / rho
    eta2 <- res$lambda2[i] / rho
         estimate(m, p, W, Theta, Z, Y, a, eta1, eta2, Sigma, n_obs,
                  rho, rho_genlasso,
                  eps, eps_genlasso,
                  maxiter, maxiter_genlasso, truncate = truncate,
                  truncate_genlasso = truncate_genlasso,
                  verbose = verbose)
  }

  # go over each pair of penalty terms
  if (n_cores > 1) { # parallel
    # est <- foreach(i = 1:(length(lambda1)*length(lambda2)),
    #                .options.snow = opts,
    #                .export = c("estimate")) %dopar% {
    #                  estimate_lambda_values(i)
    #                }
 
    est <- foreach(i = 1:(length(lambda1)*length(lambda2)),
            .multicombine = TRUE,
            .export = c("estimate")) %dopar% {
             estimate_lambda_values(i)
            }
  } else { # sequential
    est <- lapply(1:(length(lambda1) * length(lambda2)), function(i) {
      estimate_lambda_values(i)
    })
  }

  

  # Process results -----------------------------------------
  for (i in 1:(length(lambda1)*length(lambda2))) {
    global_res$Theta[[i]] <- est[[i]]$Z
    global_res$adj_matrices[[i]] <- est[[i]]$adj_matrices

    res$converged[i]    <- est[[i]]$converged
    res$value[i]        <- est[[i]]$value
    res$n_iterations[i] <- est[[i]]$n_iterations
    

    ic <- determine_information_criterion(Theta = est[[i]]$Z,
                                          adj_matrices = est[[i]]$adj_matrices,
                                          Sigma = Sigma,
                                          n_obs = n_obs,
                                          gamma = gamma_ebic)
    res$aic[i] <- ic$aic
    res$bic[i] <- ic$bic
    res$ebic[i] <- ic$ebic
    
    # Median sparsity
    res$sparsity[i] <- median(mapply(sum, est[[i]]$adj_matrices) / 2)
    res$sparsity_iqr[i] <- IQR(mapply(sum, est[[i]]$adj_matrices) / 2)

  }

  # stop the progress bar
  if (n_cores > 1 && verbose) {
    pb$terminate()
  }

  # stop the cluster
  if (n_cores > 1) {
    stopCluster(cl)
  }

  # determine whether the optimal model based on the AIC or BIC hits the border
  # of lambda1 and/or lambda2
  hit <- hits_end_lambda_intervals(res)
  global_res$hits_border_aic <- hit$hits_border_aic
  global_res$hits_border_bic <- hit$hits_border_bic

  # Collect all the results & input ---------------------------
  global_res$results  <- res

  class(global_res) <- c("cvn", "list")

  if (minimal) {
    global_res <- strip_cvn(global_res)
  }

  return(global_res)
}


