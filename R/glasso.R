#' Estimating Multiple Networks Separately 
#' 
#' A wrapper for the GLASSO in the context of CVNs. Each graph 
#' is estimated individually. There is NO smoothing between the graphs. 
#' This function relies completely on the \code{\link{glasso}} package. 
#' The output is, therefore, slightly different than for the 
#' \code{\link{CVN}} function. 
#' 
#' @param data A list with matrices, each entry associated with a single graph. 
#'             The number of columns should be the same for each matrix. 
#'             Number of observations can differ
#' @param lambda1 Vector with different \eqn{\lambda_1} LASSO penalty terms 
#'                (Default: \code{1:2})
#' @param eps Threshold for convergence (Default: \code{1e-4}; the same as in the 
#'            \code{glasso} package)
#' @param maxiter Maximum number of iterations (Default: 10,000)
#' @param n_cores Number of cores used (Default: max. number of cores - 1, or 
#'                the total number penalty term pairs if that is less)
#' @param normalized Data is normalized if \code{TRUE}. Otherwise the data is only
#'                   centered (Default: \code{FALSE})
#' @param verbose Verbose (Default: \code{TRUE}) 
#' 
#' @return A \code{CVN} object containing the estimates for all the graphs 
#'    for different value of \eqn{\lambda_1}. General results for 
#'    the different value of \eqn{\lambda_1} can be found in the data frame
#'    \code{results}. It consists of multiple columns, namely: 
#'    \item{\code{lambda1}}{\eqn{\lambda_1} value}
#'    \item{\code{value}}{value of the negative log-likelihood function}
#'    \item{\code{aic}}{Aikake information criteration}
#'    \item{\code{id}}{The id. This corresponds to the indices of the lists}
#'    The estimates of the precision matrices and the corresponding adjacency matrices
#'    for the different values of \eqn{\lambda_1} can be found 
#'    \item{\code{Theta}}{A list with the estimated precision matrices \eqn{\{ \hat{\Theta}_i(\lambda_1) \}_{i = 1}^m}}
#'    \item{\code{adj_matrices}}{A list with the estimated adjacency matrices corresponding to the 
#'                               estimated precision matrices in \code{Theta}. The entries 
#'                               are \code{1} if there is an edge, \code{0} otherwise. 
#'                               The matrices are sparse using package \code{\link[Matrix]{Matrix}}}
#'    In addition, the input given to this function is stored in the object as well:
#'    \item{\code{Sigma}}{Empirical covariance matrices \eqn{\{\hat{\Sigma}_i\}_{i = 1}^m}}
#'    \item{\code{m}}{Number of graphs}
#'    \item{\code{p}}{Number of variables}
#'    \item{\code{n_obs}}{Vector of length \eqn{m} with number of observations for each graph}
#'   \item{\code{data}}{The \code{data}, but then normalized or centered}
#'   \item{\code{maxiter}}{Maximum number of iterations for the ADMM}
#'   \item{\code{eps}}{The stopping criterion \eqn{\epsilon}} 
#'   \item{\code{n_lambda_values}}{Total number of \eqn{\lambda_1} values}
#'   \item{\code{normalized}}{If \code{TRUE}, \code{data} was normalized. Otherwise \code{data} was only centered}
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
#' lambda1 = 1:4
#' 
#' (glasso_est <- CVN::glasso(grid, lambda1 = lambda1))
#' @export
glasso <- function(data, lambda1 = 1:2, 
                eps = 1e-4,
                maxiter = 10000, 
                n_cores = min(length(lambda1), parallel::detectCores() - 1), 
                normalized = FALSE, 
                verbose = TRUE) { 
  
  # Extract variables -------------------------------------
  m <- length(data)       # total number of graphs  
  p <- ncol(data[[1]])    # total number of variables
  n_obs <- sapply(data, function(X) nrow(X))    # no. of observations per graph 
  
  # Set-up cluster ---------------------------
  cl <- makeCluster(n_cores)
  registerDoSNOW(cl)
  opts <- list(progress = function(i) {i}) # empty place holder for foreach
  
  if (verbose) { 
    cat(sprintf("Estimating %d individual graphs WITHOUT smoothing...\n\n", m))
    cat(sprintf("Number of cores: %d\n", n_cores))
    
    # Set-up a progress bar ---------------------------------
    pb <- progress::progress_bar$new(
      format = "estimating glasso [:bar] :percent eta: :eta",
      total = length(lambda1) + 1, clear = FALSE, width= 80, show_after = 0)
    pb$tick()
    
    #pb       <- txtProgressBar(max = length(lambda1)*length(lambda2), style = 3)
    progress <- function(i) pb$tick()
    opts     <- list(progress = progress) # used by doSNOW
  }
  
  # Center or normalize data -------------------------------
  data <- lapply(data, function(X) scale(X, scale = normalized))
  
  # Compute the empirical covariance matrices --------------
  Sigma <- lapply(1:m, function(i) cov(data[[i]])*(n_obs[i] - 1) / n_obs[i])
  
  # initialize results list ------------
  global_res <- list(
    Theta    = lapply(1:length(lambda1), function(i) list()),          
    adj_matrices = lapply(1:length(lambda1), function(i) list()),   
    Sigma    = Sigma,
    m        = m, 
    p        = p, 
    n_obs    = n_obs,
    data     = data,
    maxiter  = maxiter, 
    eps      = eps,
    n_lambda_values = length(lambda1), 
    normalized = normalized
  )
  
  # data frame with the results for each unique lambda1 value
  res <- data.frame(expand.grid(lambda1 = lambda1, 
                                value = NA, 
                                aic = NA))
  res$id <- 1:nrow(res)
  
  # estimate the graphs for the different values of lambda1 --------
  # go over each pair of penalty terms
  est <- foreach(i = 1:(length(lambda1)), 
                 .options.snow = opts) %dopar% {
                   lapply(1:m, function(k) { 
                      glasso::glasso(s = Sigma[[k]], rho = res$lambda1[i])
                   })
                 }

  
  # Process results -----------------------------------------
  for (i in 1:(length(lambda1))) { 
    for (k in 1:m) { 
      global_res$Theta[[i]][[k]] <- est[[i]][[k]]$wi
      
      # create the adjacency matrix given the precision matrix
      global_res$adj_matrices[[i]][[k]] <- Matrix( as.numeric( abs(est[[i]][[k]]$wi) >= 2*.Machine$double.eps), ncol = ncol(est[[i]][[k]]$wi) , sparse = TRUE)
      diag(global_res$adj_matrices[[i]][[k]]) <- 0
    }
    
    res$aic[i] <- CVN::determine_information_criterion(Theta = global_res$Theta[[i]], 
                                                       adj_matrices = global_res$adj_matrices[[i]], 
                                                       Sigma = Sigma, 
                                                       n_obs = n_obs,
                                                       type = "AIC") 
  }
  
  # stop the progress bar 
  if (verbose) { 
    pb$terminate()
  }
  
  # stop the cluster
  stopCluster(cl) 
  
  # Collect all the results & input ---------------------------
  global_res$results  <- res                  
  
  class(global_res) <- c("cvn", "cvn:glasso", "list")
  
  return(global_res)
}
