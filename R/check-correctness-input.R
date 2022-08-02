#' Check whether Input is Valid
#' 
#' Checks whether the input for the function \code{\link{CVN}} is valid.
#' This function does not return anything. The execution of the function
#' halts when an issue has been detected.
#' 
#' @param raw_data A list with matrices. The number of columns should be the 
#'                 same for each matrix 
#' @param W The \eqn{(m \times m)}-dimensional upper-triangular 
#'          weight matrix \eqn{W}
#' @param lambda1 A vector of \eqn{\lambda_1}'s LASSO penalty terms 
#' @param lambda2 A vector of \eqn{\lambda_2}'s global smoothing parameters
#' @param rho The \eqn{\rho} ADMM's penalty parameter 
#'        
#' @seealso \code{\link{CVN}}
#' @export
check_correctness_input <- function(raw_data, W, lambda1, lambda2, rho) { 
  
  if (!is.list(raw_data)) { 
    stop("raw_data input should be a list") 
  }
  
  if (length(raw_data) == 0) { 
    stop("raw_data is an empty list") 
  }
  
  # total number of graphs
  m <- length(raw_data) 
  
  ncols_per_dataset <- sapply(raw_data, function(X) ncol(X))
  
  if (m > 1 && var(ncols_per_dataset) != 0) { 
    stop("The number of columns (variables) should 
              be the same for all datasets in raw_data") 
  }
  
  # number of variables 
  p <- ncols_per_dataset[1]
  
  if (p < 2) { 
    stop("Just one variable in the dataset. At least two are required") 
  } 
  
  # number of observations for each dataset
  n_obs <- sapply(raw_data, function(X) nrow(X))
  
  if (any(n_obs < 2)) { 
    stop("Not enough observations") 
  }
  
  if (!is.matrix(W)) { 
    stop("W must be a matrix") 
  }
  
  if (nrow(W) != m) { 
    stop("Number of rows in W should be equal to m, i.e., the number of graphs")  
  }
  
  if (ncol(W) != m) { 
    stop("Number of columns in W should be equal to m, i.e., the number of graphs")  
  }
  
  if (any(W > 1) || any(W < 0)) { 
    stop("The values in the weight matrix must lie between in the interval [0,1]") 
  }
  
  if (length(lambda1) >= 1) { 
    stop("lambda1 cannot be an empty vector") 
  }
  
  if (length(lambda2) >= 1) { 
    stop("lambda2 cannot be an empty vector") 
  }
  
  if (any(lambda1 < 0)) { 
    stop("Values for lambda1 should be positive") 
  }
  
  if (any(lambda2 < 0)) { 
    stop("Values for lambda2 should be positive") 
  }
  
  if (rho < 0) { 
    stop("rho should be strictly positive") 
  }
  
}



