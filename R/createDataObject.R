
createMatrix <- function(nrow = 5, ncol = 7) { 
  matrix(rep(1, nrow * ncol), nrow = nrow) 
}

raw_data <- list(
  createMatrix(nrow = 4), 
  createMatrix(nrow = 8), 
  createMatrix(nrow = 5)
  )

#lambda1, lambda2, rho = 1,


# if (lambda1 <= 0) { 
#   stop("lambda1 has to be > 0") 
# }
# 
# if (lambda2 < 0) { 
#   stop("lambda2 has to be >= 0") 
# }
# 
# if (rho <= 0) { 
#   stop("rho has to be > 0") 
# }

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
#' @param W The \code{m x m}-dimensional uppertriangular weight matrix
#' @param normalized Data is normalized if TRUE. Otherwise the data is just 
#'                   centered (Default: TRUE)
#' @param Verbose Verbose (Default: FALSE) 
createCVNDataObject <- function(raw_data, W, normalized = TRUE, verbose = TRUE) { 
  
  if (verbose) { 
    cat("Creating CVN Data Object...")
  }
  
  ### check the correctness of the datasets 
  
  if (!is.list(raw_data)) { stop("raw_data input should be a list") }
  if (length(raw_data) == 0) { stop("raw_data is an empty list") }
  
  # total number of graphs
  m <- length(raw_data) 
  
  ncols_per_dataset <- sapply(raw_data, function(X) ncol(X))
  if (var(ncols_per_data) != 0) { 
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
  
  if (any(n_obs) < 2) { 
    stop("Not enough observations") 
  }
  
  if (verbose) { 
    cat("\u2713 | Raw datasets valid" )
    cat(sprintf("\tNumber of graphs: %d", m))
    #cat(sprintf("\tNumber of observations: %s", m))
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
  
  
  
  
  
  
  
  list(
      n_obs = n_obs, 
      m = m, 
      p = p, 
      W = W,
    ) 
  
}