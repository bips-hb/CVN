#' Structural Hamming Distance 
#' 
#' Returns the structural Hamming distance between multiple
#' graphs
#' 
#' @param adj_matrices A list of adjacency matrices
#' 
#' @return Matrix of Hamming distances
#' @importFrom Matrix Matrix
#' @export
hamming_distance_adj_matrices <- function(adj_matrices) { 
  
  # number of graphs
  m <- length(adj_matrices)
  
  # initial distance matrices 
  distance_matrix <- matrix(rep(0, m^2), ncol = m) 
  
  # go over all pair of graphs 
  for (i in 1:(m-1)) { 
    for (j in (i+1):m) { 
      dist <- as.matrix(adj_matrices[[i]] != adj_matrices[[j]])
      dist <- sum( dist | t(dist) ) / 2
        
      distance_matrix[i,j] <- dist
      distance_matrix[j,i] <- dist
    }
  }
  
  class(distance_matrix) <- "cvn:distancematrix"
  return(distance_matrix)
}

#' Structural Hamming Distance for a \code{cvn} Object
#' 
#' Returns the structural Hamming distances
#' 
#' @param cvn A \code{cvn} object 
#'            created by the \code{\link[CVN]{CVN}} function
#' @param verbose If \code{TRUE}, shows a progress bar
#' 
#' @return A list of symmetric matrices. Each matrix contains the structural 
#'         Hamming distances between the different graphs. Each item in the 
#'         list corresponds to one \eqn{(\lambda_1, \lambda_2)} pair
#' @seealso \code{\link[CVN]{CVN}}
#' @export
hamming_distance <- function(cvn, verbose = TRUE) { 
  
  if (!("cvn" %in% class(cvn))) { 
    stop("input must be a 'cvn' object") 
  }
  
  if (verbose) { 
    cat(sprintf("Determining Hamming distances between the graphs...\n\n")) 
    # progress bar for setting up the edges for the individual graphs
    pb <- progress::progress_bar$new(
      format = "Computing Hamming distance [:bar] :percent eta: :eta",
      total = cvn$n_lambda_values + 1, clear = FALSE, width= 80, show_after = 0)
    pb$tick()
  }
  
  empty_matrix <- matrix(rep(0, cvn$m^2), ncol = cvn$m) 
  distances <- lapply(1:cvn$n_lambda_values, 
                      function(i) {empty_matrix})
  
  # go over all lambda value combinations
  for (k in 1:cvn$n_lambda_values) {
    distances[[k]] <- hamming_distance_adj_matrices(cvn$adj_matrices[[k]])
    
    if (verbose) { 
      pb$tick() 
    }
  }
  
  results <- list(
    m = cvn$m, 
    p = cvn$p, 
    W = cvn$W,
    distances = distances,
    results = cvn$results
  )
  
  # stop the progress bar
  if (verbose) { 
    pb$terminate() 
  }
  
  class(results) <- c("cvn:distances", "list")
  return(results)
}