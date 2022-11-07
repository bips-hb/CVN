#' Structural Hamming Distance 
#' 
#' Returns the structural Hamming distance between multiple
#' graphs
#' 
#' @param adj_matrices A list of adjacency matrices
#' 
#' @return Matrix of Hamming distances
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
  
  return(distance_matrix)
}

#' Structural Hamming Distance for a \code{cvn} Object
#' 
#' Returns the structural Hamming distances
#' 
#' @param cvn A \code{cvn} or \code{cvn:glasso} object 
#'            created by either the \code{\link{CVN::CVN}} or the 
#'            \code{\link{CVN::glasso}} function
#' 
#' @return A list of symmetric matrices. Each matrix contains the structural 
#'         Hamming distances between the different graphs. Each item in the 
#'         list corresponds to one \eqn{(\lambda_1, \lambda_2)} pair
#' @export
hamming_distance <- function(cvn) { 
  
  if (!(class(cvn) %in% c("cvn", "cvn:glasso"))) { 
    stop("input must be a 'cvn' or 'cvn:glasso' object") 
  }
  
  empty_matrix <- matrix(rep(0, cvn$m^2), ncol = cvn$m) 
  distances <- lapply(1:cvn$n_lambda_values, 
                      function(i) {empty_matrix})
  
  # go over all lambda value combinations
  for (k in 1:cvn$n_lambda_values) {
    distances[[k]] <- hamming_distance_adj_matrices(cvn$adj_matrices[[k]])
  }
  
  return(distances)
}