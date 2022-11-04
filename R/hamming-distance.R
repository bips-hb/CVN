#' Structural Hamming Distance 
#' 
#' Returns the structural Hamming distance between two
#' graphs
#' 
#' @param adj_matrix1,adj_matrix2 Two adjacency matrices
#' 
#' @return Hamming distance
#' @export
hamming_distance_adj_matrices <- function(adj_matrix1, adj_matrix2) { 
  diff <- as.matrix(adj_matrix1 != adj_matrix2)
  return( sum( diff | t(diff) ) / 2 )
}

#' Structural Hamming Distance for a \code{cvn} Object
#' 
#' Returns the structural Hamming distances
#' 
#' @param cvn A \code{cvn} or \code{cvn:glasso} object 
#'            created by either the \code{\link{CVN::CVN}} or the 
#'            \code{\link{CVN::glasso}} function
#' 
#' @return A symmetric matrix with the structural Hamming distances between 
#'         the different graphs 
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
    # go over all pair of graphs 
    for (i in 1:(m-1)) { 
      for (j in (i+1):m) { 
        dist <- hamming_distance_adj_matrices(cvn$adj_matrices[[k]][[i]], 
                                              cvn$adj_matrices[[k]][[j]])
        
        distances[[k]][i,j] <- dist
        distances[[k]][j,i] <- dist
      }
    }
  }
  
  return(distances)
}