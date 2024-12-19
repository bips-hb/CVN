#' The Core Graph
#' 
#' Finds the 'core graph', i.e., those edges that are 
#' present in all of the estimated graphs 
#' 
#' @param cvn A \code{cvn} object  
#' 
#' @return A list of adjacency matrix, one for each 
#'         value of (lambda1, lambda2)
#' @export
find_core_graph <- function(cvn) { 
  if (!("cvn" %in% class(cvn))) { 
    stop("cvn argument should be of the type cvn")  
  }

  # Louis code includes all edges that are not shared by any other graph  
  # # go over all lambda values
  lapply(cvn$adj_matrices, function(adj_matrices) {
    # sum up all adjacency matrices. If it is one for all m of
    # them, the sum of the cells sums up to m
      Reduce('+', adj_matrices) == cvn$m
    })
  
  # The real core graph
  # lapply(cvn$adj_matrices, function(adj_matrices) {
  #   # Sum up all adjacency matrices
  #   summed_matrix <- Reduce('+', adj_matrices)
  # 
  #   # Create a binary matrix by comparing the summed matrix to obj$m
  #   binary_matrix <- ifelse(summed_matrix == fit6$m, 1, 0)
  #   # Ensure the result is a matrix
  #   return(matrix(binary_matrix, nrow = nrow(summed_matrix), ncol = ncol(summed_matrix)))
  # })
}

