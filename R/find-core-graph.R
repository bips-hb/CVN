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
    stop("cvn argument should be of the type cvn or cvn:glasso")  
  }
  
  # go over all lambda values
  lapply(cvn$adj_matrices, function(adj_matrices) { 
    # sum up all adjacency matrices. If it is one for all m of 
    # them, the sum of the cells sums up to m
      Reduce('+', adj_matrices) == cvn$m
    })
}