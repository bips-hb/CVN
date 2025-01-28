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

  # The real core graph
  lapply(cvn$adj_matrices, function(adj_matrices) {
    # Sum up all adjacency matrices
    summed_matrix <- Reduce('+', adj_matrices)

    # Create a binary matrix by comparing the summed matrix to obj$m
    binary_matrix <- ifelse(summed_matrix == cvn$m, 1, 0)
    # Ensure the result is a matrix
    return(matrix(binary_matrix, nrow = nrow(summed_matrix), ncol = ncol(summed_matrix)))
  })
}




#' Unique Edges
#' 
#' Function that finds edges that are only present in one subgraph of the CVN model.
#'
#' @param cvn  A \code{cvn} object  
#'
#' @return A list of adjacency matrix, one for each 
#'         value of (lambda1, lambda2), showing the edges that are unique for each
#'         single subgraph. 
#' @export
#'
#' @examples
#' data(grid)
#' 
#' W <- create_weight_matrix("grid", 3, 3)
#' 
#' # lambdas:
#' lambda1 = 1  # can also be lambda1 = 1:2 
#' lambda2 = .2
#' 
#' fit <- CVN(grid, W, lambda1 = lambda1, lambda2 = lambda2)
#'                  
#' find_unique_edges(fit)
find_unique_edges <- function(cvn) {
  
  if (!("cvn" %in% class(cvn))) { 
    stop("cvn argument should be of the class 'cvn'")  
  }  

  lapply(cvn$adj_matrices, function(adj_matrices, m = length(cvn$adj_matrices[[1]])) {
    
    unique_edges <- vector("list", m)

    for (i in seq_len(m)) {                                        
      unique_matrix <- Reduce("+", adj_matrices[-i]) == 0 
      unique_edges[[i]] <- adj_matrices[[i]] * unique_matrix # Filter unique edges
    }
    unique_edges
  })

}



#' Edge Overview
#' 
#' Function that gives an overview about the number of edges in each subgraph over all
#' fitted CVN models. 
#' 
#' @param cvn  A \code{cvn} object  
#'
#' @return A data.frame.Shows the number of edges in each subgraph, the number of 
#' edges present in all graphs (core edges) and the number of edges that are unique
#' for a subgraph
#' 
#' @export
#'
#' @examples
#' data(grid)
#' 
#' W <- create_weight_matrix("grid", 3, 3)
#' 
#' # lambdas:
#' lambda1 = 1  # can also be lambda1 = 1:2 
#' lambda2 = .2
#' 
#' fit <- CVN(grid, W, lambda1 = lambda1, lambda2 = lambda2)
#'                  
#' cvn_edge_summary(fit)
#'                  
cvn_edge_summary <- function(cvn){
  
  if (!inherits(cvn, "cvn")) stop("The input object is not of class 'cvn'")

  # How many edges are in each graph
  edges <- sapply(cvn$adj_matrices, function(adj_matrices){mapply(sum, adj_matrices)})

  #edges <- mapply(sum, cvn$adj_matrices[[1]])
  
  # Find the core graph between all single graphs can be queried by
  core_graph <- find_core_graph(cvn)

  # How many edges are in the core graph?
  core_edges <- mapply(sum, core_graph) / 2
  
  # Find unique edges for each subggraph
  unique_edges <- find_unique_edges(cvn)
  
  # How many edges are unique in each graph?
  unique_edges <- sapply(unique_edges, function(x) {mapply(sum, x)})
  
  # Together
  if(length(core_edges) == 1){
     result <- data.frame(edges, core_edges, unique_edges)
  } else {
    result <- data.frame(t(edges), core_edges, t(unique_edges))
    names(result) <- c(paste0("E(g", seq_len(nrow(edges)),")"), "E(core)", 
                       paste0("E(g", seq_len(nrow(edges)),"_u)"))
  }
  result
}


