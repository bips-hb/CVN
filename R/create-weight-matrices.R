#' Different Weight Matrices
#' 
#' This function generates different weight matrices for two external covariates, \eqn{W}. 
#' There are several types: 
#' \describe{
#' \item{\code{full} }{All graphs are fully connected with weight 1}
#' \item{\code{glasso} }{All graphs are disconnected with weight 0. This 
#'              mimicks the GLASSO, where each graph is estimated independently}
#' \item{\code{grid} }{A weight matrix for a \eqn{k \times l} grid}
#' \item{\code{uniform-random} }{Fully-connected, but the entries are drawn 
#'                    from a uniform distribution}
#' }
#' 
#' @param type The type of weight matrix
#' @param k Number of categories in the first external covariate
#' @param l Number of categories in the second external covariate
#' @param plot If TRUE, the weight matrix is plotted (default: FALSE)
#' 
#' @importFrom stats runif
#' 
#' @return Weight matrix
#' 
#' @examples
#' W <- create_weight_matrix(type="grid", k=3, l=2)
#' 
#' @export
create_weight_matrix <- function(type = c("full", "glasso", "grid", "uniform-random"), 
                                 k, l, plot = FALSE) { 
  
  if(k %% 1 != 0 | k < 1) stop("k must be a positive integer")
  if(l %% 1 != 0 | l < 1) stop("l must be a positive integer")

  # check correctness input
  if (!type[1] %in% c("full", "glasso", "grid", "uniform-random")) {
    stop("type unknown") 
  }
  
  m <- k * l
  
  switch(type[1], 
         "full"   = {W <- matrix(1, nrow = m, ncol = m)},
         "glasso" = {W <- matrix(0, nrow = m, ncol = m)},
         "grid"   = {W <- create_grid_adjacency(k = k, l = l)
         }, 
         "uniform-random" = {
           W <- matrix(runif(m*m), ncol = m)
           W <- W %*% t(W)
           W <- W / max(W)
           diag(W) <- 0
         }
  ) 
  
  diag(W) <- 0 

  if(plot){
    plot_weight_matrix(W, k, l)
  }

return(W)

}





#' Create grid adjacency
#' 
#' An algorithm which creates a weight matrix in grid structure
#'
#' @param k Number of categories in the first external covariate
#' @param l Number of categories in the second external covariate
#'
#' @return A matrix
#' @keywords internal
#'
#' @examples
#' create_grid_adjacency(k = 2, l = 3)

create_grid_adjacency <- function(k, l) {
  # Total number of nodes
  n <- k * l
  
  # Initialize an n x n adjacency matrix filled with zeros
  adjacency_matrix <- matrix(0, nrow = n, ncol = n)
  
  # Function to calculate the node index at (row, col) in the grid
  node_index <- function(row, col) {
    return((row - 1) * k + col)
  }
  
  # Loop through all nodes in the grid
  for (row in 1:l) {
    for (col in 1:k) {
      current_node <- node_index(row, col)
      
      # Connect to the right neighbor (if within bounds)
      if (col < k) {
        right_node <- node_index(row, col + 1)
        adjacency_matrix[current_node, right_node] <- 1
        adjacency_matrix[right_node, current_node] <- 1  # Symmetric
      }
      
      # Connect to the bottom neighbor (if within bounds)
      if (row < l) {
        bottom_node <- node_index(row + 1, col)
        adjacency_matrix[current_node, bottom_node] <- 1
        adjacency_matrix[bottom_node, current_node] <- 1  # Symmetric
      }
    }
  }
  
  return(adjacency_matrix)
}


