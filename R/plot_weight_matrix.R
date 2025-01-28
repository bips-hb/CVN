#' Plot weight matrix as a grid
#' 
#' This function requires the igraph package to be installed.
#'
#' @param W The weight matrix which equals an adjacency matrix
#' @param k Number of categories in the first external covariate
#' @param l Number of categories in the second external covariate
#'
#' @return An igraph object
#' 
#' @seealso \code{\link{hd_weight_matrix}}
#' 
#' @export
#'
#' @examples
#' # Requires to have the igraph package to be installed
#' W <- create_weight_matrix(type = "grid", k = 2, l = 3, plot = FALSE)
#' plot_weight_matrix(W, k=2, l=3)
#' 
plot_weight_matrix <- function(W, k, l){

  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("The 'igraph' package is required for this function. Please install it.",
         call. = FALSE)
  }
  
  W <- round(W, 1)
  
  g <- igraph::graph_from_adjacency_matrix(W, mode = "undirected", 
                                           weighted = TRUE, diag = FALSE)
  plot(g, layout = igraph::layout_on_grid(g, width = k, height = l), 
          vertex.size = 20, vertex.label.cex = 1,
          edge.width = igraph::E(g)$weight * 5,
          vertex.color = "skyblue") 
}
