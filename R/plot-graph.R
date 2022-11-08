


#' @export 
create_nodes_visnetwork <- function(n_nodes, labels = 1:n_nodes) {
  nodes <- data.frame(id = 1:data$p)
  nodes$title <- labels
  return(nodes)
}

#' Create a \code{data.frame} for the Edges for \code{visNetwork} 
#' 
#' In order to visualize a graph, we need to create a 
#' \code{data.frame} that can be used by the \code{visNetwork} package. 
#' This function returns the needed \code{data.frame} given a
#' adjacency matrix. 
#' 
#' @param adj_matrix A symmetric adjacency matrix
#' 
#' @return Data frame that be used as input for \code{visNetwork} 
#' @export
create_edges_visnetwork <- function(adj_matrix) { 
  # set the lower diagonal to zero so that we do not count the same edge twice
  adj_matrix[lower.tri(adj_matrix)] <- 0
  
  # find the indices of the non-zero elements in the adjacency matrix
  edges <- as.data.frame(which(adj_matrix != 0, arr.ind = T)) 
  
  # rename the columns of the data.frame to 'from' and 'to' which 
  # is required for visNetwork
  edges %>% 
    rename(
      from = row,
      to = col
    ) %>% 
    arrange(from)
}

#' Add Attributes to Subset of Edges for \code{visNetwork}
#' 
#' A subset of edges can be assign a different thickness 
#' or color. 
#' 
#' @param edges A data.frame create by \code{\link{create_edges_visnetwork}}
#' @param subset_edges A list with the elements \code{from} and \code{to}. Both 
#'              \code{from} and \code{to} are vectors of the same length 
#'              denoting the different edges 
#' @param width Vector with two values. The first is assigned to the 
#'              edges in the subset given by \code{subset_edges}. The second
#'              value is assigned to the rest. If \code{width = c(NA,NA)}, 
#'              no width is assigned
#' @param color Vector with two values. The first is assigned to the 
#'              edges in the subset given by \code{subset_edges}. The second
#'              value is assigned to the rest. If \code{color = c(NULL,NULL)}, 
#'              no color is assigned
#'              
#' @return A data frame that can be used by the \code{visNetwork} package
#' @export
set_attributes_to_edges_visnetwork <- function(edges, 
                                               subset_edges, 
                                               width = c(NA, NA), 
                                               color = c(NULL, NULL)) { 
  
  # add IDs to the edges data.frame in order to keep track of the 
  # edges. The 'id' column is removed later on
  edges$id <- 1:nrow(edges)
  
  # filter out the edges that are in the subset_edges list
  ids <- edges %>% filter(
    from %in% subset_edges$from, 
    to %in% subset_edges$to
  )
  
  # find the id codes of the edges in the 'in_group', i.e., 
  # edges in the subset_edges. The 'out_group' is the rest
  in_group <- ids$id
  out_group <- edges$id[-in_group]
  
  # Setting the width of the in- and out group differently 
  # if width is given (not NA)
  if (!is.na(width[1])) { 
    edges$width[1:nrow(edges)] <- NA
    edges$width[in_group] <- width[1]
    edges$width[out_group] <- width[2]
  }
  
  # Setting the color of the in- and out group differently 
  # if color is given (not NULL)
  if (!is.null(color[1])) { 
    edges$color[1:nrow(edges)] <- NA
    edges$color[in_group] <- color[1]
    edges$color[out_group] <- color[2]
  }
  
  # remove the initially added 'id' column
  edges %>% select(-id)
}

#' A \code{visNetwork} plot
#' 
#' Creates a \code{visNetwork} plot given a list of 
#' nodes and edges. The nodes data frame can be 
#' created with \code{\link{create_nodes_visnetwork}}; 
#' the edges with \code{create_edges_visnetwork}.  
#' In order to highlight edges, you can use 
#' \code{\link{set_attributes_to_edges_visnetwork}}. 
#' 
#' @return A \code{visNetwork} plot
#' @export
visnetwork <- function(nodes, 
                       edges, 
                       node_titles = 1:nrow(adj_matrix), 
                       title = "", 
                       igraph_layout = "layout_in_circle") { 
  visNetwork(nodes, edges, width = "100%", main = list(text = title)) %>% 
    visIgraphLayout(layout = igraph_layout) %>%
    visOptions(highlightNearest = list(enabled = T, hover = T))
}


nodes <- create_nodes_visnetwork(nrow(A), labels = raw_input$gene_labels)
edges <- create_edges_visnetwork(A)

subset_edges <- list(from = c(1, 114), to = c(42, 128))
edges <- set_attributes_to_edges_visnetwork(edges, subset_edges, width = c(8, 2), color = c("red", "blue"))

visnetwork(nodes, edges)