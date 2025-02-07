#' Nodes for the \code{visNetwork} package
#'
#' Creates a data frame that can be used for the
#' \code{visNetwork} package.
#'
#' @param n_nodes Number of nodes in the graph
#' @param labels The labels for the individual nodes
#'             (Default: \code{1:n_nodes})
#'
#' @return Data frame with two columns: \code{id} and \code{title}  
#' @examples 
#' nodes <- create_nodes_visnetwork(n_nodes = 5, labels = LETTERS[1:5])
#'
#' adj_matrix <- matrix(c(0, 1, 0, 1, 0,
#'                        1, 0, 1, 0, 0,
#'                        0, 1, 0, 0, 0,
#'                        1, 0, 0, 0, 1,
#'                        0, 0, 0, 1, 0), ncol = 5)
#'
#' edges <- create_edges_visnetwork(adj_matrix)
#' 
#' shared_edges <- data.frame(from = c(1,2), to = c(4, 3))
#'
#' edges <- set_attributes_to_edges_visnetwork(edges,
#'                                             subset_edges = shared_edges,
#'                                             width = c(3, .5),
#'                                             color = c("red", "blue"))
#'
#' visnetwork(nodes, edges) 
#' @export

create_nodes_visnetwork <- function(n_nodes, labels = 1:n_nodes) {
  nodes <- data.frame(id = 1:n_nodes)
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
#' @examples 
#' nodes <- create_nodes_visnetwork(n_nodes = 5, labels = LETTERS[1:5])
#'
#' adj_matrix <- matrix(c(0, 1, 0, 1, 0,
#'                        1, 0, 1, 0, 0,
#'                        0, 1, 0, 0, 0,
#'                        1, 0, 0, 0, 1,
#'                        0, 0, 0, 1, 0), ncol = 5)
#'
#' edges <- create_edges_visnetwork(adj_matrix)
#' 
#' shared_edges <- data.frame(from = c(1,2), to = c(4, 3))
#'
#' edges <- set_attributes_to_edges_visnetwork(edges,
#'                                             subset_edges = shared_edges,
#'                                             width = c(3, .5),
#'                                             color = c("red", "blue"))
#'
#' visnetwork(nodes, edges) 
#' @export

create_edges_visnetwork <- function(adj_matrix) {

  # needs to be of the type 'matrix' when using the function which
  adj_matrix <- as.matrix(adj_matrix)

  # set the lower diagonal to zero so that we do not count the same edge twice
  adj_matrix[lower.tri(adj_matrix)] <- 0

  # find the indices of the non-zero elements in the adjacency matrix
  edges <- as.data.frame(which(adj_matrix != 0, arr.ind = T))

  # rename the columns of the data.frame to 'from' and 'to' which
  # is required for visNetwork
  from <- NULL
  rename(edges, from = row, to = col) %>% 
    arrange(from)
}

#' Add Attributes to Subset of Edges for \code{visNetwork}
#' 
#' A subset of edges can be assign a different thickness 
#' or color. 
#' 
#' @param edges A data.frame create by \code{\link[CVN]{create_edges_visnetwork}}
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
#' @examples 
#' nodes <- create_nodes_visnetwork(n_nodes = 5, labels = LETTERS[1:5])
#'
#' adj_matrix <- matrix(c(0, 1, 0, 1, 0,
#'                        1, 0, 1, 0, 0,
#'                        0, 1, 0, 0, 0,
#'                        1, 0, 0, 0, 1,
#'                        0, 0, 0, 1, 0), ncol = 5)
#'
#' edges <- create_edges_visnetwork(adj_matrix)
#' 
#' shared_edges <- data.frame(from = c(1,2), to = c(4, 3))
#'
#' edges <- set_attributes_to_edges_visnetwork(edges,
#'                                             subset_edges = shared_edges,
#'                                             width = c(3, .5),
#'                                             color = c("red", "blue"))
#'
#' visnetwork(nodes, edges)
#' @export
set_attributes_to_edges_visnetwork <- function(edges,
                                               subset_edges,
                                               width = c(NA, NA),
                                               color = c(NULL, NULL)) {

  # check whether the subset list is not empty
  if (length(subset_edges$from) == 0) {
    if (!is.na(width[1])) {
      edges$width <- width[2]
    }

    if (!is.null(color[1])) {
      edges$color <- color[2]
    }
    return(edges)
  }

  # add IDs to the edges data.frame in order to keep track of the
  # edges. The 'id' column is removed later on
  edges$id <- 1:nrow(edges)

  # filter out the edges that are in the subset_edges list
  from <- to <- NULL
  ids <- semi_join(edges, subset_edges, by = join_by(from, to))

  # find the id codes of the edges in the 'in_group', i.e.,
  # edges in the subset_edges. The 'out_group' is the rest
  in_group <- ids$id
  out_group <- edges$id[-in_group]

  # Setting the width of the in- and out group differently
  # if width is given (not  NA)
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
  # id <- NULL
  # edges %>% select(-c(id))
  select(edges, -"id")
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
#' @param nodes A data frame with rows equal to the number of nodes with columns "id" and "title". 
#'              Can be created with \code{\link{create_nodes_visnetwork}}.
#' @param edges A data frame with columns "from" and "to. Each row represents an edge between two nodes (integer)
#' @param node_titles Vector with title of the nodes (Default: \code{1:p}) 
#' @param title A list with \code{n_lambda_values} vectors. Each vector is of the 
#'         lenght \code{m}. Regulates the titles of the graphs (Default: no title)
#' @param igraph_layout igraph layout (default: layout_in_circle)
#' 
#' @return A \code{visNetwork} plot
#' @seealso \code{CVN}, \code{\link{create_nodes_visnetwork}}
#' @examples 
#' nodes <- create_nodes_visnetwork(n_nodes = 5, labels = LETTERS[1:5])
#'
#' adj_matrix <- matrix(c(0, 1, 0, 1, 0,
#'                        1, 0, 1, 0, 0,
#'                        0, 1, 0, 0, 0,
#'                        1, 0, 0, 0, 1,
#'                        0, 0, 0, 1, 0), ncol = 5)
#'
#' edges <- create_edges_visnetwork(adj_matrix)
#' 
#' shared_edges <- data.frame(from = c(1,2), to = c(4, 3))
#'
#' edges <- set_attributes_to_edges_visnetwork(edges,
#'                                             subset_edges = shared_edges,
#'                                             width = c(3, .5),
#'                                             color = c("red", "blue"))
#'
#' visnetwork(nodes, edges)
#' @export
visnetwork <- function(nodes,
                       edges,
                       node_titles = nodes$id,
                       title = "",
                       igraph_layout = "layout_in_circle") {

  if (!(requireNamespace("igraph", quietly = TRUE))) {
    stop("The igraph package is required for visNetwork, please install it first.")
  }

  visNetwork(nodes, edges, width = "100%",
             main = list(text = title)) %>%
    visIgraphLayout(layout = igraph_layout) %>%
    visOptions(highlightNearest = list(enabled = T, hover = T))
}


#' All \code{visNetwork} plots for a CVN object
#'
#' Creates all \code{visNetwork} plots, see \code{\link{visnetwork}},
#' for all graphs in a \code{cvn} object
#'
#' @param cvn A \code{cvn} object, see \code{\link{CVN}}
#' @param node_titles Vector with title of the nodes (Default: \code{1:p})
#' @param titles A list with \code{n_lambda_values} vectors. Each vector is of the
#'         lenght \code{m}. Regulates the titles of the graphs (Default: no title)
#' @param show_core_graph Shall the core graph be visualized (Default = TRUE)
#' @param width Edge width of the core graph
#' @param color String vector with two colors. The first color marks the edges in the
#'              core graph (Default: c("red", "blue")) 
#' @param igraph_layout igraph layout (default: layout_in_circle)
#' @param verbose Verbose (Default: \code{TRUE}) 
#' 
#' @seealso \code{\link{CVN}}, \code{\link{visnetwork}}
#' 
#' @return The cvn input which is extended by the list element 'plots'
#' @export
#' @examples 
#' path <- system.file("cvnfit.RData", package = "CVN")
#' load(path)
#' fit_plot <- visnetwork_cvn(fit)
#' fit_plot$plots[[1]][[1]]
visnetwork_cvn <- function(cvn, 
                           node_titles = 1:cvn$p, 
                           titles = lapply(1:cvn$n_lambda_values, 
                                           function(i) sapply(1:cvn$m, function(j) "")), 
                           show_core_graph = TRUE, 
                           width = c(3,1), 
                           color = c("red", "blue"), 
                           igraph_layout = "layout_in_circle", 
                           verbose = TRUE) {

  if(!inherits(cvn, "cvn") & !inherits(cvn, "cvn_interpolated")) {
    stop("The input object must be of class 'cvn' or 'cvn_interpolated'.")
  }
  if (inherits(cvn, "cvn_interpolated")) {
    warning("Works only if the original cvn was combined with an interpolated cvn using 
            the function 'combine_cvn_interpolated'.")
  }
  if (!(length(node_titles) == cvn$p)) {
    stop("The number of node labels does not correspond to the number of nodes")
  }
  
  if (verbose) {
    cat(sprintf("Creating visNetwork plots for the CVN...\n\n"))
    cat(sprintf("Number of graphs:                  %d\n", cvn$m))
    cat(sprintf("Number of different lambda values: %d\n", cvn$n_lambda_values))
    cat(sprintf("Creating nodes...\n"))
  }
  
  res <- list(
    m = cvn$m,
    p = cvn$p,
    W = cvn$W,
    results = cvn$results
  )
  
  nodes <- create_nodes_visnetwork(n_nodes = cvn$p, labels = node_titles)

  # the edges that are constant in the different graphs are
  # displayed differently
  if (show_core_graph) {
    
    # get the core graphs for the different values of (lambda1, lamdba2)
    if (verbose) {
      cat(sprintf("Determining the 'core graphs'...\n"))
    }
    
    core_graphs <- find_core_graph(cvn)
    
    if (verbose) {
      cat(sprintf("Create the subset of edges in the core graphs...\n\n"))
    }
    
    subset_edges <- lapply(core_graphs, function(adj_matrix) {
      create_edges_visnetwork(adj_matrix)
    })
  }
  
  # Set-up a progress bars ---------------------------------
  if (verbose) {
    # progress bar for setting up the edges for the individual graphs
    pb_edges <- progress::progress_bar$new(
      format = "Creating edge lists [:bar] :percent eta: :eta",
      total = cvn$m * cvn$n_lambda_values + 1, clear = FALSE, width= 80, show_after = 0)
    pb_edges$tick()
    
    pb_plots <- progress::progress_bar$new(
      format = "Creating plots [:bar] :percent eta: :eta",
      total = cvn$m * cvn$n_lambda_values + 1, clear = FALSE, width= 80, show_after = 0)
    pb_plots$tick()
  }
  
  # create the edge dataframes for all the graphs
  # browser()
  all_edges <- lapply(1:cvn$n_lambda_values, function(i) {
    lapply(1:cvn$m, function(k) {
      # cat(sprintf("%d\t%d\n", i,k))
      edges <- create_edges_visnetwork(cvn$adj_matrices[[i]][[k]])
      # check whether there are core edges, since sometimes graphs are
      # completely empty
      if (show_core_graph && length(subset_edges[[i]]$from) != 0) {
        edges <- set_attributes_to_edges_visnetwork(edges, 
                                                    subset_edges = subset_edges[[i]],
                                                    width = width, 
                                                    color = color)
      }
      
      if (verbose) {
        pb_edges$tick()
      }
      
      return(edges)
    })
  })
  
  if (verbose) {
    pb_edges$terminate()
    cat(sprintf("\nCreate plots given the determined edges...\n\n"))
  }
  
  res$plots <- lapply(1:cvn$n_lambda_values, function(i) {
    lapply(1:cvn$m, function(k) {
      if (verbose) {
        pb_plots$tick()
      }
      return(visnetwork(nodes, all_edges[[i]][[k]],
                        title = titles[[i]][[k]],
                        igraph_layout = igraph_layout))
    })
  })
  
  # stop the progress bar
  if (verbose) {
    pb_plots$terminate()
  }
  
  return(res)

}
