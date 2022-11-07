#' Heat Map of the Distances between Graphs 
#' 
#' Returns a heat map of the distance matrix for 
#' a particular CVN 
#'
#' @param distance_matrix Symmetric matrix with distances 
#' @param absolute If \code{FALSE}, rescaled to [0,1]
#' @param title Title plot (Default is none)
#' @param legend_label Title of the legend (Default: "Hamming Distance")
#' @param add_counts_to_cells If \code{TRUE}, counts from the matrix are 
#'                  added to the plot (Default: \code{TRUE})
#' @param add_ticks_labels If \code{TRUE}, the number corresponding to the graph
#'                   is add to the plot (Default: \code{TRUE})
#' @param t Distance between tick labels and x-axis (Default: -6)
#' @param r Distance between tick labels and y-axis (Default: -8)
#' 
#' @return A heatmap plot                   
#' @export 
plot_hamming_distances <- function(distance_matrix, absolute = TRUE,
                                   title = "", 
                                   legend_label = "Hamming Distance", 
                                   add_counts_to_cells = TRUE, 
                                   add_ticks_labels = TRUE, 
                                   t = -6, 
                                   r = -8) { 
   
  if (!absolute) { 
    distance_matrix <- distance_matrix / max(distance_matrix) 
    
    if (legend_label == "Hamming Distance") { # user did not choose a specific legend label
      legend_label <- "Relative Hamming Distance" 
    }
  }
  
  m <- nrow(distance_matrix)
  distance_matrix <- t(apply(distance_matrix, 2, rev)) # rotate matrix
  
  colnames(distance_matrix) <- sapply(1:m, function(i) as.character(i))
  rownames(distance_matrix) <- sapply(1:m, function(i) as.character(i))
  
  data <- reshape2::melt(distance_matrix)
  
  p <- ggplot(data = data, aes(x = Var1, y = Var2, fill = value)) + 
    geom_tile() +
    ggtitle(title) + 
    xlab("") + 
    ylab("") + 
    theme(axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank()) + 
    scale_fill_continuous(name = legend_label)  
    
  if (add_ticks_labels) { 
    p <- p + 
      scale_x_continuous(breaks = 1:m) + 
      scale_y_continuous(breaks = 1:m) + 
      theme(axis.text.x = element_text(margin = margin(t = t)), 
            axis.text.y = element_text(margin = margin(r = r)))
      
  } else { 
    p <- p + theme(axis.text.x = element_blank(), 
                   axis.text.y = element_blank())  
  }
  
  if (add_counts_to_cells) { 
    p <- p + geom_text(aes(label=value)) 
  }

  return(p)
}