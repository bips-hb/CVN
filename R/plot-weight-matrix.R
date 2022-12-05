#' Plot Weight Matrix
#' 
#' Returns a heat map of a weight matrix
#'
#' @param W Symmetric weight matrix 
#' @param title Title plot (Default is none)
#' @param legend_label Title of the legend (Default: "weight")
#' @param add_counts_to_cells If \code{TRUE}, counts from the matrix are 
#'                  added to the plot (Default: \code{TRUE})
#' @param add_ticks_labels If \code{TRUE}, the number corresponding to the graph
#'                   is add to the plot (Default: \code{TRUE})
#' @param t Distance between tick labels and x-axis (Default: -6)
#' @param r Distance between tick labels and y-axis (Default: -8)
#' 
#' @return A heatmap plot                   
#' @export
plot_weight_matrix <- function(W,
                               title = "",
                               legend_label = "weight",
                               add_counts_to_cells = TRUE,
                               add_ticks_labels = TRUE,
                               t = -6,
                               r = -8) {

  m <- nrow(W)
  
  W <- t(apply(W, 2, rev)) # rotate matrix
  
  colnames(W) <- sapply(1:m, function(i) as.character(i))
  rownames(W) <- sapply(1:m, function(i) as.character(i))
  
  data <- reshape2::melt(W)
  data$Var1 <- m - data$Var1 + 1
  
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
      scale_y_reverse(breaks = 1:m) + 
      theme(axis.text.x = element_text(margin = margin(t = t)), 
            axis.text.y = element_text(margin = margin(r = r))) 
    
  } else { 
    p <- p + theme(axis.text.x = element_blank(), 
                   axis.text.y = element_blank()) +
      scale_y_reverse() 
  }
  
  if (add_counts_to_cells) { 
    p <- p + geom_text(aes(label=value)) 
  }
  
  return(p)
}