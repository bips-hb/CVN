#' Heat Map Plot of the Weight Matrix
#' 
#' Returns a heat map of the weight matrix
#'
#' @param W Weight matrix 
#' @param title Title plot (Default is none)
#' @param legend_label Title of the legend (Default: "weight")
#' @param add_counts_to_cells If \code{TRUE}, counts from the matrix are 
#'                  added to the plot (Default: \code{TRUE})
#' @param add_ticks_labels If \code{TRUE}, the number corresponding to the graph
#'                   is add to the plot (Default: \code{TRUE})
#' @param t Distance between tick labels and x-axis (Default: -6)
#' @param r Distance between tick labels and y-axis (Default: -8)
#' 
#' @importFrom Matrix Matrix
#' @import ggplot2 
#' 
#' @seealso \code{\link{plot_weight_matrix}}
#'
#' @return A heatmap plot       
#' 
#' @examples
#' W_uniform <- create_weight_matrix(type="uniform-random", 3,3)
#' W_uniform <- round(W_uniform, 2)
#' hd_weight_matrix(W_uniform)
#'             
#' @export
hd_weight_matrix <- function(W,
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
  
  data <- melt(W)
  data$Var1 <- m - data$Var1 + 1
  
  # needed for package building: visible binding for global variables Var1,...
  Var1 <- Var2 <- value <- geom_tile <- NULL
  
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