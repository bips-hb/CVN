#' Heat Map of the Distances between Graphs 
#' 
#' Returns a heat map of the distance matrix for 
#' a particular CVN 
#'
#' @param distance_matrix Symmetric matrix with distances 
#' @param absolute If \code{FALSE}, rescaled to [0,1]
#' @param limits The limits for the values of the Hamming distance 
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
plot_hamming_distances <- function(distance_matrix, 
                                   absolute = TRUE,
                                   limits = c(NA, NA),  
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
          axis.ticks.y = element_blank()) 
    
  if (is.na(limits[1])) { 
    p <- p + scale_fill_continuous(name = legend_label)  
  } else { 
    p <- p + scale_fill_continuous(name = legend_label, limits = limits)   
  }
  
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

#' Heatmaps for a CVN
#' 
#' Creates all the heatmaps for a CVN, a heatmap for each 
#' pair of \eqn{(\lambda_1, \lambda_2)}
#' 
#' @param cvn A \code{cvn} object
#' @param absolute If \code{FALSE}, rescaled to [0,1]
#' @param same_range If \code{TRUE}, all heatmaps have the same range of values 
#'                   of the Hamming distance shown (Default: TRUE)
#' @param titles Title of the plots (Default is none)
#' @param legend_label Title of the legend (Default: "Hamming Distance")
#' @param add_counts_to_cells If \code{TRUE}, counts from the matrix are 
#'                  added to the plot (Default: \code{TRUE})
#' @param add_ticks_labels If \code{TRUE}, the number corresponding to the graph
#'                   is add to the plot (Default: \code{TRUE})
#' @param t Distance between tick labels and x-axis (Default: -6)
#' @param r Distance between tick labels and y-axis (Default: -8)
#' @param verbose If \code{TRUE}, shows progress bar (Default: \code{TRUE})
#' 
#' @return List of plots 
#' @export
plot_hamming_distances_cvn <- function(cvn,
                                       absolute = TRUE,
                                       same_range = TRUE,
                                       titles = rep("", cvn$n_lambda_values),
                                       legend_label = "Hamming Distance",
                                       add_counts_to_cells = TRUE,
                                       add_ticks_labels = TRUE,
                                       t = -6,
                                       r = -8, 
                                       verbose = TRUE) {
  
  if (!("cvn" %in% class(cvn))) { 
    stop("input must be a 'cvn' object") 
  }  
  
  hamming <- CVN::hamming_distance(cvn, verbose = verbose)
  
  if (verbose) { 
    # progress bar for setting up the edges for the individual graphs
    pb <- progress::progress_bar$new(
      format = "Creating heatmaps [:bar] :percent eta: :eta",
      total = cvn$n_lambda_values + 1, clear = FALSE, width= 80, show_after = 0)
    pb$tick()
  }
  
  if (same_range) { 
    # determine the highest value of Hamming distance observed
    limits <- c(0, max(sapply(hamming$distances, function(M) max(M))))
  } else { 
    limits <- c(NA,NA) 
  }
  
  
  plots <- lapply(1:cvn$n_lambda_values, function(i) { 
    p <- CVN::plot_hamming_distances(hamming$distances[[i]], 
                                absolute = absolute, 
                                limits = limits,
                                title = titles[i],
                                legend_label = legend_label,
                                add_counts_to_cells = add_counts_to_cells,
                                add_ticks_labels = add_ticks_labels,
                                t = t,
                                r = r)
    # update progress bar
    if (verbose) { 
      pb$tick()
    }
    
    return(p)
  })
  
  # close progress bar
  if (verbose) { 
    pb$terminate()
  }
  
  return(
    list(
      m = hamming$m, 
      p = hamming$p, 
      W = hamming$W, 
      distances = hamming$distances,
      results = hamming$results,
      plots = plots 
    )
  )
}