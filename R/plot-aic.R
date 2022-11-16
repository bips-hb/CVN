#' Heat Map of Performance (AIC)  
#' 
#' Returns a heat map of the AIC 
#'
#' @param cvn 
#' @param limits The limits for the values of the Hamming distance 
#' @param title Title plot (Default is none)
#' 
#' @return A heatmap plot  
#' @export                  
plot_aic <- function(cvn,
                     limits = c(NA, NA),
                     title = "",
                     xlabel = "lambda1",
                     ylabel = "lambda2",
                     legend_label = "AIC") {
  
  p <- ggplot(data = cvn$results) + 
    geom_tile(aes(x = lambda1, y = lambda2, fill = aic)) + 
    ggtitle(title) + 
    xlab(xlabel) + 
    ylab(ylabel)  
    
  if (is.na(limits[1])) { 
     p <- p + scale_fill_continuous(name = legend_label)  
  } else { 
     p <- p + scale_fill_continuous(name = legend_label, limits = limits)   
  }
  
  return(p)
}
