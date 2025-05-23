#' Heat Map of an Information Criterion (AIC or BIC)
#' 
#' Returns a heat map of the AIC or BIC for a fitted CVN
#'
#' @param cvn Fitted CVN, see \code{\link{CVN}}
#' @param criterion The information criterion, must be either \code{'aic'} 
#'                  or \code{'bic'}. Default: \code{'bic'}
#' @param use_gammas If \code{TRUE}, plots the \eqn{\gamma}-values. Otherwise, 
#'                   the \eqn{\lambda}-values are used
#' @param show_minimum If \code{TRUE}, an orange dot is put on the point with the 
#'                     minimum value of the information criterion is. If \code{FALSE}, 
#'                     no dot is added. Default: \code{TRUE}. 
#' @param title Title plot (Default is none)
#' @param xlabel Label for the \eqn{x}-axis. Default depends on \code{use_gammas}. 
#'               If \code{use_gammas = TRUE}, then the label is 'gamma1'. Otherwise, 
#'               'lambda1'
#' @param ylabel Label for the \eqn{x}-axis. Default depends on \code{use_gammas}. 
#'               If \code{use_gammas = TRUE}, then the label is 'gamma1'. Otherwise, 
#'               'lambda1'
#' @param legend_label Title for the legend. Default depends on \code{criterion}. 
#'                     If \code{'aic'}, then the label is 'AIC'. 
#'                     Otherwise, 'BIC'. 
#' @param limits The limits for the values of the Hamming distance 
#' 
#' @return A heatmap plot  
#' 
#' @import ggplot2
#' @importFrom dplyr %>% 
#' 
#' @examples 
#' \donttest{
#' data(grid)
#' m <- 9
#' 
#' W <- create_weight_matrix(type="grid", 3, 3)
#' cvn <- CVN(grid, W, 
#'            lambda1 = 1:2, lambda2 = .5, n_cores = 1)
#' cvn$results          
#' 
#' # The smaller the IC value, the better the fit. 
#' plot_information_criterion(cvn, criterion = "aic")
#' plot_information_criterion(cvn, criterion = "bic")
#' plot_information_criterion(cvn, criterion = "ebic")
#' }
#' 
#' @export                  
plot_information_criterion <- function(cvn,
                                       criterion = c('bic', 'aic', 'ebic'),
                                       use_gammas = TRUE,
                                       show_minimum = TRUE,
                                       title = "",
                                       xlabel = NULL,
                                       ylabel = NULL,
                                       legend_label = NULL,
                                       limits = c(NA, NA)) {

  # needed for package building: visible binding for global variables Var1,...
  gamma1 <- gamma2 <- geom_tile <- geom_point <- NULL
  
  # get the criterion. Tolower, since we also accept 'BIC' and 'AIC'
  criterion <- tolower(criterion[1])
  criterion <- match.arg(criterion)
  
  if (!(criterion %in% c('aic', 'bic', 'ebic'))) {
    stop("Selected criterion must be either 'aic', 'bic or ebic'")
  }
  
  # determine the x-, y- and legend labels
  if (is.null(xlabel)) {
    xlabel <- ifelse(use_gammas, "gamma1", "lambda1")
  }
  
  if (is.null(ylabel)) {
    ylabel <- ifelse(use_gammas, "gamma2", "lambda2")
  }
  
  if (is.null(legend_label)) {
    if (criterion == 'aic') {
      legend_label <- "AIC"
    } else if (criterion == 'bic') {
      legend_label <- "BIC"
    } else {
      legend_label <- "eBIC"
    }
  }
  
  # initialize the plot
  p <- ggplot(data = cvn$results)
  
  # add the aes 
  if (criterion == "aic") {
    if (use_gammas) {
      p <- p + geom_tile(aes(x = gamma1, y = gamma2, fill = aic))
    } else {
      p <- p + geom_tile(aes(x = lambda1, y = lambda2, fill = aic))
    }
  }
  
  if (criterion == "bic") {
    if (use_gammas) {
      p <- p + geom_tile(aes(x = gamma1, y = gamma2, fill = bic))
    } else {
      p <- p + geom_tile(aes(x = lambda1, y = lambda2, fill = bic))
    }
  }
  
  if (criterion == "ebic") {
    if (use_gammas) {
      p <- p + geom_tile(aes(x = gamma1, y = gamma2, fill = ebic))
    } else {
      p <- p + geom_tile(aes(x = lambda1, y = lambda2, fill = ebic))
    }
  }
  
  # add the labels and title
  p <- p + ggtitle(title) + xlab(xlabel) + ylab(ylabel)  
  
  # add a dot where the minimum value is 
  if (show_minimum) {
    
    # determine the minimum for either the AIC, BIC or eBIC
    if (criterion == 'aic') {
      minimum <- cvn$results %>% slice(which.min(aic))
    } else if(criterion == 'bic'){
      minimum <- cvn$results %>% slice(which.min(bic))
    } else { 
      minimum <- cvn$results %>% slice(which.min(ebic))
    }
    
    # use either the gammas or lambdas
    if (use_gammas) {
      p <- p + geom_point(minimum, color = "orange", size = 4, 
                          mapping = aes(x = gamma1, y = gamma2, fill = NULL))
    } else {
      p <- p + geom_point(minimum, color = "orange", size = 4, 
                          mapping = aes(x = lambda1, y = lambda2, fill = NULL))
    }
  }
 
  # add the legend and add the limits if there were set  
  if (is.na(limits[1])) { 
     p <- p + scale_fill_continuous(name = legend_label)  
  } else { 
     p <- p + scale_fill_continuous(name = legend_label, limits = limits)   
  }
  
  return(p)
}
