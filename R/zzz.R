#' CVN
#'
#' An estimator for graphical models changing with multiple discrete
#' external covariates
#'
#' @keywords internal
#' @author Louis Dijkstra
#' @import Rcpp
#' @useDynLib CVN
"_PACKAGE"

#' Data for a grid of graphs (3 x 3)
#'
#' Data generated for 9 graphs in total, organized in a grid of 
#' (3x3). See the package \code{CVNSim} for more information 
#' on how the grid is constructed: \url{https://github.com/bips-hb/CVNSim}
#' 
#' @name grid
#' @usage data(grid)
#' @docType data
#' @keywords datasets
#' @format A list of length 9
#' @references \url{https://github.com/bips-hb/CVNSim}
#'
"grid"

utils::globalVariables("aic")
utils::globalVariables("bic")
utils::globalVariables("optim")
utils::globalVariables("lambda1")
utils::globalVariables("lambda2")
