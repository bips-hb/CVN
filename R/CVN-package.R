#' Covariate-varying Networks
#'
#' Inferring high-dimensional Gaussian graphical networks that 
#' change with multiple discrete covariates. 
#'
#' @name CVN-package
#' @references Dijkstra L, Godt A, Foraita R 
#' \emph{Inferring High-Dimensional Dynamic Networks Changing with Multiple Covariates (2024), Arxiv},
#' \url{https://arxiv.org/abs/2407.19978}.
#' @keywords graphical models 
#' @examples
#'
#' data(grid)
#' W <- create_weight_matrix(type = "grid", k=3, l=3, plot = FALSE)
#' 
#' cvn <- CVN(grid, W, lambda1 = 1, lambda2 = 1:2, 
#'            n_cores = 1,
#'            eps = 1e-2, maxiter = 1000, verbose = TRUE)
"_PACKAGE"


utils::globalVariables(c(
  "from", "to", "id",
  "lambda1", "lambda2",
  "gamma1", "gamma2",
  "Theta", "Sigma", "data",
  "Var1", "Var2",
  "value",
  "aic", "bic", "ebic"
))

## usethis namespace: start
#' @import ggplot2
#' @import Rcpp
## @import visNetwork
#' @importFrom crayon green red 
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#' @importFrom doParallel registerDoParallel
#' @importFrom dplyr %>%  
#' @importFrom dplyr arrange
#' @importFrom dplyr filter  
#' @importFrom dplyr join_by
#' @importFrom dplyr rename
#' @importFrom dplyr select 
#' @importFrom dplyr semi_join
#' @importFrom dplyr slice  
#' @importFrom glasso glasso
#' @importFrom Matrix Matrix 
#' @importFrom methods is
#' @importFrom parallel clusterEvalQ
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom reshape2 melt
#' @importFrom utils combn 
#' @importFrom stats cov
#' @importFrom stats IQR
#' @importFrom stats median
#' @importFrom stats optim
#' @importFrom stats runif
#' @importFrom stats var
#' @importFrom visNetwork visNetwork 
#' @importFrom visNetwork visIgraphLayout 
#' @importFrom visNetwork visOptions
#' @useDynLib CVN
## usethis namespace: end
NULL
