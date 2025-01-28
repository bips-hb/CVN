#' @aliases CVN-package NULL
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import ggplot2
#' @import Rcpp
#' @import visNetwork
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#' @importFrom Matrix Matrix
#' @importFrom snow makeCluster
#' @importFrom snow stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom stats cov
#' @importFrom stats optim
#' @importFrom stats runif
#' @importFrom stats var
#' @importFrom utils combn
#' @importFrom dplyr %>% filter
#' @importFrom crayon green red
#' @useDynLib CVN
## usethis namespace: end
NULL


# Silence global variable warning
# Ideally this would be solved more cleanly, but this suffices
utils::globalVariables(c(
  "from", "to", "id",
  "lambda1", "lambda2",
  "gamma1", "gamma2",
  "Theta", "Sigma", "data",
  "Var1", "Var2",
  "value",
  "aic", "bic"
))


#' Covariate-varying Networks
#'
#' Inferring high-dimensional Gaussian graphical networks that 
#' change with multiple discrete covariates. 
#'
#' @name CVN-package
#' @docType package
#' @aliases CVN-package
#' @author Louis Dijkstra\cr Maintainer and contributors:
#' Lukas Burk, Ronja Foraita <foraita@@leibniz-bips.de>
#' @references Dijkstra L, Godt A, Foraita R 
#' \emph{Inferring High-Dimensional Dynamic Networks Changing with Multiple Covariates (2024), Arxiv},
#' \url{https://arxiv.org/abs/2407.19978}.
#' @keywords graphical models 
#' @rdname CVN-package
#' @examples
#'
#' data(grid)
#' W <- create_weight_matrix(type = "grid")
#' 
#' cvn <- CVN(grid, W, lambda1 = 1, lambda2 = 1:2, 
#'            eps = 1e-3, maxiter = 1000, verbose = TRUE)
NULL
