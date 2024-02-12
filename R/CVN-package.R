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
