#' @export
aug_genlasso <- function(y, W, m, c, eta1, eta2, a, rho, max_iter, eps) { 
  aug_genlassoRcpp(y, W, m, c, eta1, eta2, a, rho, max_iter, eps) 
}