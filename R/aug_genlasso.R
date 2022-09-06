#' @export
aug_genlasso <- function(y, W, m, c, lambda1, lambda2, global_rho, a, rho, max_iter, eps) { 
  aug_genlassoRcpp(y, W, m, c, lambda1, lambda2, global_rho, a, rho, max_iter, eps) 
}