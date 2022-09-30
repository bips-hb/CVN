#' @export
genlasso_wrapper <- function(y, W, m, c, eta1, eta2, a, rho, max_iter, eps, truncate) { 
  genlassoRcpp(y, W, m, c, eta1, eta2, a, rho, max_iter, eps, truncate) 
}

#' @export
determine_steps <- function(m) { 
  determine_steps_Rccp(m)
}

#' @export
show_indices <- function(m, W) { 
  show_indices_Rcpp(m, W) 
}