#' Wrapper for update-Z 
#' 
#' Wrapper for the updateZRcpp function
#' See \code{\link{updateZ}} for more 
#' information
#' 
#' @export
updateZ_wrapper <- function(m, p, nrow_D, 
                            Theta, Y, W, eta1, eta2, a, 
                            rho_genlasso, maxiter_genlasso, eps_genlasso, 
                            truncate_genlasso) { 
  
  updateZRcpp(m, p, nrow_D, 
               Theta, Y, W, eta1, eta2, a, 
               rho_genlasso, maxiter_genlasso, eps_genlasso, 
               truncate_genlasso)  
  
}
  