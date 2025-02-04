#' The \eqn{Y}-update step for the ADMM
#' 
#' Returns the updated value of \eqn{Y} for the ADMM given the previously 
#' updated values of \eqn{\Theta} and \eqn{Z}
#' 
#' @param Theta A list with matrices with the current values of \eqn{\Theta}
#' @param Z A list with matrices with the current values of \eqn{Z}
#' @param Y A list with matrices with the current values of \eqn{Y} 
#' 
#' @return A list with matrices with the new values of \eqn{Y}
#' @noRd        
#' @keywords internal
updateY <- function(Theta, Z, Y) { 
   # Y_new = Y_old + Theta_new - Z_new
  mapply(function(theta, z, y) {y + (theta - z)}, 
         Theta, Z, Y, SIMPLIFY = FALSE)
}
