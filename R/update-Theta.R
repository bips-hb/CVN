#' The \eqn{\Theta}-update step for the ADMM
#' 
#' Returns the updated value of \eqn{\Theta} for the ADMM given the previously 
#' updated values of \eqn{Z} and \eqn{Y}
#' 
#' @param Theta A list with matrices with the current values of \eqn{\Theta}
#' @param Z A list with matrices with the current values of \eqn{Z}
#' @param Y A list with matrices with the current values of \eqn{Y} 
#' 
#' @return A list with matrices with the new values of \eqn{\Theta}
#'        
#' @keywords internal
updateTheta <- function(Theta, Z, Y) { 
  # Y_new = Y_old + Theta_new - Z_new
  #mapply(function(y, theta, z) {y + theta - z}, 
  #       Y_old, Theta_new, Z_new, SIMPLIFY = FALSE)
}