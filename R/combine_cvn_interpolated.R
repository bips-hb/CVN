#' Combine CVN object with interpolated CVN object
#' 
#' Combine an interpolated CVN to an original CVN
#'
#' @param cvn A CVN fit with \eqn{m} graphs
#' @param cvn_interpolated A interpolated CVN from a fitted CVN model with \eqn{m} graphs
#'
#' @return A 'cvn_interpolated' object (see \code{\link{interpolate}}).
#' @examples
#' path <- system.file("cvnfit.rda", package = "CVN")
#' load(path)
#' 
#' interpolate <- interpolate(fit, c(0,0,0,0,0,0,0,0.5,0.5), truncate = 0.05)
#' fit10 <- combine_cvn_interpolated(fit, interpolate)
#' 
#' @export
#'
combine_cvn_interpolated <- function(cvn, cvn_interpolated){
  
  if (!inherits(cvn, "cvn")) {
    stop("The input object is not of class 'cvn'.")
  }
  if (!inherits(cvn_interpolated, "cvn_interpolated")) {
    stop("The input object is not of class 'cvn_interpolated'.")
  }
  if (cvn$m != cvn_interpolated$m) {
    stop("m must be the same in both objects.")
  }
  if (cvn$p != cvn_interpolated$p) {
    stop("p must be the same in both objects.")
  }
  if (!all(cvn$W == cvn_interpolated$W)) {
    stop("W must be the same in both objects.")
  }
  if (length(cvn$adj_matrices) > 1) {
    stop("Extract one CVN from cvn before merging.")
  }
  if (length(cvn_interpolated$adj_matrices) > 1) {
    stop("Extract one CVN from cvn_interpolated before merging.")
  }
  if (!all(cvn$results[, c(2,3)] == cvn_interpolated$results)) {
    stop("Lambda values in cvn and cvn_interpolated are different; should not be merged.")
  }

  # combine
  cvn_interpolated_combined <- list(
    adj_matrices = list(c(cvn$adj_matrices[[1]], cvn_interpolated$adj_matrices)),
    m = cvn$m + 1, 
    p = cvn$p, 
    W = cvn$W,
    weights = cvn_interpolated$weights, 
    truncate = c(cvn$truncate, cvn_interpolated$truncate), 
    n_lambda_values = cvn$n_lambda_values, 
    results = cvn$results
  )
  
  # Set the class back to "irgendwas"
  class(cvn_interpolated_combined) <- c("cvn_interpolated", "list")
  
  return(cvn_interpolated_combined)
  
}