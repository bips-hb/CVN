#' Plot Function for CVN Object Class
#' 
#' Custom plot method for CVN objects.
#' 
#' @param x \code{cvn} object
#' @param ... Additional arguments to pass to \code{\link{visnetwork_cvn}}
#' @aliases plot.cvn visnetwork_cvn
#' @method plot cvn 
#' @seealso \code{CVN}
#'
#' @export
plot.cvn <- function(x, ...) {
  if (!inherits(x, "cvn")){
    stop("Object must be of class 'cvn'.")
  }
  visnetwork_cvn(x, ...)
}


