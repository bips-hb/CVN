#' Plot Function for CVN Object Class
#' 
#' Custom plot method for CVN objects.
#' 
#' @param x \code{cvn} object
#' @param ... Additional arguments to pass to \code{\link{CVN}}
#' @aliases plot.cvn
#' @method plot cvn 
#' @seealso \code{CVN}
#'
#' @export
plot.cvn <- function(x, ...) {
  visnetwork_cvn(x)
}

