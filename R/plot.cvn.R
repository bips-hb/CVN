#' Plot Function for CVN Object Class
#' 
#' Custom plot method for CVN objects.
#' 
#' @param x \code{cvn} object
#' @param ... Additional arguments to pass to \code{\link{visnetwork_cvn}}
#' @aliases plot.cvn 
#' @method plot cvn 
#' @seealso \code{CVN}
#' 
#' @examples
#' \donttest{
#' path <- system.file("cvnfit.rda", package = "CVN")
#' load(path)
#' fit <- plot(fit) 
#' fit$plots[[1]][[4]]
#' }
#' @export
plot.cvn <- function(x, ...) {
  if (!inherits(x, "cvn")){
    stop("Object must be of class 'cvn'.")
  }
  visnetwork_cvn(x, ...)
}


