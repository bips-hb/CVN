#' Print contents of CVN object
#' 
#' @title Print CVN
#' @param x  Object of classs 'CVN'
#' @param ... Additional arguments to pass
#' @seealso \code{\link{CVN}}
#' @examples 
#' path <- system.file("cvnfit.RData", package = "CVN")
#' load(path)
#' 
#' print(fit)
#' 
#' @export
print.cvn <- function(x, ...) {  
  cat(sprintf("Covariate-varying Network (CVN)\n\n"))

  if (all(x$results$converged)) {
    cat(green(sprintf("\u2713 all converged\n\n")))
  } else {
    cat(red(sprintf("\u2717 did not converge (maxiter of %d not sufficient)\n\n", x$maxiter)))
  }

  # print following variables
  cat(sprintf("Number of graphs (m)    : %d\n", x$m))
  cat(sprintf("Number of variables (p) : %d\n", x$p))
  cat(sprintf("Number of lambda pairs  : %d\n\n", x$n_lambda_values))

  # If x is an object returned by CVN::glasso() it doesn't have a $W element and this errors
  # Causes example in glasso.R to error during R CMD check
  # TODO: Check if glasso() behaves correctly, possibly adjust print method for CVN:glasso subclass?
  if (!is.null(x$W)) {
    cat(sprintf("Weight matrix (W):\n"))
    print(Matrix(x$W, sparse = T))
  }

  cat(sprintf("\n"))
  print(x$results)
}