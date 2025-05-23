#' Extract CVN
#' 
#' Function that extracts one CVN if the CVN object contains more than one. 
#' Its helpful when just one tuning parameter constellation wants to be, e.g., plotted.
#'
#' @param cvn A CVN object, see \code{\link{CVN}}
#' @param id  Integer; Which id from the CVN object should be extracted
#'
#' @return A CVN object
#' @export
#'
#' @examples
#' \dontrun{
#' # Example code of usage (not run because of longer running time)
#' 
#' data(grid)
#' W <- create_weight_matrix("grid", 3, 3)
#' fit <-  CVN(grid, W, n_cores = 1,
#'             lambda1 = c(0.5, 1), 
#'             lambda2 = c(0.1, 0.5)) 
#' (fit2 <- extract_cvn(fit, id = 2))
#' }
#' 
#' 
extract_cvn <- function(cvn, id){
  
  if (!inherits(cvn, "cvn")) {
    stop("The input object is not of class 'cvn'")
  }
  if(!(id %in% cvn$results$id)){
     stop("id must be a number within cvn$results$id")
  }
  
  # Subset the class objects
  new_cvn <- list(
    Theta = list(cvn$Theta[[id]]),
    adj_matrices = list(cvn$adj_matrices[[id]]),
    Sigma = cvn$Sigma,
    m = cvn$m,
    p = cvn$p,
    n_obs = cvn$n_obs,
    data = cvn$data,
    W = cvn$W,
    maxiter = cvn$maxiter,
    rho = cvn$rho,
    eps = cvn$eps,
    truncate = cvn$truncate,
    rho_genlasso = cvn$rho_genlasso,
    eps_genlasso = cvn$eps_genlasso,
    maxiter_genlasso = cvn$maxiter_genlasso,
    truncate_genlasso = cvn$truncate_genlasso,
    n_lambda_values = 1,
    normalized = cvn$normalized,
    warstart = cvn$warmstart,
    minimal = cvn$minimal,
    results = cvn$results[id,]
    )
  
  # Set the class back to "cvn"
  class(new_cvn) <- c("cvn", "list")
  
  return(new_cvn)
  
}