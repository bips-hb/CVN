#' Extract CVN
#' 
#' Function that extracts one CVN if the CVN object contains more than one. 
#' Its helpful when just one tuning parameter constellation wants to be, e.g., plotted.
#'
#' @param cvn A CVN object
#' @param nr  Integer; Which id from the CVN object should be extracted
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
#' fit <-  CVN(grid, W, 
#'             lambda1 = c(0.5, 1), 
#'             lambda2 = c(0.1, 0.5)) 
#' (fit2 <- extract_cvn(fit, 2))
#' }
#' 
#' 
extract_cvn <- function(cvn, nr){
  
  if (!inherits(cvn, "cvn")) {
    stop("The input object is not of class 'cvn'")
  }
  if(!(nr %in% cvn$results$id)){
     stop("nr must be a number within cvn$results$id")
  }
  
  if ('Theta' %in% names(cvn)) {
    cvn <- within(cvn, rm(Theta))
  }
  
  # Subset the class objects
  new_cvn <- list(
    Theta = list(cvn$Theta[[nr]]),
    adj_matrices = list(cvn$adj_matrices[[nr]]),
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
    results = cvn$results[nr,]
    )
  
  # Set the class back to "irgendwas"
  class(new_cvn) <- "cvn"
  
  return(new_cvn)
  
}