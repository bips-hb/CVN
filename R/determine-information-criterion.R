#' @references 
#' Danaher, P., Wang, P., & Witten, D. M. (2014). The joint graphical 
#' lasso for inverse covariance estimation across multiple classes. 
#' Journal of the Royal Statistical Society. Series B: Statistical 
#' Methodology, 76(2), 373–397. https://doi.org/10.1111/rssb.12033
#' @keywords internal
determine_information_criterion <- function (Theta, adj_matrices, Sigma, n_obs, gamma = 0.5) 
{
  if(gamma < 0 | gamma > 1) stop("gamma value must be between 0 and 1")
  
  m <- length(Theta)
  E <- mapply(sum, adj_matrices)              # edges in each subgraph
  
  # Scale Matrices for calculating determinant
  T_scaled <- lapply(Theta, function(x) x / max(abs(x)))
  Sigma_scaled <- lapply(Sigma, function(x) x / max(abs(x)))
  
  # loglik
  loglik <- sapply(1:m, function(i){
    n_obs[i] * (sum(diag(Sigma[[i]] %*% Theta[[i]]))) -
      n_obs[i] * determinant(T_scaled[[i]], logarithm = TRUE)$modulus})
  
  aic <- sum(
    sapply(1:m, function(i) {loglik[i] + 2 * E[i]})      
  )
  
  bic <- sum(
    sapply(1:m, function(i) {loglik[i] + log(n_obs[i]) * E[i]})      
  )
  
  p <- nrow(T_scaled[[1]])
  ebic <- sum(
    sapply(1:m, function(i) {loglik[i] + log(n_obs[i]) * E[i] + 4 * gamma * E[i] * log(p)})    
  )
  data.frame(aic, bic, ebic)
}



#' Information Criteria for a \code{cvn} object
#' 
#' Determines a given information criteria for a \code{cvn} object, 
#' see \code{\link{CVN}}. 
#' 
#' @param cvn A CVN object, see \code{\link{CVN}}
#' @param gamma The gamma value for the eBIC (default: 0.5)
#'
#' @export
determine_information_criterion_cvn <- function (cvn, gamma = 0.5) 
{
  if (!("cvn" %in% class(cvn))) {
    stop("cvn parameter must be of type 'cvn'")
  }
  sapply(1:cvn$n_lambda_values, function(i) {
    determine_information_criterion(cvn$Theta[[i]],
                                    cvn$adj_matrices[[i]], 
                                    cvn$Sigma, 
                                    cvn$n_obs, 
                                    gamma = gamma)
  })
}
