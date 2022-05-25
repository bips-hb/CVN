#' Create matrix D to be used for the Generalized LASSO
#'
#' Generates a matrix D to be used for the generalized LASSO. 
#' We solve a generalized LASSO problem for each egde \code{(s,t)} 
#' for each update step for \code{Z}. 
#' 
#' @param W The \code{m x m}-dimensional uppertriangular weight matrix
#' @param lambda1 The \eqn{\lambda_1} LASSO penalty term 
#' @param lambda2 The \eqn{\lambda_2} global smoothing parameter 
#' @param rho The \eqn{\rho} ADMM's penalty parameter (Default: 1)
#' 
#' @return A \code{(m*(m+1)/2) x m}-matrix 
#'         
#' @references Tibshirani, R. J., & Taylor, J. (2011). 
#'             The solution path of the generalized lasso. 
#'             Annals of Statistics, 39(3), 1335–1371. 
#'             https://doi.org/10.1214/11-AOS878
#' @export
createMatrixD <- function(W, lambda1, lambda2, rho = 1) { 
  
  eta1 <- lambda1 / rho 
  eta2 <- lambda2 / rho 
    
  # number of matrices that need to be estimated
  m <- nrow(W)
    
  # get all unique pairs, i.e., (1,2), (1,3) ... (m - 2, m), (m - 1, m)
  all_unique_pairs <- combn(1:m, 2) 
  
  # number of unique pairs 
  n_unique_pairs <- m*(m-1)/2 
  
  # matrix that will contain all pairs
  C <- matrix(rep(0, n_unique_pairs*m), ncol = m)
  
  # fill in each row of C independently
  sapply(1:n_unique_pairs, function(i) { 
      indices <- all_unique_pairs[,i] 
      combination <- rep(0, m) 
      combination[indices] <- c(1, -1) * W[indices[1], indices[2]]
      C[i, ] <<- combination
    })
  
  # identity matrix 
  I <- diag(m)
  
  # combinating the identity matrix with C
  D <- rbind(eta1 * I, eta2 * C)
  
  return(D)
}