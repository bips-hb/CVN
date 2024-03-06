#' Update Z using FLSA package
#' 
#' This function updates Z matrix using FLSA package.
#' 
#' @param m Number of graphs
#' @param p Number of covariates
#' @param Theta List of matrices
#' @param Y List of matrices
#' @param connListObj connList object for FLSA package
#' @param connList_is_null Logical value indicating whether all elements of connList are NULL
#' @param eta1 Sparsity parameter
#' @param eta2 Smoothness parameter
#' 
#' @return Updated Z matrix
#' @export
updateZ_flsa_package <- function(m, p, Theta, Y, connListObj, connList_is_null, 
                                 eta1, eta2) {
  # update the diagonal
  Z <- lapply(1:m, function(i) {
    diag(diag(Theta[[i]]) - diag(Y[[i]]))}
  )
  
  for (s in 1:(p-1)) {
    for (t in (s+1):p) {
      # get the y vector 
      y <- sapply(1:m, function(i) Theta[[i]][s,t] + Y[[i]][s,t])
      if(connList_is_null) {
        fit <- as.vector(flsa::flsa(y, lambda1 = eta1, lambda2 = 0))
      } else {
        fit <- as.vector(flsa::flsa(y, lambda1 = eta1, lambda2 = eta2, connListObj = connListObj))
      }
      
      sapply(1:m, function(i) {
        Z[[i]][s,t] <<- fit[i]
        Z[[i]][t,s] <<- fit[i]
      })
    }
  }
  
  # # update the off-diagonal
  # apply(combn(p,2), 2, function(pair) {
  #   s <- pair[1]
  #   t <- pair[2]
  #   
  #   # get the y vector 
  #   y <- sapply(1:m, function(i) Theta[[i]][s,t] + Y[[i]][s,t])
  #   
  #   if(connList_is_null) {
  #     fit <- as.vector(flsa::flsa(y, lambda1 = eta1, lambda2 = 0))
  #   } else {
  #     fit <- as.vector(flsa::flsa(y, lambda1 = eta1, lambda2 = eta2, connListObj = connListObj))
  #   }
  #   
  #   sapply(1:m, function(i) {
  #     Z[[i]][s,t] <<- fit[i]
  #     Z[[i]][t,s] <<- fit[i]
  #   })
  # })
  
  return(Z)
}
