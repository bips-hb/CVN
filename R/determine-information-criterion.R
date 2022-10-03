#' @references 
#' Danaher, P., Wang, P., & Witten, D. M. (2014). The joint graphical 
#' lasso for inverse covariance estimation across multiple classes. 
#' Journal of the Royal Statistical Society. Series B: Statistical 
#' Methodology, 76(2), 373–397. https://doi.org/10.1111/rssb.12033
#' @export
determine_information_criterion <- function(Theta, adj_matrices, Sigma, n_obs,
                                            type = c("AIC", "BIC")) { 
  
  # number of graphs
  m <- length(Theta)
  
  # number of non-zero edges
  E <- mapply(sum, adj_matrices) 
  
  if (type[1] == "AIC" || type[1] == "aic") {
    return(sum(
      sapply(1:m, function(i) { 
        n_obs[i] * ( sum(diag( Sigma[[i]] %*% Theta[[i]] )) - log(det(Theta[[i]])) ) + 2*E[i] 
      })
    ))
  } else if (type[1] == "BIC" || type[1] == "bic") { 
    warning("BIC needs to be implemented") 
    return(sum(
      sapply(1:m, function(i) { 
        n_obs[i] * ( sum(diag( Sigma[[i]] %*% Theta[[i]] )) - log(det(Theta[[i]])) ) + 2*E[i] 
      })
    )) 
  } else { 
    stop("type should be 'AIC' or 'BIC'") 
  }
}