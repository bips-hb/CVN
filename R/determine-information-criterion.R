#' @references 
#' Danaher, P., Wang, P., & Witten, D. M. (2014). The joint graphical 
#' lasso for inverse covariance estimation across multiple classes. 
#' Journal of the Royal Statistical Society. Series B: Statistical 
#' Methodology, 76(2), 373–397. https://doi.org/10.1111/rssb.12033
#' @keywords internal
determine_information_criterion <- function(Theta, adj_matrices, Sigma, n_obs,
                                            type = c("AIC", "BIC")) { 
  
  m <- length(Theta)             # number of graphs
  E <- mapply(sum, adj_matrices) # number of non-zero edges

  # Scale Matrices for calculating determinant
  T_scaled <- lapply(Theta, function(x) x / max(abs(x)))
  Sigma_scaled <- lapply(Sigma, function(x) x / max(abs(x)))

  if (type[1] == "AIC" || type[1] == "aic") {
    return(sum(
      sapply(1:m, function(i) { 
        as.numeric(n_obs[i] * ( sum(diag( Sigma[[i]] %*% Theta[[i]] ))) - n_obs[i] * determinant(T_scaled[[i]], logarithm = TRUE)$modulus + 2*E[i])
      })
    ))
  } else if (type[1] == "BIC" || type[1] == "bic") { 
    return(sum(
      sapply(1:m, function(i) { 
        as.numeric(n_obs[i] * ( sum(diag( Sigma[[i]] %*% Theta[[i]] ))) - n_obs[i] * determinant(T_scaled[[i]], logarithm = TRUE)$modulus + 2*log(n_obs[i]) * E[i])
      })
    )) 
  } else { 
    stop("type should be 'AIC' or 'BIC'") 
  }
}

#' Information Criteria for a \code{cvn} object
#' 
#' Determines a given information criteria for a \code{cvn} object, 
#' see \code{\link{CVN}}. 
#' 
#' @keywords internal
determine_information_criterion_cvn <- function(cvn, type = c("AIC", "BIC")) { 
  if (!("cvn" %in% class(cvn))) { 
    stop("cvn parameter must be of type 'cvn'") 
  }
  
  # compute for each different lambda1 and lambda2 pair
  # the information criterion 
  sapply(1:cvn$n_lambda_values, function(i) { 
       determine_information_criterion(cvn$Theta[[i]], 
                                       cvn$adj_matrices[[i]], 
                                       cvn$Sigma, 
                                       cvn$n_obs,
                                       type = type)
    })
}
