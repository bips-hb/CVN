#' Hitting the end points of \eqn{(\lambda_1, \lambda_2)} interval
#' 
#' One often selected the optimal model for the \eqn{(\lambda_1, \lambda_2)}-values
#' based on the AIC and BIC. 
#' This function checks and warns when the optimal value lies on the border of the 
#' values \eqn{(\lambda_1, \lambda_2)} takes. 
#' 
#' @param results Results of the \code{\link{CVN}} function
#' 
#' @return List with two values: 
#'    \item{\code{hits_border_aic}}{If \code{TRUE}, hits the border for the AIC}
#'    \item{\code{hits_border_bic}}{If \code{TRUE}, hits the border for the BIC}
#' @export
hits_end_lambda_intervals <- function(results) { 
  
  hits_border_aic <- FALSE
  hits_border_bic <- FALSE
  
  lambda1 <- results$lambda1
  lambda2 <- results$lambda2
  
  # if there are less than 3 entries for both lambda1 and lambda2, there is no point
  # in checking whether the optimal estimate based on the AIC or the BIC lies on 
  # the edge of the parameter search space
  if(length(lambda1) >= 3 && length(lambda2) >= 3) { 
    
    # write the standard message and fill in the gaps later
    message_warning <- function(criterion = c("BIC", "AIC"), 
                                lambda = c("lambda1", "lambda2"), 
                                hits = c("smallest", "largest")) {
      sprintf("In case you are selecting your model on the basis of the %s: Note that the maximum %s is achieved at the %s value of the %s interval. We suggest you increase the %s interval", 
              criterion[1], 
              criterion[1], 
              hits[1], 
              lambda[1],
              lambda[1])
    }
    
    # For the BIC score ----------------------------------------------------------
    temp <- cvn$results %>% filter(bic == max(bic))
    
    # if lambda1 hits one of the left border
    if (temp$lambda1 == min(lambda1)) {
      warning(message_warning("BIC", "lambda1", "smallest"))
      hits_border_bic <- TRUE
    }
    
    # if lambda1 hits one of the right border
    if (temp$lambda1 == max(lambda1)) {
      warning(message_warning("BIC", "lambda1", "largest"))
      hits_border_bic <- TRUE
    }
    
    # if lambda1 hits one of the left border
    if (temp$lambda2 == min(lambda2)) {
      warning(message_warning("BIC", "lambda2", "largest"))
      hits_border_bic <- TRUE
    }
    
    # if lambda1 hits one of the left border
    if (temp$lambda2 == max(lambda2)) {
      warning(message_warning("BIC", "lambda2", "largest"))
      hits_border_bic <- TRUE
    }
    
    # For the AIC score ----------------------------------------------------------
    temp <- cvn$results %>% filter(aic == max(aic))
    
    # if lambda1 hits one of the left border
    if (temp$lambda1 == min(lambda1)) {
      warning(message_warning("AIC", "lambda1", "smallest"))
      hits_border_aic <- TRUE
    }
    
    # if lambda1 hits one of the right border
    if (temp$lambda1 == max(lambda1)) {
      warning(message_warning("AIC", "lambda1", "largest"))
      hits_border_aic <- TRUE
    }
    
    # if lambda1 hits one of the left border
    if (temp$lambda2 == min(lambda2)) {
      warning(message_warning("AIC", "lambda2", "largest"))
      hits_border_aic <- TRUE
    }
    
    # if lambda1 hits one of the left border
    if (temp$lambda2 == max(lambda2)) {
      warning(message_warning("AIC", "lambda2", "largest"))
      hits_border_aic <- TRUE
    }
  } 
  
  return(
    list(
      hits_border_aic = hits_border_aic,
      hits_border_bic = hits_border_bic
    )
  )
}