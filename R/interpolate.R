#' Interpolation of a Graph
#' 
#' Estimates a graph for which there are no observation based on 
#' a previously fitted CVN model
#' 
#' @param cvn A CVN fit with \eqn{m} graphs
#' @param weights A vector of length \eqn{m} with the regression coefficients
#' @param truncate Truncation value. When a value in the precision matrix is 
#'                 considered 0. If \code{NULL}, the same truncation is 
#'                 used as for the fitted CVN model (Default)
#' 
#' @return A 'cvn_interpolated' object, which is a list with 
#'    \item{\code{adj_matrices}}{A list of adjacency matrix. One for each pair 
#'                               of \eqn{(\lambda_1, \lambda_2)} values. 
#'                               The entries are \code{1} if there is an edge, \code{0} otherwise. 
#'                               The matrices are sparse using package \code{\link[Matrix]{Matrix}}} 
#'   \item{\code{m}}{Number of graphs}
#'   \item{\code{p}}{Number of variables}
#'   \item{\code{weights}}{The weights used for interpolation}
#'   \item{\code{truncate}}{Truncation value}
#'   \item{\code{n_lambda_values}}{Total number of \eqn{(\lambda_1, \lambda_2)} value combinations}
#'   \code{results}. It consists of two columns:
#'   \item{\code{lambda1}}{\eqn{\lambda_1} value}
#'   \item{\code{lambda2}}{\eqn{\lambda_2} value}
#'   
#' @export
interpolate <- function(cvn, weights, truncate = NULL) { 
  
  
  # Check correctness input -----------------------------------------
  if (!("cvn" %in% class(cvn))) { 
    stop("cvn should be of the type 'cvn'")  
  }
  
  if (!("Theta" %in% names(cvn))) { 
    stop("The estimated Theta matrices are missing. Did you use the option minimal while fitting?")  
  }
  
  if (!is.vector(weights) || !is.numeric(weights)) { 
    stop("weights should be a numeric vector")
  }
  
  if (length(weights) != cvn$m) { 
    stop("weights should be a vector of length m")
  }

  # If truncate value is not specified, same is used as for the CVN fit
  if (is.null(truncate)) { 
    truncate <- cvn$truncate  
  }
  
  # Function for interpolation -------------------------------------------------
  inter <- function(x, y, w, lambda1, lambda2) { 
    lambda1 * abs(x) + lambda2 * sum(abs(w*y - x)) 
  }

  # Go over all lambda1 and lambda2 pairs --------------------------------------
  adj_matrices <- lapply(1:cvn$n_lambda_values, function(i) {
    lambda1 <- cvn$results$lambda1[i]
    lambda2 <- cvn$results$lambda2[i]
    
    # initial the matrix
    adj_matrix <- Matrix(0, nrow = cvn$p, ncol = cvn$p, sparse = TRUE)

    # go over all potential edges (s,t)
    combn(cvn$p, 2, function(pair) {
      y <- sapply(cvn$Theta[[i]], function(Theta_i) Theta_i[pair[1],pair[2]])
      
      # Check if y is numeric
      if (!is.numeric(y)) {
        stop("y must be a numeric vector")
      }      
      
      fit <- optim(
        par = 0,
        fn = inter,
        y = y,
        w = weights,
        lambda1 = lambda1,
        lambda2 = lambda2,
        method = "Brent",
        lower = min(weights * y) - 1,
        upper = max(weights * y) + 1
      )

      edge_exists <- abs(fit$value) >= truncate
      adj_matrix[pair[1], pair[2]] <<- as.numeric(edge_exists)
      adj_matrix[pair[2], pair[1]] <<- as.numeric(edge_exists)
    })
    return(adj_matrix)
  })

  cvn_interpolated <- 
    list(
         adj_matrices = adj_matrices,
         m = cvn$m, 
         p = cvn$p, 
         W = cvn$W,
         weights = weights, 
         truncate = truncate, 
         n_lambda_values = cvn$n_lambda_values, 
         results = cvn$results %>% select(lambda1,lambda2)
  )
  
  # Set the class back to "irgendwas"
  class(cvn_interpolated) <- "cvn_interpolated"
  
  return(cvn_interpolated)
}
  
  