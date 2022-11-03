
generator <- function(m, w) { 
  return(matrix(rep(w, m*m), ncol = m)) 
}

param <- expand.grid(
  m = 9,
  w = 1:10 / 10
)



#' 
#' @export 
parameterized_weight_matrix <- function(generator, param) { 
  
  if (!("m" %in% names(param))) { 
    stop("m must be a parameter in the param list") 
  }
  
  weight_matrix <- list(
    generator = generator, 
    param = param 
  )
  
  class(weight_matrix) <- "parameterizedW"
  return(weight_matrix)
}

#' Print Function for the Parameterized Weight Matrix Object
#'
#' @export
print.parameterizedW <- function(weight_matrix, ...) { 
  cat("Parameterized Weight Matrix\n\n")
  cat("Generator function ----------\n")
  print(weight_matrix$generator) 
  cat("\nParameters ------------------\n")
  print(weight_matrix$param)
  cat("\n")
}