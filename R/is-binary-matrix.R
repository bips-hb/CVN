#' Check if a matrix is binary
#'
#' This function checks whether a given matrix is binary, i.e., all its elements are either 0 or 1.
#'
#' @param mat A matrix to be tested.
#'
#' @return A logical value indicating whether the matrix is binary or not.
#'
#' @examples
#' binary_mat <- matrix(c(0, 1, 1, 0, 1, 0, 0, 1), nrow = 2)
#' is_binary_matrix(binary_mat)
#'
#' @export
is_binary_matrix <- function(mat) {
  if (!is.matrix(mat)) {
    stop("Input must be a matrix")
  }
  
  if (any(mat != 0 & mat != 1)) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}
