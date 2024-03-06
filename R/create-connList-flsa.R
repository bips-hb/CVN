#' Create a connListObj for FLSA package
#' 
#' This function creates a connListObj for the FLSA package given a weight matrix.
#' 
#' @param W A weight matrix.
#' @return A connListObj for FLSA package and a logical value indicating whether 
#'         all elements of connList are NULL.
#' @export
create_connListObj <- function(W) {
  m <- nrow(W)
  
  connList <- vector("list", m)
  class(connList) <- "connListObj"
  
  lapply(1:m, function(i) {
    indices <- which(W[i, ] != 0)
    if (length(indices) != 0) {
      connList[[i]] <<- as.integer(which(W[i, ] != 0) - 1)
    }
  })
  
  names(connList) <- as.character(0:(m-1))
  
  connList_is_null <- all(sapply(connList, function(l) is.null(l)))
  
  return(list(connList = connList, connList_is_null = connList_is_null))
}