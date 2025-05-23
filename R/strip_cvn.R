#' Strip CVN
#'
#' Function that removes most of the items to make the CVN object
#' more memory sufficient. This is especially important when the
#' graphs are rather larger
#'
#' @param cvn Object of class 'CVN'
#'
#' @return Reduced CVN where \code{Theta}, \code{data} and \code{Sigma}
#'         are removed
#' @export
strip_cvn <- function(cvn) {

  if ('Theta' %in% names(cvn)) {
    cvn <- within(cvn, rm(Theta))
  }

  if ('data' %in% names(cvn)) {
    cvn <- within(cvn, rm(data))
  }

  if ('Sigma' %in% names(cvn)) {
    cvn <- within(cvn, rm(Sigma))
  }

  # set the variable keeping track of whether the cvn is
  # striped to TRUE
  cvn$minimal <- TRUE

  return(cvn)
}

