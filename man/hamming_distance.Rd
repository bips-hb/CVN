% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hamming-distance.R
\name{hamming_distance}
\alias{hamming_distance}
\title{Structural Hamming Distance for a \code{cvn} Object}
\usage{
hamming_distance(cvn, verbose = TRUE)
}
\arguments{
\item{cvn}{A \code{cvn} object
created by the \code{\link[CVN]{CVN}} function}

\item{verbose}{If \code{TRUE}, shows a progress bar}
}
\value{
A list of symmetric matrices. Each matrix contains the structural
Hamming distances between the different graphs. Each item in the
list corresponds to one \eqn{(\lambda_1, \lambda_2)} pair
}
\description{
Returns the structural Hamming distances
}
\seealso{
\code{\link[CVN]{CVN}}
}
