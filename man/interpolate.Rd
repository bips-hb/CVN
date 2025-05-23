% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interpolate.R
\name{interpolate}
\alias{interpolate}
\title{Interpolation of a Graph}
\usage{
interpolate(cvn, weights, truncate = NULL)
}
\arguments{
\item{cvn}{A CVN fit with \eqn{m} graphs}

\item{weights}{A vector of length \eqn{m} with the regression coefficients}

\item{truncate}{Truncation value. When a value in the precision matrix is
considered 0. If \code{NULL}, the same truncation is
used as for the fitted CVN model (Default)}
}
\value{
A 'cvn_interpolated' object, which is a list with
\item{\code{adj_matrices}}{A list of adjacency matrix. One for each pair
of \eqn{(\lambda_1, \lambda_2)} values.
The entries are \code{1} if there is an edge, \code{0} otherwise.
The matrices are sparse using package \code{\link[Matrix]{Matrix}}}
\item{\code{m}}{Number of graphs}
\item{\code{p}}{Number of variables}
\item{\code{weights}}{The weights used for interpolation}
\item{\code{truncate}}{Truncation value}
\item{\code{n_lambda_values}}{Total number of \eqn{(\lambda_1, \lambda_2)} value combinations}
\code{results}. It consists of two columns:
\item{\code{lambda1}}{\eqn{\lambda_1} value}
\item{\code{lambda2}}{\eqn{\lambda_2} value}
}
\description{
Estimates a graph for which there are no observation based on
a previously fitted CVN model
}
\examples{
\donttest{
path <- system.file("cvnfit.rda", package = "CVN")
load(path)

interpolate <- interpolate(fit, c(0,0,0,0,0,0,0,0.5,0.5), truncate = 0.05)
fit10 <- combine_cvn_interpolated(fit, interpolate)
fit10 <- visnetwork_cvn(fit10)


if (requireNamespace("htmltools", quietly = TRUE)) {
htmltools::browsable(htmltools::tagList(
  htmltools::div(style = "display: grid; grid-template-columns: repeat(3, 1fr); 
                 column-gap: 1px; row-gap: -10px",
    fit10$plots[[1]][8:10])))
} else {
 fit10$plots[[1]][8:10]
} 
}
             
}
