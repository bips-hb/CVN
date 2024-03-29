% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create-matrix-D.R
\name{create_matrix_D}
\alias{create_matrix_D}
\title{Create matrix \eqn{D} to be used for the Generalized LASSO}
\usage{
create_matrix_D(W, lambda1, lambda2, rho = 1, remove_zero_row = TRUE)
}
\arguments{
\item{W}{The \eqn{(m \times m)}-dimensional upper-triangular 
weight matrix \eqn{W}}

\item{lambda1}{The \eqn{\lambda_1} LASSO penalty term}

\item{lambda2}{The \eqn{\lambda_2} global smoothing parameter}

\item{rho}{The \eqn{\rho} ADMM's penalty parameter (Default: \code{1})}

\item{remove_zero_row}{If \code{TRUE}, rows with zeros are removed. 
(Default: \code{TRUE})}
}
\value{
A \eqn{((m \cdot (m+1)/2) \times m)}-dimensional matrix
}
\description{
Generates a matrix \eqn{D} to be used for the generalized LASSO. 
We solve a generalized LASSO problem for each edge \eqn{(s,t)} 
for each update step for \eqn{Z}.
}
\examples{
m <- 4 # number of graphs
W <- matrix(1, nrow = m, ncol = m) 

# penalty terms:
lambda1 <- .2
lambda2 <- .4
rho <- 1

CVN::create_matrix_D(W, lambda1, lambda2, rho) 
#      [,1] [,2] [,3] [,4]
# [1,]  0.2  0.0  0.0  0.0
# [2,]  0.0  0.2  0.0  0.0
# [3,]  0.0  0.0  0.2  0.0
# [4,]  0.0  0.0  0.0  0.2
# [5,]  0.4 -0.4  0.0  0.0
# [6,]  0.4  0.0 -0.4  0.0
# [7,]  0.4  0.0  0.0 -0.4
# [8,]  0.0  0.4 -0.4  0.0
# [9,]  0.0  0.4  0.0 -0.4
# [10,]  0.0  0.0  0.4 -0.4
}
\references{
Tibshirani, R. J., & Taylor, J. (2011). 
            The solution path of the generalized lasso. 
            Annals of Statistics, 39(3), 1335–1371. 
            https://doi.org/10.1214/11-AOS878
}
