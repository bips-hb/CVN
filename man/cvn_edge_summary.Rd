% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find-core-graph.R
\name{cvn_edge_summary}
\alias{cvn_edge_summary}
\title{Edge Overview}
\usage{
cvn_edge_summary(cvn)
}
\arguments{
\item{cvn}{A \code{cvn} object}
}
\value{
A data.frame showing the number of edges in each subgraph, the number of
edges present in all graphs (core edges) and the number of edges that are unique
for a subgraph
}
\description{
Function that gives an overview about the number of edges in each subgraph over all
fitted CVN models.
}
\examples{
data(grid)
W <- create_weight_matrix("grid", 3, 3)

# lambdas:
lambda1 = c(1, 1.5)
lambda2 = .2

fit <- CVN(grid, W, lambda1 = lambda1, lambda2 = lambda2, n_cores = 1)

# Edge summary for a list of CVN models
cvn_edge_summary(fit)

# Edge summary for a single CVN
fit2 <- extract_cvn(fit, id = 2)
cvn_edge_summary(fit2)                 

}
