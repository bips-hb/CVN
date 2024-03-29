% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot-graph.R
\name{create_edges_visnetwork}
\alias{create_edges_visnetwork}
\title{Create a \code{data.frame} for the Edges for \code{visNetwork}}
\usage{
create_edges_visnetwork(adj_matrix)
}
\arguments{
\item{adj_matrix}{A symmetric adjacency matrix}
}
\value{
Data frame that be used as input for \code{visNetwork}
}
\description{
In order to visualize a graph, we need to create a 
\code{data.frame} that can be used by the \code{visNetwork} package. 
This function returns the needed \code{data.frame} given a
adjacency matrix.
}
