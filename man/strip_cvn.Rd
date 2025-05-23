% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/strip_cvn.R
\name{strip_cvn}
\alias{strip_cvn}
\title{Strip CVN}
\usage{
strip_cvn(cvn)
}
\arguments{
\item{cvn}{Object of class 'CVN'}
}
\value{
Reduced CVN where \code{Theta}, \code{data} and \code{Sigma}
are removed
}
\description{
Function that removes most of the items to make the CVN object
more memory sufficient. This is especially important when the
graphs are rather larger
}
