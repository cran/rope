% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rope.R
\name{plotrope}
\alias{plotrope}
\title{Plot rope results}
\usage{
plotrope(result, data, types = c("global"), ...)
}
\arguments{
\item{result}{An object returned by \code{rope} or \code{explore}}

\item{data}{Matrix of variable presence counts. One column for each variable,
one row for each parameter value (e.g. levels of regularization).}

\item{types}{List of names of plots to draw (alternatives \code{'global'},
\code{'q-values'} or \code{'fits'})}

\item{...}{Pass level=v for a vector v of indices when drawing the fits plot
to only plot for penalization levels corresponding to v}
}
\description{
Plot rope results
}
