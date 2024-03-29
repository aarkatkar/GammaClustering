% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gcluster.R
\name{gcluster}
\alias{gcluster}
\title{Perform gamma clustering}
\usage{
gcluster(gm, n, maxiter = 1000, method = "affinity")
}
\arguments{
\item{gm}{A gmatrix object}

\item{n}{The maximum number of clusters to find}

\item{maxiter}{The maximum number of iterations of the algorithm}

\item{method}{Use the "distance" or "affinity" matrix}
}
\value{
The results of gamma clustering
}
\description{
This function loads a file as a matrix. It assumes that the first column
contains the rownames and the subsequent columns are the sample identifiers.
Any rows with duplicated row names will be dropped with the first one being
kepted.
}
