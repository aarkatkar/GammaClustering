% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gcpath.R
\name{gc_path}
\alias{gc_path}
\title{Perform gamma clustering}
\usage{
gc_path(n0, n1, gm, clusters)
}
\arguments{
\item{n0}{The starting phenotype}

\item{n1}{The ending phenotype}

\item{gm}{A gmatrix object}

\item{clusters}{Cluster labels obtained from gcluster}
}
\value{
The maximum likelihood path connecting the phenotypes
}
\description{
This function loads a file as a matrix. It assumes that the first column
contains the rownames and the subsequent columns are the sample identifiers.
Any rows with duplicated row names will be dropped with the first one being
kepted.
}
