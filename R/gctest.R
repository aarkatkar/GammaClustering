#' Perform gamma clustering
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#'
#' @param n_initial The lower bound number of clusters
#' @param n_final The upper bound number of clusters
#' @param gm A gmatrix object
#' @param method Use the "distance" or "affinity" matrix
#' @return A plot of the mean confidence for each cluster number
#' @export
gctest <- function(n_initial, n_final, gm,method="affinity"){
  inertia <- c()
  for (n in c(n_initial:n_final)){
    clustered <- gcluster(gm, n, method=method)
    score <- clustered$inertia
    inertia<- c(inertia,score)
  }
  n <- c(n_initial:n_final)
  plot(n, inertia)
}
