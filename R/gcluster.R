#' Perform gamma clustering
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#'
#' @param gm A gmatrix object
#' @param n The maximum number of clusters to find
#' @param maxiter The maximum number of iterations of the algorithm
#' @param method Use the "distance" or "affinity" matrix
#' @return The results of gamma clustering
#' @export
gcluster <- function(gm,n, maxiter=1000, method="affinity"){
  N <- length(gm$affinity)^0.5
  bld <- cluster::pam(gm$distance, n, diss=TRUE, do.swap=FALSE)
  selected <- bld$medoids

  m <- matrix(rep(0, N*n), nrow=n)
  for (i in 1:n){
    m[i,selected[i]]=1
  }
  clusters <- apply(m, 2, which.max)

  d <- m %*% gm$affinity
  upper <- exp(d-max(d))
  total <- rowSums(t(upper))
  m <- t(t(upper)/total)

  if (method=="affinity"){
    pfunc <- FALSE
    pmat <- gm$affinity
  }else{
    pfunc <- TRUE
    pmat <- gm$distance
  }

  for (i in 1:maxiter){
    d <- m %*% pmat
    total <- rowSums(m)-m
    p0 <- pgamma(d, total, lower.tail=pfunc)
    p <- -log(ifelse(p0 <=1e-300,1e-300,p0))
    upper <- exp(p-max(p))
    total <- rowSums(t(upper))
    m <- t(t(upper)/total)
    clusters0 <- apply(m, 2, which.max)
    if (all(clusters0 == clusters) & i > 10){
      break
    }
    clusters <- clusters0
  }


  a <- apply(m, 2, max)
  b <- c(apply(p, 2, function(p){sort(p)[length(p)-2]}))
  score <- a
  return(list(labels=clusters, inertia=mean(score)))
}
