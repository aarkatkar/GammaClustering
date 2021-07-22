#' Load a Matrix
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#'
#' @param x A matrix of fold changes with shape (n_genes, n_samples)
#' @param partial Replace pearson correlations with partial correlations
#' @param squared Square correlation values
#' @param normalization Perform quantile normalization on x
#' @param permuted Iterations of permutation testing for p values
#' @return An object containing processed matrices
#' @export
gmatrix <- function(x, partial=FALSE, squared=TRUE, normalization=FALSE,
                    permuted=0){
  geneNames <- rownames(x)
  if (normalization==TRUE){
    x <- preprocessCore::normalize.quantiles(x)
  }
  oCor <- as.matrix(cor(t(x), method="pearson"))
  if (partial == TRUE){
    corMat <- as.matrix(GeneNet::ggm.estimate.pcor(t(x)))
  } else{
    corMat <- oCor
  }
  n <- dim(x)[2]
  t <- oCor *sqrt((n-2)/(1-oCor^2))
  if (squared==TRUE){
    corMat <- corMat^2
    p <- pt(-abs(t), n-2)*2
  } else{
    p <- pt(t, n-2, lower.tail=FALSE)
  }
  oCor <- oCor^2
  diag(oCor) <-0

  if (permuted>0){
    oCor2 <- gcperm(x, iter=permuted, partial=partial, squared=squared)(corMat)
    tCor <- matrix(oCor2, nrow=dim(corMat)[1])
    p <- 1-tCor
  }

  odist <- -log(1-p)
  diag(odist) <- 0

  diag(corMat)=1.01
  rankMat <- 1-rank(corMat, ties.method="min")/(length(corMat)-length(diag(corMat)))
  N <- dim(corMat)[1]
  dim(rankMat) <- c(N,N)
  diag(rankMat) <- 1
  A <- -log(rankMat+1/(length(corMat)-length(diag(corMat))))
  diag(rankMat) <- 0
  diag(A) <- 0

  D <- -log(1-rankMat)
  diag(corMat) <- 0
  gm <- list(affinity=A, distance=D, odist=odist, pvals=p, corMat=corMat, oCor=oCor, geneNames=geneNames)
  return(gm)
}

gcperm <- function(x, iter=10, partial=FALSE, squared=TRUE){
  sampled <- c()
  if (partial==TRUE){
    for (i in 1:iter){
      pMat <- apply(x, 1, sample)
      pCorMat <- GeneNet::ggm.estimate.pcor(pMat)
      sampled <-c(sampled, pCorMat)
    }

  } else{
    for (i in 1:iter){
      pMat <- apply(x, 1, sample)
      pCorMat <- cor(pMat, method="pearson")
      sampled <-c(sampled, pCorMat)
    }
  }
  if (squared==TRUE){
    sampled <- sampled^2
  }

  return(ecdf(sampled))
}





