#' Perform gamma clustering
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#'
#' @param n0 The starting phenotype
#' @param n1 The ending phenotype
#' @param gm A gmatrix object
#' @param clusters Cluster labels obtained from gcluster
#' @return The maximum likelihood path connecting the phenotypes
#' @export
gc_path <- function(n0, n1, gm, clusters){
  s <- which(gm$geneNames==n0)
  t <- which(gm$geneNames==n1)

  g3 <- igraph::graph_from_adjacency_matrix(gm$odist, weighted=TRUE)
  sp <- igraph::shortest_paths(g3, s, t, output="both")
  sp_vert <- sp$vpath[[1]]
  g3df <- data.frame(from=c(1:(length(sp_vert)-1)), to=c(2:length(sp_vert)))
  g3df <- igraph::graph.data.frame(g3df)

  n_phen <- length(gm$geneNames)-length(clusters)
  hpal <- c("white", scales::hue_pal()(max(clusters)))
  spal <- c("none", rep("circle", max(clusters)))
  geneNames <- gm$geneNames[sp_vert]
  geneClusters <- c(rep(0, n_phen),clusters)[sp_vert]+1
  plot(g3df, layout=igraph::layout_as_tree, vertex.label=geneNames,
       vertex.color=hpal[geneClusters],
       edge.arrow.size=0.5,vertex.shape=spal[geneClusters])
}
