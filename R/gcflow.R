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
#' @return The maximum information flow path connecting the phenotypes
#' @export
gc_max_flow <- function(n0, n1, gm, clusters){
  s <- which(gm$geneNames==n0)
  t <- which(gm$geneNames==n1)

  g1 <- igraph::graph_from_adjacency_matrix(gm$pvals, weighted=TRUE, mode="directed", add.colnames=NA)
  g2 <- igraph::graph_from_adjacency_matrix(gm$oCor, weighted=TRUE,mode="directed", add.colnames=NA)
  g3 <- igraph::graph_from_adjacency_matrix(gm$odist, weighted=TRUE)
  cutoff <- 1-exp(-max(igraph::E(g3)[igraph::shortest_paths(g3, s, t, output="epath")$epath[[1]]]$weight))

  gcweights <- igraph::E(g1)$weight
  igraph::edge_attr(g1, "capacity") <- -0.5*log(1-igraph::E(g2)$weight)
  igraph::edge_attr(g1, "weight") <- 1
  drop <- which(gcweights > cutoff)
  g1_filter <- igraph::delete_edges(g1, igraph::E(g1)[drop])
  mflow <- gcflow(g1_filter, s, t)$flows
  drop <- which(mflow==0)
  igraph::edge_attr(g1_filter, "weight") <- mflow
  g1_flow <- igraph::delete_edges(g1_filter, igraph::E(g1_filter)[drop])
  gcvertices <- igraph::degree(g1_flow)
  drop <- which(gcvertices == 0)
  gflow <- igraph::delete_vertices(g1_flow, igraph::V(g1_flow)[drop])
  drop <- which(gcvertices != 0)
  n_phen <- length(gm$geneNames)-length(clusters)
  hpal <- c("white", scales::hue_pal()(max(clusters)))
  spal <- c("none", rep("circle", max(clusters)))
  geneNames <- gm$geneNames[drop]
  geneClusters <- c(rep(0, n_phen),clusters)[drop]+1
  plot(gflow, layout=igraph::layout_as_tree, vertex.label=geneNames,
       vertex.color=hpal[geneClusters], edge.label=round(mflow[which(mflow!=0)],3),
       edge.arrow.size=0.5,vertex.shape=spal[geneClusters])
}

gcflow <- function(g, v, to){
  igraph::edge_attr(g, "flow") <- 0
  flow <- 0
  # BFS
  while(TRUE){
    drop <- which(igraph::E(g)$flow >= igraph::E(g)$capacity)
    g0 <- igraph::delete_edges(g, drop)
    drop <- which(igraph::E(g)$flow < igraph::E(g)$capacity)
    edrop <- igraph::E(g)[drop]
    sp <- igraph::shortest_paths(g0, v, to, weights=NA,output="epath", mode="out")
    aug_edges <- sp$epath[[1]]
    if(length(aug_edges)==0){
      break
    }
    rem <- igraph::E(g0)$capacity[aug_edges]-igraph::E(g0)$flow[aug_edges]
    igraph::E(g0)$flow[aug_edges] <- min(rem) + igraph::E(g0)$flow[aug_edges]
    igraph::E(g)$flow[edrop] <- igraph::E(g0)$flow
    flow<- flow + min(rem)
  }
  return(list(flow=flow, flows=igraph::E(g)$flow))
}
