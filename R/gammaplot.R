#' Visualize gamma clustering
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#'
#' @param gm A gmatrix object
#' @param clusters A list of length(n_genes) labeling each gene with a cluster
#' @return A UMAP plot color-coded by cluster membership
#' @export
gammaplot <- function(gm, clusters){
  colorList <- scales::hue_pal()(max(clusters))
  u <- as.data.frame(umap::umap(gm$distance, input="dist")$layout)
  colnames(u) <- c("UMAP1", "UMAP2")
  u$clusters <- as.factor(clusters)
  ggplot2::ggplot() + ggplot2::coord_cartesian() +
    ggplot2::scale_x_continuous() +
    ggplot2::scale_y_continuous() +
    ggplot2::scale_color_hue() +
    ggplot2::layer(data=u,
                   mapping=ggplot2::aes(x=UMAP1, y=UMAP2, color=clusters),
                   stat="identity",
                   geom="point",
                   position=ggplot2::position_jitter()
    )
}
