#' Visualizes tSNE maps with colored clusters by kmeans.
#'
#' @param tsne The data matrix from a tSNE dimension reduction call
#' @param k The number of clusters. If NULL, will perform Gap statistics.
#' @param use.pal The RColorBrewer palette name to use for cluster colors.
#' @param add.center Whether cluster centers should be added to the plot.
#'
#' @return The augmented input data matrix. Clusters are added in .cluster column.
#'
#' @export
vistSNE <- function(tsne, k = NULL, use.pal = "Dark2", add.center = T) {
  if(is.null(k)) {
    cg <- cluster::clusGap(tsne, kmeans, 10, d.power = 2)
    k <- maxSE(cg$Tab[,3],cg$Tab[,4],method="firstSEmax",.25)
  }
  if(k < 2) {
    stop("Cannot cluster tSNE into less than 2 clusters. Set k to at least 2.")
  }
  km <- kmeans(tsne, centers = k)
  pal <- RColorBrewer::brewer.pal(k, use.pal)

  plot(tsne[,1],tsne[,2],xlab = "tSNE dimension 1", ylab = "tSNE dimension 2", col = pal[km$cluster], pch = 20)
  if(add.center) {
    points(km$centers[,1],km$centers[,2], pch = 25, col = pal)
    text(km$centers[,1]+0.05,km$centers[,2], labels = 1:k)
  }

  return(invisible(broom::augment(km, tsne)))
}
