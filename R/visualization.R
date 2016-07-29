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

#' A function to color a tSNE map based on a vector.
#'
#' @param tsne The data matrix from a tSNE dimension reduction call
#' @param data A numerical vector on which the coloring is based
#' @param use.pal The RColorBrewer color palette
#'
#' @return NULL
#'
#' @export
colortSNE <- function(tsne, data, use.pal = "RdBu") {
  if(length(data) != nrow(tsne)) {
    stop("Data is not same length as tsne.")
  }
  pal <- colorRampPalette(RColorBrewer::brewer.pal(10,use.pal))
  color <- rev(pal(100))[as.numeric(cut(data,breaks=100))]
  plot(tsne[,1],tsne[,2],xlab = "tSNE dimension 1", ylab = "tSNE dimension 2", col = color, pch = 20)
  return(invisible(NULL))
}

#' Creates a sorted bar plot.
#'
#' @param data The data matrix or data.frame.
#' @param marker A numerical index or rowname.
#' @param decreasing Whether the expression values should be sorted decreasing.
#' @param cell.names A vector of alternative cell names.
#' @param cell.groups A factor of cell groups. If specified, cells will be colored accordingly.
#' @param threshold A lower bound for expression. Values below will be shown as 10^-2.
#' @param y.unit A character giving the unit of expression.
#'
#' @return A ggplot2 figure.
#'
#' @export
visMarkerBar <- function(data, marker, decreasing = TRUE, cell.names = NULL, cell.groups = NULL, threshold = 10^-2, y.unit = NULL) {
  if(!is.null(cell.names) & ncol(data) == length(cell.names)) {
    colnames(data) <- cell.names
  }
  if(!is.null(y.unit)) {
    y.label <- paste("(",y.unit,")",sep="")
  } else {
    y.label <- ""
  }
  if(is.null(cell.groups)) {
    cell.groups <- rep("none", ncol(data))
  }
  e <- data[marker,]
  if(!is.null(e)) {
    order <- order(e, decreasing = decreasing)
    data <- data.frame(cell = factor(colnames(e), levels = colnames(e)[order]),
                       expr = ifelse(as.numeric(e)<threshold, 10^-2, as.numeric(e)),
                       group = cell.groups)
    p <- ggplot2::ggplot(data, aes(x=cell,y=expr, fill=group)) +
      ggplot2::geom_bar(stat="identity") +
      ggplot2::scale_fill_brewer(type = "qual", palette = "Paired") +
      ggplot2::ylab(paste("Expression level",y.label)) +
      ggplot2::xlab("Individual single-cells sorted by expression") +
      ggplot2::ggtitle(marker) +
      ggplot2::theme(axis.text.x  = element_text(angle=90,hjust=0.5,vjust=0.5,size=8), panel.background=element_rect(fill = "#FFFFFF"))
    return(p)
  }
  else {
    warning(paste("Data did not contain a row",marker))
    return(NULL)
  }
}
