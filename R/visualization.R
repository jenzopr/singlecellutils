#' Plots a silhouette plot
#'
#' @param object A SinglecellExperiment object.
#' @param use_dimred A character string indicating which dimension reduction to use.
#' @param clusters A character string indicating which annotation to use as clusters.
#' @param na.rm Remove NA values from clusters.
#'
#' @return A ggplot object
#'
#' @importFrom rlang .data
#' @export
plotSilhouette <- function(object, use_dimred, clusters, na.rm = TRUE) {
  if ( !methods::is(object, "SingleCellExperiment") ) {
    stop("Object must be of class SingleCellExperiment")
  }

  if (!use_dimred %in% SingleCellExperiment::reducedDimNames(object)) {
    stop(paste("Object must contain a reducedDim named", use_dimred))
  }

  if (!clusters %in% colnames(SummarizedExperiment::colData(object))) {
    stop(paste("Object must contain a column", clusters, "in colData"))
  }

  if (!is.factor(SummarizedExperiment::colData(object)[, clusters])) {
    cl <- factor(SummarizedExperiment::colData(object)[, clusters])
  } else {
    cl <- SummarizedExperiment::colData(object)[, clusters]
  }

  if (na.rm) {
    object <- object[, !is.na(cl)]
    cl <- cl[!is.na(cl)]
  }

  s <- cluster::silhouette(as.numeric(cl), dist = stats::dist(SingleCellExperiment::reducedDim(object, use_dimred)))
  df <- data.frame(cell = factor(colnames(object), levels = colnames(object)), silhouette = s[, "sil_width"], cluster = factor(s[,"cluster"], levels = unique(cl[order(cl)]), ordered = T))

  df$cell <- factor(df$cell, levels = df$cell[order(df$cluster, df$silhouette)])

  ggplot2::ggplot(data = df, ggplot2::aes(.data$cell, .data$silhouette, color = .data$cluster, fill = .data$cluster)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::coord_flip()
}
