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

#' Creating a ComplexHeatmap from expression or reducedDim data of a \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#'
#' This function serve as a wrapper for the \code{\link[ComplexHeatmap]{Heatmap}} function and uses \code{features} to show expression values or columns from the \code{use_dimred} slot.
#'
#' @param object A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object.
#' @param features A character vector, vector of indices or a named list thereof. In case of a list, rows are split accordingly.
#' @param exprs_values String indicating which assay contains the data that should be used for plotting.
#' @param use_dimred  A character string indicating which dimension reduction to use, or NULL.
#' @param split_by Character vector indicating by which columns should be split. In case of length of one, it determines a column name of the \code{\link[SummarizedExperiment]{colData}} slot.
#' @param rownames A string indicating a column name of \code{\link[SummarizedExperiment]{rowData}} slot, used as alternative rownames.
#' @param scale A logical, indicating whether data should be scaled before plotting.
#' @param col A vector of colors if the color mapping is discrete or a color mapping function if the matrix is continuous. See \code{\link[ComplexHeatmap]{Heatmap}}.
#' @param ... Additional arguments passed on to the \code{\link[ComplexHeatmap]{Heatmap}} function.
#'
#' @return A \code{\link[ComplexHeatmap]{Heatmap-class}} object.
#'
#' @export
plotComplexHeatmap <- function(object, features, exprs_values = "normcounts", use_dimred = NULL, split_by = NULL, rownames = NULL, scale = FALSE, col = NULL, ...) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Package ComplexHeatmap needed for this function to work. Please install it.", call. = FALSE)
  }

  # If a list of features, enable split
  heatmap_split <- NULL
  if(is.list(features)) {
    # Make sure features is a named list
    if(is.null(names(features))) {
      names(features) <- as.character(1:length(features))
    }
    heatmap_split <- rep(names(features), sapply(features, length))
    features <- unlist(features)

    # Remove duplicated entries
    duplicated_features <- duplicated(features)
    if(any(duplicated_features)) {
      warning(paste(sum(duplicated_features), "duplicated features have been removed from features."))
    }

    features <- features[!duplicated_features]
    heatmap_split <- heatmap_split[!duplicated_features]
  }

  # Plot expression values or reduced dim?
  if(!is.null(exprs_values)) {
    assertive.sets::is_subset(exprs_values, SummarizedExperiment::assayNames(object))

    if(is.character(features)) {
      feature_match <- match(features, rownames(object))
      object <- object[na.omit(feature_match), ]
    } else {
      object <- object[features, ]
    }

    object %>%
      SummarizedExperiment::assay(i = exprs_values) %>%
      as.matrix() -> data

    # Add symbols instead of gene IDs
    if(!is.null(rownames) && rownames %in% colnames(SummarizedExperiment::rowData(object))) {
      rownames(data) <- SummarizedExperiment::rowData(object)[, rownames]
    }
  } else {
    if(!is.null(use_dimred)) {
      assertive.sets::is_subset(use_dimred, SingleCellExperiment::reducedDimNames(object))

      if(is.character(features)) {
        feature_match <- features
      } else {
        feature_match <- colnames(SingleCellExperiment::reducedDim(object, use_dimred))[features]
      }

      object %>%
        SingleCellExperiment::reducedDim(use_dimred) %>%
        as.data.frame() %>%
        dplyr::select_(dots = feature_match) %>%
        as.matrix() %>%
        t() -> data
    } else {
      stop("Both, exprs_values and use_dimred, cannot be NULL.")
    }
  }

  # Scale to z-score
  if(scale) {
    data <- t(scale(t(data), scale = F))
    heatmap_color <- circlize::colorRamp2(breaks = seq(from = -10, to = 10, length.out = 9), colors = rev(RColorBrewer::brewer.pal(9, "RdBu")))
  } else {
    heatmap_color <- circlize::colorRamp2(breaks = seq(from = 0, to = max(data), length.out = 99), colors = viridis::magma(99))
  }

  # Heatmap color handling
  if(is.null(col)) {
    col <- heatmap_color
  }

  ComplexHeatmap::Heatmap(matrix = data,
                          col = col,
                          split = heatmap_split,
                          ...)
}
