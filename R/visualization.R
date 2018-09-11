#' A function to color observations on a map based on values of a vector.
#'
#' @param data The data matrix with observation coordinates in the first two columns.
#' @param colour_by Either a numerical vector on which the coloring is based, or a numeric or character of length 1 indicating the column of the data matrix from which to take the values, or a factor.
#' @param shape_by An optional factor on which point shape is based.
#' @param palette The color palette.
#' @param na.colour The color of \code{NA}s in \code{colour_by}.
#' @param na.shape The shape of \code{NA}s in \code{shape_by}
#' @param na.cex The cex of \code{NA} observations.
#' @param xlab,ylab Labels of x and y axis.
#' @param pch Default point shape.
#' @param cex Default point size.
#' @param ... Other parameters passed to lattice::xyplot
#'
#' @return A lattice xyplot. In case \code{colour_by} or \code{shape_by} is a factor, a legend is drawn.
#'
#' @export
colorAMap <- function(data, colour_by, shape_by = NULL, palette = RColorBrewer::brewer.pal(8, "Dark2"), na.colour = "darkgrey", na.shape = 4, na.cex = 1, xlab = "Dimension 1", ylab = "Dimension 2", pch = rep(16, nrow(data)), cex = rep(1, nrow(data)), ...) {
  if (length(colour_by) != nrow(data) & length(colour_by) > 1) {
    stop("Values are not of same length as data and more than one value given")
  }

  if (length(colour_by) == 1 & ncol(data) <= length(colour_by)) {
    v <- as.numeric(data[, colour_by])
  }

  if (is.factor(colour_by)) {
    levels <- levels(colour_by)
    v <- as.numeric(colour_by)
  } else {
    v <- as.numeric(colour_by)
  }

  if (is.factor(shape_by)) {
    symbols <- c(16, 15, 17, 18, 19, 0:14)
    pch <- symbols[as.numeric(shape_by)]
    legend.pch <- symbols[as.numeric(levels(shape_by))]
  } else {
    if (is.null(pch)) {
      pch <- rep(16, length(v))
    }
    legend.pch <- pch
  }

  if (is.factor(colour_by) & length(unique(colour_by)) <= length(palette)) {
    warning("Less unique data values than colors in palette.. reducing palette.")
    palette <- palette[1:length(unique(stats::na.omit(colour_by)))]
  }
  color <- palette[as.numeric(cut(v, breaks = length(palette)))]

  if (is.null(cex) | length(cex) != length(v)) {
    cex <- rep(1, length(v))
  }

  if (any(is.na(colour_by))) {
    color[is.na(colour_by)] <- na.colour
    pch[is.na(colour_by)] <- na.shape
    cex[is.na(colour_by)] <- na.cex
  }
  if (any(is.na(pch))) {
    pch[is.na(pch)] <- na.shape
  }

  #
  # Construct legend key in case of factor
  #
  color.key <- NULL
  shape.key <- NULL
  legend.key <- list(corner=c(1,1), rep = FALSE)

  if (is.factor(shape_by)) {
    shape.key <- list(points = list(col = "darkgrey", pch = legend.pch),
                      text = list(levels(shape_by)))
    legend.key <- c(legend.key, shape.key)
  }
  if (is.factor(colour_by)) {
    color.key <- list(points = list(col = palette[1:length(levels)], pch = 16),
                      text = list(levels))
    legend.key <- c(legend.key, color.key)
  }

  # If legend is only of length 2 (nothing added), set to NULL
  if (length(legend.key) == 2) {
    legend.key <- NULL
  }

  p <- lattice::xyplot(data[, 2] ~ data[, 1], xlab = xlab, ylab = ylab, col = color, key = legend.key, pch = pch, cex = cex, ...)
  return(p)
}

#' Creates a sorted bar plot of expression values.
#'
#' @param object A SingleCellExperiment object.
#' @param exprs_values String indicating which assay contains the data that should be used to perform thinning.
#' @param feature A numerical index or character vector giving the rowname.
#' @param colour_by Character giving the name of a column in \code{colData} to color by.
#' @param threshold A lower bound for expression. Values below will be shown as 10^-2.
#' @param decreasing Whether the expression values should be sorted decreasing.
#'
#' @return A ggplot2 figure.
#'
#' @importFrom rlang .data
#' @export
visMarkerBar <- function(object, exprs_values, feature, colour_by = NULL, threshold = 10^-2, decreasing = TRUE) {
  e <- SummarizedExperiment::assay(object, i = exprs_values)
  if (!is.null(e)) {
    order <- order(e, decreasing = decreasing)

    if (!is.null(colour_by)) {
      data <- data.frame(cell = factor(colnames(object), levels = colnames(object)[order]),
                         expr = ifelse(as.numeric(e) < threshold, 10^ - 2, as.numeric(e)),
                         group = SummarizedExperiment::colData(object)[, colour_by])
      p <- ggplot2::ggplot(data, ggplot2::aes(x = .data$cell, y = .data$expr, fill = .data$group)) +
        ggplot2::scale_fill_brewer(type = "qual", palette = "Paired")
    } else {
      data <- data.frame(cell = factor(colnames(object), levels = colnames(object)[order]),
                         expr = ifelse(as.numeric(e) < threshold, 10^ - 2, as.numeric(e)))
      p <- ggplot2::ggplot(data, ggplot2::aes(x = .data$cell, y = .data$expr))
    }
    p <- p +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::ylab("Expression level") +
      ggplot2::xlab("Individual single-cells sorted by expression") +
      ggplot2::ggtitle(feature) +
      ggplot2::theme(axis.text.x  = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 8), panel.background = ggplot2::element_rect(fill = "#FFFFFF"))
    return(p)
  }
  else {
    warning(paste("Data did not contain a row", feature))
    return(NULL)
  }
}

#' Visualizes SOM components
#'
#' @param som The self-organizing map object.
#' @param codes The codes that should be displayed.
#' @param titles The title for each component.
#' @param what A keyword for what to plot.
#' @param aggregate.FUN A function to use for aggregation of codes.
#' @param scale Which scaling should be applied.
#' @param palette The colour palette used for scaling. If NULL, will be chosen automatically.
#' @param omit.title Whether the title should be omitted during plot.
#' @param omit.fill Whether the colour fill should be omitted during plot.
#' @param coords.fix Whether ratio of coordinates should be fixed to 1.
#'
#' @return A ggplot2 object or a list of ggplot2 objects, in case of multiple maps.
#'
#' @importFrom rlang .data
#' @export
visSOM <- function(som, codes=NULL, titles=NULL, what=c("code", "population", "aggregate", "overexpression", "underexpression"), aggregate.FUN = mean, scale = c("none", "minmax", "zscore"), palette = NULL, omit.title = FALSE, omit.fill = FALSE, coords.fix = TRUE) {
  if (!methods::is(som, "kohonen")) {
    stop("Supplied object som is not of class kohonen.")
  }
  if (is.null(codes)) {
    codes <- 1:ncol(som$codes)
  }
  if (!is.null(titles) && !is.null(codes) && length(titles) != length(codes)) {
    warning("Length of titels does not recapitulate length of codes.")
    titles <- rep(x = titles, length.out = length(som$codes))
  }

  scaling <- match.arg(scale)
  som$codes <- switch(scaling,
                      minmax = apply(som$codes, 2, function(x) 2 * (x + .Machine$double.eps - min(x)) / (max(x) - min(x)) - 1),
                      zscore = scale(som$codes),
                      som$codes)

  action <- match.arg(what)
  data <- switch(action,
                 code = as.matrix(som$codes[, codes]),
                 population = as.matrix(table(som$unit.classif)),
                 aggregate = as.matrix(apply(som$codes[, codes], 1, aggregate.FUN)),
                 overexpression = as.matrix(apply(apply(som$codes[, codes], 2, function(x) ifelse(x > stats::quantile(x, probs = 0.98), x, NA)), 1, mean, na.rm=T)),
                 underexpression = as.matrix(apply(apply(som$codes[, codes], 2, function(x) ifelse(x < stats::quantile(x, probs = 0.02), x, NA)), 1, mean, na.rm=T)))

  if (is.null(palette)) {
    palette <- switch(scaling,
                      minmax = viridis::magma(99),
                      zscore = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(10, "RdBu")))(99),
                      viridis::viridis(99))
  }
  num_maps <- ncol(data)
  grid <- as.data.frame(som$grid$pts)
  hexgrid <- hexcoords2(grid$x, grid$y, width = 1)

  res <- lapply(1:ncol(data), function(i) {
    df <- data.frame(hexgrid, value = rep(data[, i], each = 6))
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y, group = .data$group, fill = .data$value)) +
      ggplot2::geom_polygon()

    if (!omit.fill) p <- p + ggplot2::scale_fill_gradientn(colours = palette)
    if (coords.fix) p <- p + ggplot2::coord_fixed(ratio = 1)
    if (!omit.title) p <- p + ggplot2::ggtitle(titles[i])
    return(p)
  })
  if (length(res) == 1) {
    return(res[[1]])
  } else {
    return(res)
  }
}

#' Calculates corners from a hexagon center.
#'
#' @param x The hexagon centers x coordinate
#' @param y The hexagon centers y coordinate
#' @param size The size of the hexagon as height/2.
#'
#' @return A dataframe with x and y coordinates.
hex_corner <- function(x, y, size) {
  angle_deg <- 60 * 0:5 + 30
  angle_rad <- pi / 180 * angle_deg
  data.frame(x = x + size * cos(angle_rad),
             y = y + size * sin(angle_rad))
}

#' Calculates hexagon coordinates from hexagon centers.
#'
#' @param x The x coordinates of the centers
#' @param y The y coordinates of the centers
#' @param height The hexagons height. Precedes width if given.
#' @param width The hexagons width. Will be calculated from height by default.
#'
#' @return A dataframe with the hexagons x and y coordinates, as well a group variable indicating individual hexagons.
#'
#' @export
hexcoords2 <- function(x, y, height = NULL, width = sqrt(3) / 2 * height) {
  if (is.null(height) & is.null(width)) {
    stop("height and width can not be NULL at the same time.")
  }
  if (is.null(height)) height <- width * (2 / sqrt(3))
  size <- height / 2

  X <- grDevices::xy.coords(x, y)
  x <- cbind(X$x, X$y)
  res_l <- apply(x, 1, function(d, s) hex_corner(d[1], d[2], size = s), s = size)
  res <- do.call("rbind", res_l)
  res$group <- rep(1:nrow(x), each = 6)
  res
}

#' Draws a smoothScatter-Plot with hvg fit line and highlighted variable genes.
#'
#' @param data The normalized count table.
#' @param hvg.fit A list result from the hvg() function.
#' @param n Number of highly variable genes to highlight.
#'
#' @export
hvg.plot <- function(data, hvg.fit, n = 500) {
  means <- rowMeans(data)
  vars <- apply(data, 1, stats::var)
  cv2 <- vars / means^2

  xg <- exp(seq(min(log(means[means > 0])), max(log(means)), length.out = 1000))

  p <- lattice::xyplot(log(cv2) ~ log(means), xlab = "Average normalized read count", ylab = "Squared coefficient of variation",
              panel = function(x, y, ...) {
                lattice::panel.xyplot(x, y, ...)
                vfit <- hvg.fit$a1 / xg + hvg.fit$a0
                lattice::panel.lines(log(xg), log(vfit), col = 2)
                lattice::panel.lines(log(xg), log(vfit * stats::qchisq(0.975, hvg.fit$df) / hvg.fit$df), lty = 2, col = 2)
                lattice::panel.lines(log(xg), log(vfit * stats::qchisq(0.025, hvg.fit$df) / hvg.fit$df), lty = 2, col = 2)
                lattice::panel.points(x[hvg.fit$order[1:n]], y[hvg.fit$order[1:n]], col = 2)
              }, col = "black", hvg.fit = hvg.fit, xg = xg, n = n)
  return(p)
}

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
