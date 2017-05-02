#' Visualizes tSNE maps with colored clusters by pam.
#'
#' @param tsne The data matrix from a tSNE dimension reduction call
#' @param k The number of clusters. If NULL, will perform Gap statistics.
#' @param use.pal The RColorBrewer palette name to use for cluster colors.
#' @param add.center Whether cluster centers should be added to the plot.
#' @param medoids length-k vector of integer indices specifying intial medoids.
#' @param ... Other parameters passed to cluster::pam.
#'
#' @return The augmented input data matrix. Clusters are added in .cluster column.
#'
#' @export
vistSNE <- function(tsne, k = NULL, use.pal = "Dark2", add.center = T, medoids = NULL, ...) {
  if (!is.null(medoids)) {
    k <- length(medoids)
  }
  if (is.null(k)) {
    cg <- cluster::clusGap(tsne, cluster::pam, 10, d.power = 2)
    k <- cluster::maxSE(cg$Tab[, 3], cg$Tab[, 4], method = "firstSEmax", .25)
  }
  if (k < 2) {
    stop("Cannot cluster tSNE into less than 2 clusters. Set k to at least 2.")
  }
  cl <- cluster::pam(x = tsne, k = k, medoids = medoids, ...)
  pal <- RColorBrewer::brewer.pal(k, use.pal)

  plot(tsne[, 1], tsne[, 2], xlab = "tSNE dimension 1", ylab = "tSNE dimension 2", col = pal[cl$clustering], pch = 20)
  if (add.center) {
    points(cl$medoids[, 1], cl$medoids[, 2], pch = 25, col = pal)
    text(cl$medoids[, 1] + 0.05, cl$medoids[, 2], labels = 1:k)
  }

  return(invisible(broom::augment(cl, tsne)))
}

#' A function to color observations on a map based on values of a vector.
#'
#' @param data The data matrix with observation coordinates in the first two columns.
#' @param values A numerical vector on which the coloring is based, or a numeric or character of length 1 indicating the column of the data matrix from which to take the values.
#' @param palette The color palette
#' @param title The name of the data, shown in plot.main
#' @param outlier.col The color of outlier points (\code{NA} in \code{values})
#' @param ... Other parameters passed to lattice::xyplot
#'
#' @return A lattice xyplot.
#'
#' @export
colorAMap <- function(data, values, palette = RColorBrewer::brewer.pal(8, "Dark2"), title = NULL, outlier.col = "darkgrey", xlab = "Dimension 1", ylab = "Dimension 2", ...) {
  if (length(values) != nrow(data) & length(values) > 1) {
    stop("Values are not of same length as data and more than one value given")
  }
  if (!is.null(title)) {
    main <- title
  } else {
    main <- ""
  }

  if (length(values) == 1 & ncol(data) <= length(values)) {
    values <- as.numeric(data[, values])
  }

  if (length(unique(values)) <= length(palette)) {
    warning("Less unique data values than colors in palette.. reducing palette.")
  }
  color <- palette[as.numeric(cut(values, breaks = length(palette)))]

  if (any(is.na(values))) {
    color[is.na(values)] = outlier.col
  }

  p <- lattice::xyplot(data[, 2] ~ data[, 1], xlab = xlab, ylab = ylab, main = main, col = color, ...)
  return(p)
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
  if (!is.null(cell.names) & ncol(data) == length(cell.names)) {
    colnames(data) <- cell.names
  }
  if (!is.null(y.unit)) {
    y.label <- paste("(", y.unit, ")", sep = "")
  } else {
    y.label <- ""
  }
  if (is.null(cell.groups)) {
    cell.groups <- rep("none", ncol(data))
  }
  e <- data[marker, ]
  if (!is.null(e)) {
    order <- order(e, decreasing = decreasing)
    data <- data.frame(cell = factor(colnames(e), levels = colnames(e)[order]),
                       expr = ifelse(as.numeric(e) < threshold, 10^ - 2, as.numeric(e)),
                       group = cell.groups)
    p <- ggplot2::ggplot(data, aes(x = cell, y = expr, fill = group)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::scale_fill_brewer(type = "qual", palette = "Paired") +
      ggplot2::ylab(paste("Expression level", y.label)) +
      ggplot2::xlab("Individual single-cells sorted by expression") +
      ggplot2::ggtitle(marker) +
      ggplot2::theme(axis.text.x  = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 8), panel.background = element_rect(fill = "#FFFFFF"))
    return(p)
  }
  else {
    warning(paste("Data did not contain a row", marker))
    return(NULL)
  }
}

#' Visualizes SOM components
#'
#' @param som The self-organizing map object.
#' @param code The codes that should be displayed.
#' @param titles The title for each component.
#' @param what A keyword for what to plot.
#' @param aggregate.FUN A function to use for aggregation of codes.
#' @param scale Which scaling should be applied.
#'
#' @return The matrix of color codes (invisible)
#'
#' @export
visSOM <- function(som, codes=NULL, titles=NULL, what=c("code", "population", "aggregate", "overexpression", "underexpression"), aggregate.FUN = mean, scale = c("none", "minmax", "zscore")) {
  if (!is(som, "kohonen")) {
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
                 overexpression = as.matrix(apply(apply(som$codes[, codes], 2, function(x) ifelse(x > quantile(x, probs = 0.98), x, NA)), 1, mean, na.rm=T)),
                 underexpression = as.matrix(apply(apply(som$codes[, codes], 2, function(x) ifelse(x < quantile(x, probs = 0.02), x, NA)), 1, mean, na.rm=T)))
  num_maps <- ncol(data)
  rows <- som$grid$ydim
  cols <- som$grid$xdim

  par(mar = c(1, 1, 3, 1))
  palette <- rev(colorRampPalette(RColorBrewer::brewer.pal(9, "Spectral"))(99))

  cutpoints <- seq(min(som$codes[, codes], na.rm=T), max(som$codes[, codes], na.rm=T), length.out = 99)
  colors <- apply(data, 2, function(x) {
    c <- cut(x, breaks = cutpoints, labels = FALSE)
    ifelse(is.na(c), "#FFFFFF", palette[c])
  })

  if (num_maps > 1) {
    par(mfrow = c(3, 3))
  }
  for (i in 1:num_maps) {
    shift <- 0.5
    plot(0, 0, type = "n", axes = FALSE, xlim = c(0, cols), ylim = c(0, rows), xlab = "", ylab = "", asp = 1, main = titles[i])
    for (row in 0:(rows - 1)) {
      for (col in 1:cols)
        Hexagon(col + shift, row, col = colors[row * cols + col, i])
      shift <- ifelse(shift, 0, 0.5)
    }
  }

  par(mfrow = c(1, 1))
  return(invisible(colors))
}

#' Creates a polygon
#'
#' @param x The polygons x
#' @param y The polygons y
#' @param unitcell Scaling
#' @param col The color of the polygon
#'
#' @return A colored polygon
Hexagon <- function (x, y, unitcell = 1, col = col) {
  polygon(c(x, x, x + unitcell / 2, x + unitcell, x + unitcell,
            x + unitcell / 2), c(y + unitcell * 0.125,
                               y + unitcell * 0.875,
                               y + unitcell * 1.125,
                               y + unitcell * 0.875,
                               y + unitcell * 0.125,
                               y - unitcell * 0.125),
          col = col, border = NA)
}

#' Creates a violin plot for selected genes
#'
#' @param data The expression data (as data.frame)
#' @param gene.names A character vector of associated gene names (optional)
#'
#' @return A ggplot2 figure.
#'
#' @export
visMarkerViolin <- function(data, gene.names = NULL, group = NULL) {
  if (is.null(gene.names)) {
    gene.names <- row.names(data)
  }
  if (is.null(group)) {
    group <- rep(1, ncol(data))
  }
  data$gene <- factor(gene.names, levels = gene.names[order(gene.names, decreasing = T)])
  m <- reshape2::melt(data, id.vars = c("gene"))
  colnames(m) <- c("gene", "cell", "value")
  m$group <- gsub(group, "\\1", m$cell)

  p <- ggplot2::ggplot(m, ggplot2::aes(gene, value)) +
    ggplot2::geom_violin(trim = T, scale = "width") +
    ggplot2::geom_jitter(alpha = 0.5, ggplot2::aes(color = group), width = 0.75) +
    ggplot2::coord_flip() +
    ggplot2::xlab("") +
    ggplot2::ylab("log2( Transcripts per Million )") #+
    #ggplot2::scale_colour_manual(values = cols)
  return(p)
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
  vars <- apply(data, 1, var)
  cv2 <- vars / means^2

  xg <- exp(seq(min(log(means[means > 0])), max(log(means)), length.out = 1000))

  p <- lattice::xyplot(log(cv2) ~ log(means), xlab = "Average normalized read count", ylab = "Squared coefficient of variation",
              panel = function(x, y, ...) {
                lattice::panel.xyplot(x, y, ...)
                vfit <- hvg.fit$a1 / xg + hvg.fit$a0
                lattice::panel.lines(log(xg), log(vfit), col = 2)
                lattice::panel.lines(log(xg), log(vfit * qchisq(0.975, hvg.fit$df) / hvg.fit$df), lty = 2, col = 2)
                lattice::panel.lines(log(xg), log(vfit * qchisq(0.025, hvg.fit$df) / hvg.fit$df), lty = 2, col = 2)
                lattice::panel.points(x[hvg.fit$order[1:n]], y[hvg.fit$order[1:n]], col = 2)
              }, col = "black", hvg.fit = hvg.fit, xg = xg, n = n)
  return(p)
}

#' Plot cells along their fraction of identity to a given state
#'
#' @param identity A vector giving the cells identity value or a list from calcFractionOfIdentity call
#' @param data A numerical vector according which the cells are colored
#' @param data.name A name that will be shown as plot.main
#' @param data.discrete Logical indicating whether data is discrete (like groups)
#' @param pal The color palette
#' @param state If identity is a list, chooses the given state for ordering
#' @param decreasing Logical, whether ordering should be decreasing (the default)
#'
#' @export
colorIdentity <- function(identity, data=NULL, data.name=NULL, data.discrete=FALSE, pal=viridis::viridis(99), state=NULL, decreasing = TRUE) {
  if (typeof(identity) != "double" & typeof(identity) != "list") {
    stop("Parameter identity should be of type double or list.")
  }
  if ( typeof(identity) == "double") {
    id <- as.matrix(identity)
    state <- 1
  } else {
    if ( is.null(state)) {
      warning("When identity is a list, state cannot be NULL. Taking state 1 instead..")
      state <- 1
    }
    id <- identity$identity
  }
  ind <- 1:nrow(id)
  order <- order(id[, state], decreasing = decreasing)

  colors <- "black"
  if ( !is.null(data) ) {
    colors <- pal[findInterval(data, seq(0, 1, length.out = length(pal) + 1), all.inside = TRUE)][order]
    if ( data.discrete ) {
      colors <- pal[as.numeric(data)][order]
    }
  }

  plot(ind, id[order, state], col = colors, ylim = c(0, 1), pch = 16, bty = "n", axes = F, xlab = "pseudotime", ylab = "fraction of identity", main = data.name)
  axis(1, labels = F)
  axis(2, labels = T)
}
