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
  if(!is.null(medoids)) {
    k <- length(medoids)
  }
  if(is.null(k)) {
    cg <- cluster::clusGap(tsne, cluster::pam, 10, d.power = 2)
    k <- cluster::maxSE(cg$Tab[,3],cg$Tab[,4],method="firstSEmax",.25)
  }
  if(k < 2) {
    stop("Cannot cluster tSNE into less than 2 clusters. Set k to at least 2.")
  }
  cl <- cluster::pam(x=tsne, k = k, medoids = medoids, ...)
  pal <- RColorBrewer::brewer.pal(k, use.pal)

  plot(tsne[,1],tsne[,2],xlab = "tSNE dimension 1", ylab = "tSNE dimension 2", col = pal[cl$clustering], pch = 20)
  if(add.center) {
    points(cl$medoids[,1],cl$medoids[,2], pch = 25, col = pal)
    text(cl$medoids[,1]+0.05,cl$medoids[,2], labels = 1:k)
  }

  return(invisible(broom::augment(cl, tsne)))
}

#' A function to color a tSNE map based on a vector.
#'
#' @param tsne The data matrix from a tSNE dimension reduction call
#' @param data A numerical vector on which the coloring is based
#' @param data.name The name of the data, shown in plot.main
#' @param use.pal The RColorBrewer color palette
#'
#' @return NULL
#'
#' @export
colortSNE <- function(tsne, data, data.name = NULL, use.pal = "RdBu") {
  if(length(data) != nrow(tsne)) {
    stop("Data is not same length as tsne.")
  }
  if(!is.null(data.name)) {
    main = data.name
  } else {
    main = ""
  }
  pal <- colorRampPalette(RColorBrewer::brewer.pal(10,use.pal))
  color <- rev(pal(100))[as.numeric(cut(data,breaks=100))]
  plot(tsne[,1],tsne[,2],xlab = "tSNE dimension 1", ylab = "tSNE dimension 2", main = main, col = color, pch = 20)
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

#' Visualizes SOM components
#'
#' @param som The self-organizing map object.
#' @param code The codes that should be displayed.
#' @param titles The title for each component.
#'
#' @return The matrix of color codes (invisible)
#'
#' @export
visSOM <- function(som, code=NULL, titles=NULL) {
  if(!is(som,"kohonen")) {
    stop("Supplied object som is not of class kohonen.")
  }
  if(is.null(code)) {
    code <- 1:ncol(som$codes)
  }
  if(!is.null(titles) && length(titles) != length(code)) {
    warning("Length of titels does not recapitulate length of codes.")
    titles <- rep(x = titles, length.out = length(codes))
  }
  data <- scale(as.matrix(som$codes[,code]))
  rows <- som$grid$ydim
  cols <- som$grid$xdim

  par(mar=c(1,1,3,1))
  palette <- rev(colorRampPalette(RColorBrewer::brewer.pal(9,"Spectral"))(50))

  colors <- apply(data, 2, function(x) {
    c <- cut(x, breaks=50, labels=FALSE)
    ifelse(is.na(c),"#FFFFFF",palette[c])
  })

  if(length(code) > 1) {
    par(mfrow=c(3,3))
  }
  shift <- 0.5
  for(i in 1:length(code)) {
    plot(0, 0, type = "n", axes = FALSE, xlim=c(0, cols), ylim=c(0, rows), xlab="", ylab= "", asp=1, main=titles[i])
    for(row in 0:(rows-1)) {
      for(col in 1:cols)
        Hexagon(col + shift, row, col = colors[row*cols+col,i])
      shift <- ifelse(shift, 0, 0.5)
    }
  }

  par(mfrow=c(1,1))
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
  polygon(c(x, x, x + unitcell/2, x + unitcell, x + unitcell,
            x + unitcell/2), c(y + unitcell * 0.125,
                               y + unitcell * 0.875,
                               y + unitcell * 1.125,
                               y + unitcell * 0.875,
                               y + unitcell * 0.125,
                               y - unitcell * 0.125),
          col = col, border=NA)
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
  if(is.null(gene.names)) {
    gene.names <- row.names(data)
  }
  if(is.null(group)) {
    group <- rep(1, ncol(data))
  }
  data$gene <- factor(gene.names, levels = gene.names[order(gene.names, decreasing = T)])
  m <- reshape2::melt(data, id.vars=c("gene"))
  colnames(m) <- c("gene","cell","value")
  m$group <- gsub(group,"\\1",m$cell)

  p = ggplot2::ggplot(m, ggplot2::aes(gene, value)) +
    ggplot2::geom_violin(trim=T,scale="width") +
    ggplot2::geom_jitter(alpha=0.5, ggplot2::aes(color=group), width = 0.75) +
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
  cv2 <- vars/means^2

  smoothScatter(log(means),log(cv2), xlab = "Average normalized read count", ylab = "Squared coefficient of variation")
  xg <- exp(seq(min(log(means[means>0])), max(log(means)), length.out=1000))
  vfit <- hvg.fit$a1/xg + hvg.fit$a0
  lines(log(xg),log(vfit), col="black", lwd=3)
  lines(log(xg),log(vfit * qchisq(0.975,hvg.fit$df)/hvg.fit$df),lty=2,col="black")
  lines(log(xg),log(vfit * qchisq(0.025,hvg.fit$df)/hvg.fit$df),lty=2,col="black")

  points(log(means[hvg.fit$order[1:n]]),log(cv2[hvg.fit$order[1:n]]),col=2)
}
