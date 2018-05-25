#' Adds results from a particular clustering technique to the object.
#'
#' @param object A SingleCellExperiment object.
#' @param flavor Determines which clustering technique to apply.
#' @param column Determines the column name of the \code{colData} slot to store results to.
#' @param ... Additional parameters passed to functions.
#'
#' @return A SingleCellExperiment object with modified \code{colData} slot.
#'
#' @export
add_clustering <- function(object, flavor = c("hdbscan"), column = ".cluster", ...) {
  clustering <- switch(flavor,
                       hdbscan = hdbscan(object, ...))
  object <- SingleCellExperiment::mutate(object, !!column := clustering)
  return(object)
}

#' Conveinience function to perform HDBSCAN via reticulate.
#'
#' Performs HDBSCAN clustering on a SingleCellExperiment object.
#' The algorithm is explained in detail in http://joss.theoj.org/papers/10.21105/joss.00205 and needs the python package \code{hdbscan} installed.
#'
#' @param object A SingleCellExperiment object.
#' @param use_dimred A string or integer scalar indicating the reduced dimension result in \code{reducedDims(object)} to use as input.
#' @param min_samples Measure of how conservative the clustering should to be. The larger the value of \code{min_samples}, the more conservative the clustering – more points will be declared as noise, and clusters will be restricted to progressively more dense areas.
#' @param min_cluster_size The smallest size grouping that is considered a cluster.
#' @param outlier Determines how outliers are encoded in the resulting clustering.
#' @param seed A numeric seed to initialize the random number generator.
#'
#' @return A factor with the assigned cluster.
#'
#' @export
hdbscan <- function(object, use_dimred, min_samples = 7L, min_cluster_size = 9L, outlier = 0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  h <- reticulate::import("hdbscan")
  cl <- h$HDBSCAN(min_samples = min_samples, min_cluster_size = min_cluster_size)

  labels <- cl$fit_predict(SingleCellExperiment::reducedDim(object, use_dimred)) + 1
  labels[labels == 0] <- outlier
  return(factor(labels))
}

#' Calculates a summarized som for a list of components
#'
#' @param som The SOM object
#' @param codes The components that should be sumarized
#' @param summarize.fun The function that is used to summarize each SOM node
#' @param ... Additional parameters passed to summarize.fun
#'
#' If \code{codes} is a list of component numbers, the summarization is done sequentially for each list.
#'
#' @return A summarized SOM codes matrix. Each column corresponds to a list element from the \code{codes} input.
summarizeSOM <- function(som, codes = NULL, summarize.fun = mean, ...) {
  if (!is(som, "kohonen")) {
    stop("Supplied object som is not of class kohonen.")
  }
  if (is.null(codes)) {
    codes <- 1:ncol(som$codes)
  }
  if (typeof(codes) != "list") {
    codes <- list(codes)
  }

  dat <- lapply(codes, function(c) {
    apply(som$codes[, c], 1, summarize.fun, ...)
  })
  ret <- do.call("cbind", dat)
  return(ret)
}

#' Extracts the rows of the data matrix underlying each SOM node.
#'
#' @param som The SOM object.
#' @param clustering An optional clustering of SOM nodes.
#'
#' @return Row indices or data row names for each node or each cluster as a list.
getSOMgenes <- function(som, clustering = NULL) {
  if (!is(som, "kohonen")) {
    stop("Supplied object som is not of class kohonen.")
  }
  if (is.null(clustering)) {
    clustering <- 1:nrow(som$codes)
  }
  k <- max(clustering)
  l <- lapply(1:k, function(i) {
    units <- which(clustering == i)
    ind <- which(som$unit.classif %in% units)
    if (!is.null(som$data)) {
      return(rownames(som$data)[ind])
    }
    else {
      return(ind)
    }
  })
  names(l) <- 1:k
  return(l)
}

#' Construct the neighborhood of a hexagonal r*q matrix
#'
#' @param r The number of rows
#' @param q The number of columns
#' @param even.layout Whether the layout is even
#'
#' @return A list of length n*q with their neighbors on the grid.
constructHexNeighborhood <- function(r, q, even.layout = T) {
  if (!is.logical(even.layout)) {
    stop("even.layout has to be logical!")
  }
  n <- r * q


  nl <- lapply(1:n, function(i) {
    row <- ceiling(i / q) - 1
    is_even_row <- ifelse(row %% 2 == 0, T, F)

    e <- i + 1
    w <- i - 1

    if (even.layout & is_even_row) {
      ne <- i + q + 1
      nw <- i + q
      se <- i - q + 1
      sw <- i - q
    }
    if (even.layout & !is_even_row) {
      ne <- i + q
      nw <- i + q - 1
      se <- i - q
      sw <- i - q - 1
    }
    if (!even.layout & is_even_row) {
      ne <- i + q
      nw <- i + q - 1
      se <- i - q
      sw <- i - q - 1
    }
    if (!even.layout & !is_even_row) {
      ne <- i + q + 1
      nw <- i + q
      se <- i - q + 1
      sw <- i - q
    }
    dir <- c(ne, e, se, sw, w, nw)
    dir[dir < 0] <- 0
    dir[dir > n] <- 0
    # Cell is on the right border
    if (i %% q == 0) {
      dir[dir %% q == 1] <- 0
    }
    else {
      dir[dir %% q == 0] <- 0
    }
    return(dir[dir != 0])
  })
  return(nl)
}

#' Calculates the L2 norm of a numeric vector
#'
#' @param x A numeric vector.
#'
#' @return The L2 norm of x.
L2norm <- function(x) {
  sqrt(sum(x^2))
}

#' Calculates cluster centers based on a probability threshold.
#'
#' @param c The clustering
#' @param min.probability The minimum probability
#'
#' @return The cluster centers
calculateClusterCenter <- function(c, min.probability = 0.4) {
  clusters <- unique(c$.cluster[c$.cluster > 0])
  data_columns <- which(colnames(c) == ".cluster") - 1
  centers <- sapply(clusters, function(x) {
    i <- c$.cluster == x & c$.probabilities > min.probability
    if (data_columns > 1) {
      colMeans(c[i,1:data_columns])
    } else {
      mean(c[i,1:data_columns])
    }
  })
  centers
}
