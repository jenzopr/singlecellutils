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

#' Clusters a SOM into k regions determined by k-means
#'
#' @param som The SOM object.
#' @param kmax The number of clusters. If NULL, will determine k internally.
#' @param dist.fun Function to aggregate neighborhood distances.
#' @param use.codes Perform clustering only on a subset of codes.
#' @param even.layout Whether the map layout is even (odd otherwise).
#' @param max.rounds Maximum number of cluster assignments.
#' @param ... Additional parameters.
#'
#' @return A clustering of the SOM.
clusterSOM <- function(som, kmax = NULL, dist.fun = median, use.codes = NULL, even.layout = T, max.rounds = max(100, nrow(som$codes)), ...) {
  if (!is(som, "kohonen")) {
    stop("Supplied object som is not of class kohonen.")
  }
  if (is.null(use.codes)) {
    use.codes <- 1:ncol(som$codes)
  }
  codes <- som$codes[, use.codes]
  # Get neighborhood
  neighborhood <- constructHexNeighborhood(som$grid$ydim, som$grid$xdim, even.layout = even.layout)
  # Create neighborhood distances
  n.distances <- sapply(1:nrow(codes), function(mi) {
    Ni <- neighborhood[[mi]]
    #distances <- lapply(Ni, function(mj, mi) L2norm(codes[mi,]-codes[mj,]), mi=mi)
    distances <- lapply(Ni, function(mj, mi) 1 - cor(codes[mi, ], codes[mj, ]), mi = mi)
    do.call(dist.fun, list(x = unlist(distances)))
  })
  # Find local minima
  local_minima_set <- sapply(1:nrow(codes), function(mi) {
    Ni <- c(mi, neighborhood[[mi]])
    min <- which.min(n.distances[Ni])
    Ni[min]
  })
  # Maximum quorum
  local_minima_table <- table(local_minima_set)
  local_minima <- as.numeric(names(local_minima_table)[order(local_minima_table, decreasing = T)])
  # Removing neighboring local minima in decreasing order
  kept_minima <- rep(T, length(local_minima))
  for (i in 1:length(local_minima)) {
    if (kept_minima[i]) {
      minima_neighbors <- neighborhood[[local_minima[i]]]
      kept_minima[which(local_minima %in% minima_neighbors)] <- F
    }
  }
  # Remove non-selfvoters
  cluster_seeds <- local_minima[kept_minima]
  selfvoter <- (cluster_seeds == local_minima_set[cluster_seeds])
  cluster_seeds <- cluster_seeds[selfvoter]
  # Reduce cluster.seeds further if needed
  if (!is.null(kmax)) {
    cluster_seeds <- cluster_seeds[1:min(kmax, length(cluster_seeds))]
  }
  # Prepare for clustering
  unassigned <- !(1:nrow(codes) %in% cluster_seeds)
  clustering <- as.list(cluster_seeds)
  # Clustering
  round <- 1
  while(sum(unassigned) > 0 & round < max.rounds) {
    cl_border_dists <- sapply(clustering, function(c) {
      cluster_seed <- c[1]
      cluster_neighborhood <- unique(unlist(neighborhood[c]))
      cluster_border <- cluster_neighborhood[which(unassigned[cluster_neighborhood])]
      if (length(cluster_border) == 0) {
        return(c(NA, Inf))
      }
      #border.dists <- sapply(cluster.border, function(mj,ci) L2norm(codes[ci,]-codes[mj,]), ci=cluster.seed)
      border_dists <- sapply(cluster_border, function(mj, ci) 1 - cor(codes[ci, ], codes[mj, ]), ci = cluster_seed)
      min <- which.min(border_dists)[1]
      return(c(cluster_border[min], border_dists[min]))
    })
    i <- which.min(cl_border_dists[2, ])
    new_member <- cl_border_dists[1, i]
    if (is.na(new.member)) {
      next
    }
    if (unassigned[new_member]) {
      clustering[[i]] <- c(clustering[[i]], new_member)
      unassigned[new_member] <- F
    }
    round <- round + 1
  }
  # Prepare for return
  ret <- rep(NA, nrow(codes))
  for (i in 1:length(clustering)) {
    ret[clustering[[i]]] <- i
  }
  return(ret)
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
#'
#' @export
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
#'
#' @export
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
