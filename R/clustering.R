#' Calculates a (weighted) self-organizing map
#'
#' @param data A numerical matrix or data.frame.
#' @param train Row indices of training data.
#' @param weights Additonal weights per row
#' @param num_epochs How long the training should last.
#' @param resolution A multiplier for the map size
#' @param seed A random seed (default NULL).
#' @param init A list with initialization parameters.
#' @param parallel Indicates if the parallel implementation should be used.
#'
#' @return A som object.
#'
#' @export
calcSOM <- function(data, train = NULL, weights = NULL, num_epochs = 200, resolution = 1, seed = NULL, init = NULL, init.FUN = map.init, parallel = F) {
    if (is.null(train)) {
      train <- 1:min(nrow(data), 5000)
    }
    if (num_epochs < length(train) / 25) {
        lower_bound <- floor(length(train) / 25)
        warning(paste(num_epochs, "epochs is low for", length(train), "training genes, it was set to",
            lower_bound))
        num_epochs <- lower_bound
    }

    if (!is.null(seed)) {
        set.seed(seed)
    }

    data.train <- data[train, ]
    data.rest <- data[-train, ]

    if (is.null(weights)) {
      weights <- matrix(1L, nrow = nrow(data.train), ncol = ncol(data.train))
    }

    if (is.null(init)) {
      init <- do.call(init.FUN, list(data = data.train, resolution = resolution))
    }

    if (!parallel) {
      maxr <- min(0.5 * init$h, init$w)
      test.som <- kohonen::som(data = data.train, grid = class::somgrid(init$w, init$h, "hexagonal"), rlen = num_epochs,
                               radius = c(maxr, 1), init = init$initgrid, toroidal = F)
      par(mfrow = c(2, 2))
      plot(test.som, type = "changes", main = "t")
      plot(test.som, type = "count")
      plot(test.som, type = "dist.neighbours")
      plot(test.som, main = "t")
      par(mfrow = c(1, 1))
    }
    else {
      psom <- Rsomoclu::Rsomoclu.train(data.train, nEpoch = num_epochs, nSomX = init$w, nSomY = init$h, codebook = init$initgrid, radius0 = round(min(init$w, init$h)/2), radiusN = 1, radiusCooling = "linear", scale0 = 1, scaleN = 0.01, scaleCooling = "linear")
      test.som <- Rsomoclu::Rsomoclu.kohonen(data.train, psom, n.hood = "circular")
    }
    return(test.som)
}

#' Initializes a SOM map with the first two PCA components in a linear fashion
#'
#' @param data The training data of the SOM
#' @param resolution A size factor increasing or decreasing the map resolution
#'
#' @return A list with initialization parameters
map.init <- function(data, resolution = 1) {
  units <- 5 * resolution * sqrt(nrow(data))
  phi <- (1 + sqrt(5))/2
  h <- round(sqrt(units/phi))
  w <- round(h + (h/phi))

  pcar <- prcomp(data)
  initcodes <- array(0, dim = c(h, w, ncol(data)))
  gridshape <- array(0, dim = c(h * w, ncol(data)))
  for (x in 1:w) {
    for (y in 1:h) {
      initcodes[y, x, ] <- (x / w) * pcar$rotation[, 2] + (y / h) * pcar$rotation[,1]
      gridshape[(y - 1) * w + x, ] <- initcodes[y, x, ]
    }
  }
  list(h = h, w = w, initgrid = gridshape)
}

#' Initializes a SOM map with the first two PCA components in a localized fashion
#'
#' @param data The training data of the SOM
#' @param resolution A size factor increasing or decreasing the map resolution
#'
#' @return A list with initialization parameters
map.init.local <- function(data, resolution = 1) {
  units <- 5 * resolution * sqrt(nrow(data))
  phi <- (1 + sqrt(5))/2
  h <- round(sqrt(units/phi))
  w <- round(h + (h/phi))

  pcar <- prcomp(data)
  gridshape <- array(0, dim = c(h * w, ncol(data)))

  h_pos <- as.numeric(cut(pcar$x[,1], breaks=h))
  w_pos <- as.numeric(cut(pcar$x[,2], breaks=w))
  cell_pos <- (h_pos - 1) * w + w_pos
  n <- constructHexNeighborhood(h, w, even.layout = T)

  h_genes <- rownames(datamat)[order(abs(pcar$rotation[,1]), decreasing = T)[1:100]]
  w_genes <- rownames(datamat)[order(abs(pcar$rotation[,2]), decreasing = T)[1:100]]
  genes <- unique(c(h_genes, w_genes))

  for(i in 1:ncol(data)) {
    direct <- c(cell_pos[i], n[[cell_pos[i]]])
    second <- unique(unlist(sapply(direct, function(k) { n[[k]] })))
    third <- unique(unlist(sapply(second, function(k) { n[[k]] })))
    gridshape[third,i] <- sum(1/3 * data[genes, i])
    gridshape[second,i] <- gridshape[second,i] + sum(1/2 * data[genes, i])
    gridshape[direct,i] <- gridshape[direct,i] + sum(data[genes, i])
  }
  list(h = h, w = w, initgrid = scale(gridshape))
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
#'
#' @export
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

#' A wrapper for running HDBSCAN in python
#'
#' @param data Matrix-like data to be clustered
#' @param min.cluster.size The minimum size of clusters
#' @param min.samples The number of samples in a neighbourhood for a point to be considered a core point.
#' @param kwargs A named list of parameters passed to HDBSCAN
#'
#' @return the clustered data
#'
#' @export
hdbscan <- function(data, min.cluster.size = 5, min.samples = NULL, kwargs = NULL, return.data = T) {
  data <- as.matrix(data)

  if (nrow(data) < min.cluster.size) {
    stop("Less observations than min.cluster.size")
  }
  if ((!is.null(kwargs) && !is.list(kwargs)) || (is.list(kwargs) && is.null(names(kwargs)))) {
    stop("kwargs should be a named list")
  }

  rPython::python.exec("import hdbscan")
  rPython::python.exec("import numpy as np")

  # Assign variables
  rPython::python.assign("l", data)
  rPython::python.exec("data = np.array(l)")
  rPython::python.exec("min_samples = None")
  rPython::python.exec("additional_arguments = dict()")
  rPython::python.assign("min_cluster_size", min.cluster.size)

  if (!is.null(min.samples)) {
    rPython::python.assign("min_samples", min.samples)
  }
  if (!is.null(kwargs)) {
    rPython::python.assign("additional_arguments", kwargs)
  }

  # Initialise clusterer
  rPython::python.exec("cl = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size, min_samples=min_samples, **additional_arguments)")

  # Cluster data
  rPython::python.exec("labels = cl.fit_predict(data).tolist()")
  rPython::python.exec("prob = cl.probabilities_.tolist()")
  rPython::python.exec("persistence = cl.cluster_persistence_.tolist()")
  #rPython::python.exec("outlier_score = cl.outlier_scores_.tolist()")

  # Get results
  .cluster <- rPython::python.get("labels")
  .probabilities <- rPython::python.get("prob")
  .persistence <- rPython::python.get("persistence")
  #.outlier.score <- rPython::python.get("outlier_score")

  .cluster <- .cluster + 1
  if (return.data) {
    data.frame(data, .cluster = .cluster, .probabilities = .probabilities)
  } else {
    data.frame(.cluster = .cluster, .probabilities = .probabilities)
  }
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

