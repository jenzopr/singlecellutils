#' Calculates a self-organizing map
#'
#' @param data A numerical matrix or data.frame.
#' @param sizemultiplier A size multiplier.
#' @param num_epochs How long the training should last.
#' @param train Row indices of training data.
#' @param seed A random seed (default NULL).
#' @param initrand Whether initialization should be random (not recommended).
#' @param intoroidal If a toroidal map should be created.
#'
#' @return A som object.
#'
#' @export
calcSOM <- function(data, sizemultiplier = 1, num_epochs = 200, train = NULL, seed = NULL, initrand = FALSE, intoroidal = FALSE) {
    if (num_epochs < length(train)/25) {
        lower_bound = floor(length(train)/25)
        warning(paste(num_epochs, "epochs is low for", length(train), "training genes, it was set to",
            lower_bound))
        num_epochs = lower_bound
    }

    if (!is.null(seed)) {
        set.seed(seed)
    }

    if (is.null(train)) {
        train = 1:min(nrow(datamat), 5000)
    }
    data.train = data[train, ]

    init.map = mapinit(data.train, sizemultiplier, intoroidal)
    # Handle output of mapinit

    maxr = min(0.5 * init.map$h, init.map$w)
    test.som <- kohonen::som(data = data.train, grid = class::somgrid(init.map$w, init.map$h, "hexagonal"), rlen = num_epochs,
        radius = c(maxr, 1), init = init.map$initgrid, toroidal = intoroidal)

    par(mfrow = c(2, 2))
    plot(test.som, type = "changes", main = "t")
    plot(test.som, type = "count")
    plot(test.som, type = "dist.neighbours")
    plot(test.som, main = "t")
    par(mfrow = c(1, 1))

    return(test.som)
}

#' Initializes a SOM map with two PCA components
#'
#' @param datamat A numerical matrix or data.frame.
#' @param prefactor A factor.
#' @param toroidal If toroidal.
#'
#' @return A list with intialization parameters.
mapinit <- function(datamat, prefactor = 1, toroidal = FALSE) {
    munits <- 5 * prefactor * sqrt(length(datamat[, 1]))
    pcar <- prcomp(datamat)
    ev <- pcar$sdev^2
    if (ev[1] > 5 * ev[2]) {
        h <- round(sqrt(munits * 5))
    } else {
        h <- round(sqrt(munits * ev[1]/ev[2]))
    }
    w <- round(munits/h)
    initcodes <- array(0, dim = c(h, w, length(datamat[1, ])))
    gridshape <- array(0, dim = c(h * w, length(datamat[1, ])))
    for (x in 1:w) {
        for (y in 1:h) {
            if (toroidal) {
                initcodes[y, x, ] = sin((pi * x)/w) * pcar$rotation[, 2] + sin((pi * y)/h) * pcar$rotation[,
                  1]
            } else {
                initcodes[y, x, ] = (x/w) * pcar$rotation[, 2] + (y/h) * pcar$rotation[, 1]
            }
            gridshape[(y - 1) * w + x, ] = initcodes[y, x, ]
        }
    }
    list(h = h, w = w, initgrid = gridshape)
}

#' Clusters a SOM into k regions determined by k-means
#'
#' @param som The SOM object.
#' @param kmax The number of clusters. If NULL, will determine k internally.
#' @param dist.fun Function to calculate neighborhood distances.
#' @param even.layout Whether the map layout is even (odd otherwise).
#' @param max.rounds Maximum number of cluster assignments.
#' @param ... Additional parameters.
#'
#' @return A clustering of the SOM.
#'
#' @export
clusterSOM <- function(som, kmax = NULL, dist.fun = median, even.layout = T, max.rounds = max(100,nrow(som$codes)/kmax), ...) {
  if(!is(som,"kohonen")) {
    stop("Supplied object som is not of class kohonen.")
  }
  codes <- som$codes
  # Get neighborhood
  neighborhood <- constructHexNeighborhood(som$grid$ydim, som$grid$xdim, even.layout = even.layout)
  # Create neighborhood distances
  n.distances <- sapply(1:nrow(codes), function(mi) {
    Ni <- neighborhood[[mi]]
    distances <- lapply(Ni, function(mj, mi) L2norm(codes[mi,]-codes[mj,]), mi=mi)
    do.call(dist.fun, list(x=unlist(distances)))
  })
  # Find local minima
  local.minima.set <- sapply(1:nrow(codes), function(mi) {
    Ni <- c(mi,neighborhood[[mi]])
    min <- which.min(n.distances[Ni])
    Ni[min]
  })
  # Maximum quorum
  local.minima.table <- table(local.minima.set)
  local.minima <- as.numeric(names(local.minima.table)[order(local.minima.table, decreasing = T)])
  # Removing neighboring local minima in decreasing order
  kept.minima <- rep(T,length(local.minima))
  for(i in 1:length(local.minima)) {
    if(kept.minima[i]) {
      minima.neighbors <- neighborhood[[local.minima[i]]]
      kept.minima[which(local.minima %in% minima.neighbors)] <- F
    }
  }
  # Remove non-selfvoters
  cluster.seeds <- local.minima[kept.minima]
  selfvoter <- (cluster.seeds == local.minima.set[cluster.seeds])
  cluster.seeds <- cluster.seeds[selfvoter]
  # Reduce cluster.seeds further if needed
  if(!is.null(kmax)) {
    cluster.seeds <- cluster.seeds[1:kmax]
  }
  # Prepare for clustering
  unassigned <- !(1:nrow(codes) %in% cluster.seeds)
  clustering <- as.list(cluster.seeds)
  # Clustering
  round <- 1
  while(sum(unassigned) > 0 & round < max.rounds) {
    cl.border.dists <- sapply(clustering, function(c) {
      cluster.seed <- c[1]
      cluster.neighborhood <- unique(unlist(neighborhood[c]))
      cluster.border <- cluster.neighborhood[which(unassigned[cluster.neighborhood])]
      if(length(cluster.border) == 0) { return(c(NA,Inf)) }
      border.dists <- sapply(cluster.border, function(mj,ci) L2norm(codes[ci,]-codes[mj,]), ci=cluster.seed)
      min <- which.min(border.dists)[1]
      return(c(cluster.border[min],border.dists[min]))
    })
    o <- order(cl.border.dists[2,])
    for(i in o) {
      new.member <- cl.border.dists[1,i]
      if(is.na(new.member)) { next; }
      if(unassigned[new.member]) {
        clustering[[i]] <- c(clustering[[i]], new.member)
        unassigned[new.member] <- F
      }
    }
    round <- round+1
  }
  # Prepare for return
  ret <- rep(NA, nrow(codes))
  for(i in 1:length(clustering)) {
    ret[clustering[[i]]] <- i
  }
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
  if(!is(som,"kohonen")) {
    stop("Supplied object som is not of class kohonen.")
  }
  if(is.null(clustering)) {
    clustering <- 1:nrow(som$codes)
  }
  k <- max(clustering)
  l <- lapply(1:k, function(i) {
    units <- which(clustering == i)
    ind <- which(som$unit.classif %in% units)
    if(!is.null(som$data)) {
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
  if(!is.logical(even.layout)) {
    stop("even.layout has to be logical!")
  }
  n <- r*q
  corners <- c(1,q,(n-q+1),n)

  nl <- lapply(1:n, function(i) {
    row <- ceiling(i/q)-1
    is.even.row <- ifelse(row %% 2 == 0, T, F)

    e <- i+1
    w <- i-1

    if(even.layout & is.even.row) {
      ne <- i+q+1
      nw <- i+q
      se <- i-q+1
      sw <- i-q
    }
    if(even.layout & !is.even.row) {
      ne <- i+q
      nw <- i+q-1
      se <- i-q
      sw <- i-q-1
    }
    if(!even.layout & is.even.row) {
      ne <- i+q
      nw <- i+q-1
      se <- i-q
      sw <- i-q-1
    }
    if(!even.layout & !is.even.row) {
      ne <- i+q+1
      nw <- i+q
      se <- i-q+1
      sw <- i-q
    }
    dir <- c(ne,e,se,sw,w,nw)
    dir[dir < 0] <- 0
    dir[dir > n] <- 0
    # Cell is on the right border
    if(i %% q == 0) {
      dir[dir %% q == 1] <- 0
    }
    else {
      dir[dir %% q == 0] <- 0
    }
    return(dir[dir!=0])
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
