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
    test.som <- kohonen::som(data = data.train, grid = kohonen::somgrid(init.map$w, init.map$h, "hexagonal"), rlen = num_epochs,
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
#' @param codes The codes matrix from a SOM object.
#' @param k Number of clusters. If NULL, will perform Gap statistics.
#' @param ... Additional parameters passed to hclust.
#'
#' @return A clustering of the SOM.
#'
#' @export
clusterSOM <- function(codes, k = NULL, ...) {
  if(is.null(k)) {
    cg <- cluster::clusGap(codes, kmeans, floor(sqrt(nrow(codes))), d.power = 2)
    k <- cluster::maxSE(cg$Tab[,3],cg$Tab[,4],method="firstSEmax",.25)
  }
  if(k < 2) {
    stop("Cannot cluster som codes into less than 2 clusters. Set k to at least 2.")
  }
  cutree(hclust(dist(codes), ...), k)
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
