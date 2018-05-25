#' Adds results from a particular dimension reduction technique to the object.
#'
#' @param object A SingleCellExperiment object.
#' @param flavor Determines which dimension reduction technique to apply.
#' @param slot Determines which entry of the \code{reducedDims} slot to use for reduced embedding.
#' @param ... Additional parameters passed to functions.
#'
#' @return A SingleCellExperiment object with modified \code{reducedDims} slot.
#'
#' @export
reduce_dimension <- function(object, flavor = c("umap", "som"), slot = flavor, ...) {
  embedding <- switch(flavor,
    umap = umap(object, ...),
    som = t(som(object, ...)$codes[[1]]))
  SingleCellExperiment::reducedDim(object, slot) <- embedding
  return(object)
}

#' Conveinience function to perform UMAP via reticulate.
#'
#' Performs Uniform Manifold Approximation and Projection (UMAP) on a SingleCellExperiment object.
#' The algorithm is explained in detail in https://arxiv.org/abs/1802.03426 and needs the python packages umap and numpy installed.
#'
#' @param object A SingleCellExperiment object.
#' @param exprs_values String indicating which assay contains the data that should be used to perform UMAP.
#' @param features A character vector (of feature names), a logical vector or numeric vector (of indices) specifying the features to use for UMAP. The default of NULL will use all features.
#' @param scale Logical scalar indicating whether data should be scaled before UMAP.
#' @param n_neighbors This determines the number of neighboring points used in local approximations of manifold structure. Larger values will result in more global structure being preserved at the loss of detailed local structure. In general this parameter should often be in the range 5 to 50, with a choice of 10 to 15 being a sensible default.
#' @param min_dist This controls how tightly the embedding is allowed compress points together. Larger values ensure embedded points are more evenly distributed, while smaller values allow the algorithm to optimise more accurately with regard to local structure. Sensible values are in the range 0.001 to 0.5, with 0.1 being a reasonable default.
#' @param metric This determines the choice of metric used to measure distance in the input space.
#' @param seed A numeric seed to initialize the random number generator.
#'
#' @return A matrix with the two-dimensional embedding.
#'
#' @export
umap <- function(object, exprs_values = 'norm_TPM', features = NULL, scale = T, n_neighbors = 5L, min_dist = 0.1, metric = "correlation", seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  umap <- reticulate::import("umap")
  np <- reticulate::import("numpy", convert=FALSE)
  umap_cl <- umap$UMAP(n_neighbors = as.integer(n_neighbors), min_dist = min_dist, metric = metric)

  input <- assay(object, i = exprs_values)
  if (!is.null(features)) {
    input <- SummarizedExperiment::assay(object, i = exprs_values)[features,]
  }
  if (scale) input <- t(scale(t(input)))

  embedding <- umap_cl$fit_transform(t(input))
  return(embedding)
}

#' Calculates a (weighted) self-organizing map
#'
#' @param object A SingleCellExperiment object.
#' @param exprs_values String indicating which assay contains the data that should be used to perform SOM calculation.
#' @param features A character vector (of feature names), a logical vector or numeric vector (of indices) specifying the features to use for SOM calculation. The default of NULL will use all features.
#' @param scale Logical scalar indicating whether data should be scaled before SOM calculation.
#' @param num_epochs How long the training should last.
#' @param resolution A multiplier for the map size.
#' @param seed A numeric seed to initialize the random number generator.
#' @param init A list with initialization parameters.
#' @param init.FUN A function to initialize the map.
#' @param ... Additional arguments passed to \code{kohonen::som}.
#'
#' @return A som object.
#'
#' @export
calcSOM <- function(object, exprs_values = 'norm_TPM', features = NULL, scale = T, num_epochs = 200, resolution = 1, seed = NULL, init = NULL, init.FUN = singlecellutils::map.init, ...) {
  if (!is.null(features)) {
    input <- SummarizedExperiment::assay(object, i = exprs_values)[features, ]
  } else {
    input <- assay(object, i = exprs_values)[1:min(nrow(data), 5000), ]
  }
  if (scale) input <- t(scale(t(input)))

  if (num_epochs < length(train) / 25) {
    lower_bound <- floor(length(train) / 25)
    warning(paste(num_epochs, "epochs is low for", length(train), "training genes, it was set to",
                  lower_bound))
    num_epochs <- lower_bound
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (is.null(init)) {
    init <- do.call(init.FUN, list(data = input, resolution = resolution))
  }

  maxr <- min(0.5 * init$h, init$w)
  som <- kohonen::som(X = input, grid = kohonen::somgrid(init$w, init$h, "hexagonal", toroidal = F), rlen = num_epochs,
                           radius = c(maxr, 1), init = init$initgrid, ...)
  # graphics::par(mfrow = c(2, 2))
  # plot(test.som, type = "changes", main = "t")
  # plot(test.som, type = "count")
  # plot(test.som, type = "dist.neighbours")
  # plot(test.som, main = "t")
  # graphics::par(mfrow = c(1, 1))

  return(som)
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
