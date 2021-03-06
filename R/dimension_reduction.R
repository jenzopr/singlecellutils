#' Performs dimension reduction and adds results to the given input.
#'
#' See details for available dimension reduction techniques that can be applied to the object.
#'
#' The following dimension reduction techniques are available:
#' \itemize{
#' \item \emph{UMAP} - Uniform Manifold Approximation and Projection, see \link{umap}.
#' \item \emph{approximate SVD} - see \link{approximate_svd}
#' \item \emph{randomized SVD} - see \link{randomized_svd}
#' \item \emph{Multidimensional scaling} - see \link{isomds}
#' \item \emph{Self-organizing maps} - see \link{scSOM}
#' }
#'
#' @param object A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object.
#' @param flavor Determines which dimension reduction technique to apply.
#' @param slot Determines which entry of the \code{\link[SingleCellExperiment]{reducedDims}} slot to use for reduced embedding.
#' @param ... Additional parameters passed to functions.
#'
#' @return A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object with modified \code{reducedDims} slot.
#'
#' @examples
#' \dontrun{
#' # Apply UMAP on normalized expression values, but store results in a slot called dimred_umap.
#' obj <- reduce_dimension(obj, flavor = "umap", slot = "dimred_umap", exprs_values = "norm_exprs")
#'
#' # Apply randomized SVD on normalized expression values, keeping 50 dimensions, followed by umap.
#' # Using magrittr pipes, the result of randomized SVD is stored in the rsvd slot and
#' # picked up by umap.
#' obj %<>%
#'   reduce_dimension(flavor = "rsvd", exprs_values = "norm_exprs", n_dims = 50) %>%
#'   reduce_dimension(flavor = "umap", exprs_values = NULL, use_dimred = "rsvd")
#' }
#' @export
reduce_dimension <- function(object, flavor = c("umap", "som", "asvd", "isomds", "rsvd"), slot = flavor, ...) {
  embedding <- switch(flavor,
    umap = umap(object, ...),
    som = t(scSOM(object, ...)$codes[[1]]),
    asvd = approximate_svd(object, ...),
    isomds = isomds(object, ...),
    rsvd = randomized_svd(object, ...))
  SingleCellExperiment::reducedDim(object, slot) <- embedding
  return(object)
}

#' Conveinience function to perform UMAP via reticulate.
#'
#' Performs Uniform Manifold Approximation and Projection (UMAP) on a SingleCellExperiment object.
#' The algorithm is explained in detail in \href{https://arxiv.org/abs/1802.03426}{McInnes et. al.} and requires the python packages \href{https://pypi.org/project/umap-learn/}{umap-learn} and \href{https://pypi.org/project/numpy/}{numpy} to be installed.
#'
#' @param object A SingleCellExperiment object.
#' @param exprs_values String indicating which assay contains the data that should be used to perform UMAP.
#' @param features A character vector (of feature names), a logical vector or numeric vector (of indices) specifying the features to use for UMAP. The default of NULL will use all features.
#' @param scale Logical scalar indicating whether data should be scaled before UMAP.
#' @param n_neighbors This determines the number of neighboring points used in local approximations of manifold structure. Larger values will result in more global structure being preserved at the loss of detailed local structure. In general this parameter should often be in the range 5 to 50, with a choice of 10 to 15 being a sensible default.
#' @param min_dist This controls how tightly the embedding is allowed compress points together. Larger values ensure embedded points are more evenly distributed, while smaller values allow the algorithm to optimise more accurately with regard to local structure. Sensible values are in the range 0.001 to 0.5, with 0.1 being a reasonable default.
#' @param metric This determines the choice of metric used to measure distance in the input space.
#' @param use_dimred String indicating which slot in \code{reducedDim} contains the data that should be used to perform UMAP (if \code{exprs_values = NULL}).
#' @param seed A numeric seed to initialize the random number generator.
#'
#' @return A matrix with the two-dimensional embedding.
#'
#' @export
umap <- function(object, exprs_values, features = NULL, scale = T, n_neighbors = 5L, min_dist = 0.1, metric = "correlation", use_dimred = NULL, seed = NULL) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package reticulate needed for this function to work. Please install it.", call. = FALSE)
  }

  if (!is.null(seed)) set.seed(seed)
  umap <- reticulate::import("umap")
  np <- reticulate::import("numpy", convert=FALSE)
  umap_cl <- umap$UMAP(n_neighbors = as.integer(n_neighbors), min_dist = min_dist, metric = metric)

  if (!is.null(exprs_values)) {
    input <- as.matrix(SummarizedExperiment::assay(object, i = exprs_values))
    if (!is.null(features)) {
      input <- as.matrix(SummarizedExperiment::assay(object, i = exprs_values)[features,])
    }
  } else {
    if (!is.null(use_dimred)) {
      input <- t(as.matrix(SingleCellExperiment::reducedDim(object, type = use_dimred)))
    } else {
      stop("Both exprs_values and use_dimred can't be NULL.")
    }
  }
  if (scale) input <- t(scale(t(input)))

  embedding <- umap_cl$fit_transform(t(input))
  return(embedding)
}

#' Conveinience function to run isoMDS on a SingleCellExperiment object.
#'
#' @param object A SingleCellExperiment object.
#' @param exprs_values String indicating which assay contains the data that should be used to perform isoMDS.
#' @param features A character vector (of feature names), a logical vector or numeric vector (of indices) specifying the features to use for isoMDS. The default of NULL will use all features.
#' @param method A character string specifying the method to be used to calulate a dissimilarity structure using WGCNA::cor.
#' @param ... Additional parameters to be passed on to MASS::isoMDS.
#'
#' @return A matrix with the k-dimensional embedding.
#'
#' @export
isomds <- function(object, exprs_values, features = NULL, method = "spearman", ...) {
  input <- as.matrix(SummarizedExperiment::assay(object, i = exprs_values))
  if (!is.null(features)) {
    input <- as.matrix(SummarizedExperiment::assay(object, i = exprs_values)[features,])
  }

  dissimilarity <- (1 - WGCNA::cor(input, method = method))
  return(MASS::isoMDS(dissimilarity, ...)$points)
}

#' Approximate or randomized singular value decomposition
#'
#' Approximate SVD is performed using the implementation from \code{\link[irlba]{irlba}}, randomized SVD is performed using the implementation from \code{\link[rsvd]{rsvd}}.
#'
#' @param object A SingleCellExperiment object.
#' @param exprs_values String indicating which assay contains the data that should be used to perform SVD.
#' @param n_dims The number of approximate singular values to calculate.
#' @param features A character vector (of feature names), a logical vector or numeric vector (of indices) specifying the features to use for SVD. The default of NULL will use all features.
#' @param skip A numeric vector indicating which singular values to set to zero (and remove).
#' @param seed A numeric seed to initialize the random number generator.
#' @param ... Additional arguments passed on to \code{\link[rsvd]{rsvd}} or \code{\link[irlba]{irlba}}.
#'
#' @return A matrix of dimension \code{ncol(object)} x \code{dims(object) - length(skip)}.
#' @name svd
NULL
#> NULL

#' @rdname svd
#' @export
randomized_svd <- function(object, exprs_values, n_dims, features = NULL, skip = NULL, seed = NULL, ...) {
  if (!requireNamespace("rsvd", quietly = TRUE)) {
    stop("Package rsvd needed for this function to work. Please install it.", call. = FALSE)
  }

  if (!is.null(seed)) set.seed(seed)

  input <- as.matrix(SummarizedExperiment::assay(object, i = exprs_values))
  if (!is.null(features)) {
    input <- as.matrix(SummarizedExperiment::assay(object, i = exprs_values)[features,])
  }

  svd <- rsvd::rsvd(input, k = n_dims, ...)
  diag_s <- matrix(0, nrow = n_dims, ncol = n_dims)
  diag(diag_s) = svd$d

  if (!is.null(skip)) {
    diag_s[skip,skip] <- 0
    t(diag_s %*% t(svd$v))[, -skip]
  } else {
    # return S %*% t(V)
    t(diag_s %*% t(svd$v))
  }
}

#' @rdname svd
#' @export
approximate_svd <- function(object, exprs_values, n_dims, features = NULL, skip = NULL, seed = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)

  input <- as.matrix(SummarizedExperiment::assay(object, i = exprs_values))
  if (!is.null(features)) {
    input <- as.matrix(SummarizedExperiment::assay(object, i = exprs_values)[features,])
  }

  svd <- irlba::irlba(input, nv = n_dims, nu = n_dims, ...)
  diag_s <- matrix(0, nrow = n_dims, ncol = n_dims)
  diag(diag_s) = svd$d

  if (!is.null(skip)) {
    diag_s[skip,skip] <- 0
    t(diag_s %*% t(svd$v))[, -skip]
  } else {
    # return S %*% t(V)
    t(diag_s %*% t(svd$v))
  }
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
#' @param ... Additional arguments passed to \code{\link[kohonen]{som}}.
#'
#' @return A som object.
#'
#' @export
scSOM <- function(object, exprs_values = 'norm_TPM', features = NULL, scale = T, num_epochs = 200, resolution = 1, seed = NULL, init = NULL, init.FUN = map.init, ...) {
  if (!is.null(features)) {
    input <- as.matrix(SummarizedExperiment::assay(object, i = exprs_values)[features, ])
  } else {
    input <- as.matrix(SummarizedExperiment::assay(object, i = exprs_values)[1:min(nrow(object), 5000), ])
  }
  if (scale) input <- t(scale(t(input)))

  if (num_epochs < nrow(input) / 25) {
    lower_bound <- floor(nrow(input) / 25)
    warning(paste(num_epochs, "epochs is low for", nrow(input), "training genes, it was set to",
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

  pcar <- stats::prcomp(data)
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

  pcar <- stats::prcomp(data)
  gridshape <- array(0, dim = c(h * w, ncol(data)))

  h_pos <- as.numeric(cut(pcar$x[,1], breaks=h))
  w_pos <- as.numeric(cut(pcar$x[,2], breaks=w))
  cell_pos <- (h_pos - 1) * w + w_pos
  n <- constructHexNeighborhood(h, w, even.layout = T)

  h_genes <- rownames(data)[order(abs(pcar$rotation[,1]), decreasing = T)[1:100]]
  w_genes <- rownames(data)[order(abs(pcar$rotation[,2]), decreasing = T)[1:100]]
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

#' Calculates a 2D QC map of cells based on quality features
#'
#' @param object A SingleCellExperiment object.
#' @param features A character vector with column names form \code{colData(object)} that are used to calculate the map.
#' @param slot Determines which entry of the \code{reducedDims} slot to use for QC embedding.
#' @param ... Additional parameters passed onto \code{\link[Rtsne]{Rtsne}}.
#'
#' @return A SingleCellExperiment object with modified \code{reducedDims} slot.
#'
#' @export
calculate_qc_map <- function(object, features, slot = "qcmap", ...) {
  valid_features <- features %in% colnames(SummarizedExperiment::colData(object))

  if(all(!valid_features)) {
    stop("None of features is present ib object. Please choose valid column names.")
  }
  if(!all(valid_features)) {
    warning(paste("Not all features are present in object. Using only:", paste(features[valid_features], collapse = ",")))
  }

  qcdata <- as.matrix(SummarizedExperiment::colData(object)[, features[valid_features]])
  tsne <- Rtsne::Rtsne(qcdata, ...)

  SingleCellExperiment::reducedDim(object, slot) <- tsne$Y
  return(object)
}
