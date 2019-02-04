#' Add a clustering
#'
#' Adds results from a particular clustering technique to the object.
#'
#' @param object A \linkS4class{SingleCellExperiment::SingleCellExperiment} object.
#' @param flavor Determines which clustering technique to apply.
#' @param column Determines the column name of the \code{\link[SummarizedExperiment]{colData}} slot to store results to.
#' @param ... Additional parameters passed to functions.
#'
#' @return A SingleCellExperiment object with modified \code{\link[SummarizedExperiment]{colData}} slot.
#'
#' @seealso [hdbscan()] for HDBSCAN clustering.
#' @export
add_clustering <- function(object, flavor = c("hdbscan"), column = ".cluster", ...) {
  clustering <- switch(flavor,
                       hdbscan = hdbscan(object, ...))
  object <- scater::mutate(object, !!column := clustering)
  return(object)
}

#' Conveinience function to perform HDBSCAN via reticulate.
#'
#' Performs HDBSCAN clustering on a SingleCellExperiment object.
#' The algorithm is explained in detail in \href{http://joss.theoj.org/papers/10.21105/joss.00205}{McInnes et. al.} and needs the python package \href{https://pypi.org/project/hdbscan/}{hdbscan} installed.
#'
#' @param object A SingleCellExperiment object.
#' @param use_dimred A string or integer scalar indicating the reduced dimension result in \code{reducedDims(object)} to use as input.
#' @param min_samples Measure of how conservative the clustering should to be. The larger the value of \code{min_samples}, the more conservative the clustering and more points will be declared as noise, and clusters will be restricted to progressively more dense areas.
#' @param min_cluster_size The smallest size grouping that is considered a cluster.
#' @param outlier Determines how outliers are encoded in the resulting clustering.
#' @param seed A numeric seed to initialize the random number generator.
#'
#' @return A factor with the assigned cluster.
#'
#' @export
hdbscan <- function(object, use_dimred, min_samples = 7L, min_cluster_size = 9L, outlier = 0, seed = NULL) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package reticulate needed for this function to work. Please install it.", call. = FALSE)
  }

  if (!is.null(seed)) set.seed(seed)

  h <- reticulate::import("hdbscan")
  cl <- h$HDBSCAN(min_samples = min_samples, min_cluster_size = min_cluster_size)

  labels <- cl$fit_predict(SingleCellExperiment::reducedDim(object, use_dimred)) + 1
  labels[labels == 0] <- outlier
  return(factor(labels))
}

#' Conveinience function to perform HDBSCAN via reticulate.
#'
#' Performs HDBSCAN clustering on a n-dimensional matrix.
#' The algorithm is explained in detail in http://joss.theoj.org/papers/10.21105/joss.00205 and needs the python package \code{hdbscan} installed.
#'
#' @param embedding A matrix giving the embedding to be used for clustering.
#' @param min_samples Measure of how conservative the clustering should to be. The larger the value of \code{min_samples}, the more conservative the clustering and more points will be declared as noise, and clusters will be restricted to progressively more dense areas.
#' @param min_cluster_size The smallest size grouping that is considered a cluster.
#' @param outlier Determines how outliers are encoded in the resulting clustering.
#' @param seed A numeric seed to initialize the random number generator.
#'
#' @return A factor with the assigned cluster.
#'
#' @export
hdbscan_ <- function(embedding, min_samples = 7L, min_cluster_size = 9L, outlier = 0, seed = NULL) {
  if(!is.matrix(embedding)) {
    embedding <- as.matrix(embedding)
  }

  if (!is.null(seed)) set.seed(seed)

  h <- reticulate::import("hdbscan")
  cl <- h$HDBSCAN(min_samples = min_samples, min_cluster_size = min_cluster_size)

  labels <- cl$fit_predict(embedding) + 1
  labels[labels == 0] <- outlier
  return(factor(labels))
}
