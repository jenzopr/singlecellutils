#' Thinning of a expression matrix from a SingleCellExperiment object.
#'
#' @param object A SingleCellExperiment object.
#' @param exprs_values String indicating which assay contains the data that should be used to perform thinning.
#' @param current_libsizes The current library sizes (can directly be calculated from \code{object}).
#' @param target_libsizes The target library sizes th thin expression to.
#' @param return_as String indicating the name of the assay which contains the thinned expression values.
#'
#' @return A SingleCellExperiment object with a new assay containing the thinned values.
#'
#' @export
thin_counts <- function(object, exprs_values, current_libsizes = scater::librarySizeFactors(object, exprs_values),
                       target_libsizes = min(current_libsizes), return_as = "downsampled") {
  y <- SummarizedExperiment::assay(object, exprs_values)
  n <- dim(y)[1];
  m <- dim(y)[2];

  if (any (target_libsizes > current_libsizes)) {
    stop("Error: Found a target library size that is larger than its current library size.");
  }

  if (length(target_libsizes) == 1) {
    target.lib.sizes = rep(target.lib.sizes, m);
  }

  y_thinned = y;

  for (k in 1:m) {
    keep.p =  target_libsizes[k]/current_libsizes[k]; # proportion of reads to keep
    if (keep.p == 1) {
      y_thinned[, k] <- y[,k]
    }
    else {
      y_thinned[,k] = rbinom(n, y[,k], keep.p);
    }
  }

  SummarizedExperiment::assay(object, return_as) <- y_thinned
  return(object)
}
