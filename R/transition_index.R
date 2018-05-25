#' Adds the cirtical transition index to the object.
#'
#' @param object A SingleCellExperiment object.
#' @param column Determines the column name of the \code{colData} slot to store results to.
#' @param ... Additional parameters passed to \code{computeIC}.
#'
#' @return A SingleCellExperiment object with modified \code{colData} slot.
#'
#' @export
add_transition_index <- function(object, column, ...) {
  ic <- compute_transition_index(object, ...)
  object <- SingleCellExperiment::mutate(object, !!column := ic$ic)
  return(object)
}

#' Calculate the critical transition index for cells in a SingleCellExperiment object.
#' The critical transition index is described in detail in http://www.pnas.org/content/114/9/2271.long.
#'
#' @param object A SingleCellExperiment object.
#' @param exprs_values String indicating which assay contains the data that should be used to calculate the IC.
#' @param features A character vector (of feature names), a logical vector or numeric vector (of indices) specifying the features to use for calculation. The default of NULL will use all features.
#' @param pval The p-value cutoff for correlation values to use in calculations
#'
#' @return A list with the IC value (ic), as well as the used correlation values (within (GxG) / between (CxC))
#' @export
compute_transition_index <- function(object, exprs_values, features = NULL, pval = 0.01) {
  d <- SummarizedExperiment::assay(object, exprs_values)
  if(!is.null(features)) {
    d <- SummarizedExperiment::assay(object, exprs_values)[features, ]
  }

  r <- Hmisc::rcorr(t(d), type="pearson")
  use <- (r$P < p.val & !is.na(r$P) & upper.tri(r$r))
  within <- as.vector(abs(r$r[use]))

  r <- Hmisc::rcorr(d, type="pearson")
  use <- (r$P < p.val & !is.na(r$P) & upper.tri(r$r))
  between <- as.vector(r$r[use])

  list(ic = mean(within) / mean(between), within = within, between = between)
}

#' Run non-parametric bootstrapping for critical transition index values and a grouping variable
#'
#' @param object A SingleCellExperiment object.
#' @param exprs_values String indicating which assay contains the data that should be used to calculate IC.
#' @param grouping A single character indicating the column of \code{colData} to use as grouping, or a factor of \code{length > 1}.
#' @param features A character vector (of feature names), a logical vector or numeric vector (of indices) specifying the features to use for calculation. The default of NULL will use all features.
#' @param n The subsample size for non-parametric bootstrapping.
#' @param R The bootstrapping size.
#'
#' @return A matrix with all bootstrap samples per grouping level (\code{R} rows, \code{length(levels(group))} columns)
#' @export
bootstrap_transition_index <- function(object, exprs_values, grouping, features = NULL, n = 20, R = 999, ...) {
  if(length(grouping == 1)) {
    grouping <- SummarizedExperiment::colData(object)[[grouping]]
  }

  if (!is.factor(grouping)) {
    grouping <- factor(grouping)
  }

  b.ic <- lapply(levels(grouping), function(g) {
    d <- SummarizedExperiment::assay(object, exprs_values)[, grouping == g]
    if(!is.null(features)) {
      d <- SummarizedExperiment::assay(object, exprs_values)[features, grouping == g]
    }
    boot::boot(d, function(x, p.val) compute_transition_index(x, p.val)$ic, sim="parametric", R=R, ran.gen = boot_subsample, mle=n, ...)
  })

  ic_values <- sapply(b.ic, function(l) l$t)
  colnames(ic_values) <- levels(groups)
  return(ic_values)
}

#' Helper function for non-parametric subsampling
#'
#' @param x The data vector (or matrix)
#' @param subsample_size The size of the subsample
#'
#' @return A matrix in the same dimensions as x, but subsampled to subsample_size rows.
boot_subsample <- function(x, subsample_size) {
  x[sample(x = seq_len(nrow(x)), size = subsample_size), ]
}
