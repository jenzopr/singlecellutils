#' Calculate IC values for a data matrix (genes-by-cells)
#'
#' @param data Data matrix
#' @param p.val The p-value cutoff for correlation values to use in calculations
#'
#' @return The IC value
#' @export
computeIC <- function(data, p.val = 0.01) {
  d <- data

  r <- Hmisc::rcorr(t(d), type="pearson")
  use <- (r$P < p.val & !is.na(r$P) & upper.tri(r$r))
  within <- as.vector(abs(r$r[use]))

  r <- Hmisc::rcorr(d, type="pearson")
  use <- (r$P < p.val & !is.na(r$P) & upper.tri(r$r))
  between <- as.vector(r$r[use])

  mean(within) / mean(between)
}

#' Run non-parametric bootstrapping for IC values and groups
#'
#' @param data The expression matrix (genes-by-cells)
#' @param groups A vector giving group information (will be coerced to factor)
#' @param n The subsample size for non-parametric bootstrapping
#' @param R The bootstrapping size
#'
#' @return A matrix with all bootstrap samples per grouping level (R rows, length(levels(group)) columns)
#' @export
boot.ic <- function(data, groups, n = 20, R = 999, ...) {
  if (!is.factor(groups)) {
    groups <- factor(groups)
  }

  b.ic <- lapply(levels(groups), function(g) {
    d <- data[, groups == g]
    boot::boot(d, computeIC, sim="parametric", R=R, ran.gen = boot_subsample, mle=n, ...)
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
