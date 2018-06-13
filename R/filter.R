#' Conveinience function to perform feature filtering using \code{genefilter}.
#'
#' @param object A SingleCellExperiment object.
#' @param flist A named list with gene filters. Names should be callable functions, list items should be named paramters.
#' @param exprs_values String indicating which assay contains the data that should be used for filtering.
#' @param tolerate A number indicating how many failed filters should be toleated.
#'
#' @return A SingleCellExperiment object with features passing the filter(s).
#'
#' @export
filter_features <- function(object, flist, exprs_values = "counts", tolerate = 0) {
  filter_functions <- construct_filters(flist)

  keep <- apply(SummarizedExperiment::assay(object, i = exprs_values), 1, n_filterfun(n = tolerate, filter_functions))
  SingleCellExperiment:::int_metadata(object)$feature_filter <- flist
  return(object[keep, ])
}

#' Conveinience function to perform sample filtering using \code{genefilter}.
#'
#' @param object A SingleCellExperiment object.
#' @param flist A named list with gene filters. Names should be callable functions, list items should be named paramters.
#' @param exprs_values String indicating which assay contains the data that should be used for filtering.
#' @param tolerate A number indicating how many failed filters should be toleated.
#'
#' @return A SingleCellExperiment object with samples passing the filter(s).
#'
#' @export
filter_samples <- function(object, flist, exprs_values = "counts", tolerate = 0) {
  filter_functions <- construct_filters(flist)

  keep <- apply(t(SummarizedExperiment::assay(object, i = exprs_values)), 1, n_filterfun(n = tolerate, filter_functions))
  SingleCellExperiment:::int_metadata(object)$sample_filter <- flist
  return(object[, keep])
}

#' EpOverA returns a filter function with bindings for A. This function evaluates to TRUE if the aggregate of the arguments elements are larger than A.
#'
#' @param A The value you want to exceed.
#' @param na.rm Whether NAs should be removed.
#'
#' @return A function with bindings for A.
#'
#' @export
EpOverA <- function(A = 4, na.rm = TRUE) {
  function(x) {
    if (na.rm)
      x <- x[!is.na(x)]
    sum(x) > A
  }
}

#' FanoOverA returns a filter function with bindings for A. This function evaluates to TRUE if the Fano factor of the argument is larger than A.
#'
#' @param A The value that should be exceeded.
#' @param na.rm Whether NAs should be removed.
#'
#' @return A function with bindings for A.
#'
#' @details See Weinreb et.al., https://doi.org/10.1093/bioinformatics/btx792
#'
#' @export
FanoOverA <- function(A = 2, na.rm = TRUE) {
  function(x) {
    if (na.rm)
      x <- x[!is.na(x)]
    stats::var(x) / mean(x) > A
  }
}

#' scaterIsOutlier returns a filter function with bindings for parameters for \code{\link[scater]{isOutlier}}.
#'
#' @param ... Parameters that will be bound to the \code{\link[scater]{isOutlier}} code.
#'
#' @return A function with ... bindings for \code{\link[scater]{isOutlier}}.
#'
#' @export
scaterIsOutlier <- function(...) {
  function(x) {
    scater::isOutlier(metric = x, ...)
  }
}

#' A custom \code{filterfun} that creates a n-FALSE exiting function from its input.
#'
#' @param n The number of tolerated FALSE exiting functions before the returned function returns FALSE. The behaviour of the default (n=0) resembles \code{\link[genefilter]{filterfun}}.
#' @param ... A list of filter functions.
#'
#' @return A function.
#'
#' @seealso \code{\link[genefilter]{filterfun}}.
n_filterfun <- function (n = 0, ...) {
  flist <- list(...)
  if (length(flist) == 1 && is.list(flist[[1]]))
    flist <- flist[[1]]
  f <- function(x) {
    fails <- 0
    for (fun in flist) {
      fval <- fun(x)
      if (is.na(fval) || !fval)
        fails <- fails + 1
        if (fails > n) return(FALSE)
    }
    return(TRUE)
  }
  return(f)
}

#' Conveneince function to construct a list of filter function(s) from a named list
#'
#' @param flist A named list with gene filters. Names should be callable functions, list items should be named paramters.
#'
#' @return A list.
construct_filters <- function(flist) {
  lapply(names(flist), function(f) {
    fn <- strsplit(f, "::")[[1]]
    what <- if(length(fn)==1) {
      get(fn[[1]], envir=parent.frame(), mode="function")
    } else {
      get(fn[[2]], envir=asNamespace(fn[[1]]), mode="function")
    }
    do.call(what = what, args = flist[[f]])
  })
}
