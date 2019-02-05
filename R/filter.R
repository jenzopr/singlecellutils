#' Functions to filter features or samples using \code{\link[genefilter]{genefilter}}.
#'
#' The functions will assemble a (combined) filtering function from the inputs and then apply it to features/samples from the given \code{\link[SingleCellExperiment]{SingleCellExperiment}} object.
#'
#' The list of filters should be of the format \code{list("package::function" = list(...))}. The filter functions described here will try to fetch the correct function from the associated package namespace before evaluation.
#' In case several filters are given, the filter functions allow to tolerate up to \code{tolerate} \emph{failing} filters. This is useful to implement \emph{OR} conjuctions between filters.
#'
#' @param object A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object.
#' @param filters A named list. Names should be callable (filter) functions, list items should be named paramters.
#' @param exprs_values String indicating which \code{\link[SummarizedExperiment]{assay}} contains the data that should be used for filtering.
#' @param tolerate A number indicating how many failed filters should be toleated.
#'
#' @return A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object with features or samples passing the filter(s).
#' @name filter
#'
#' @examples
#' \dontrun{
#' # Keep only features with an aggregated expression above 1000 counts.
#' obj <- filter_features(obj, list("EpOverA" = list(A = 1000)))
#'
#' # Keep only features with an aggregated expression above 1000 counts and expression
#' # in at least in 10 samples.
#' obj <- filter_features(obj, list("EpOverA" = list(A = 1000), "genefilter::kOverA" = list(k=10, A=0)))
#'
#' # Keep only features with an aggregated expression above 1000 counts or expression
#' # in at least in 10 samples (tolerate one being not met).
#' obj <- filter_features(obj, list("EpOverA" = list(A = 1000), "genefilter::kOverA" = list(k=10, A=0)), tolerate = 1)
#'
#' # Keep only features with expression over 10 in at least 5 samples, using the fpkm assay.
#' obj <- filter_features(obj, list("genefilter::kOverA" = list(k=5, A=10)), exprs_values = "fpkm")
#'
#' # Keep only samples that express at least 1000 features.
#' obj <- filter_samples(obj, list("genefilter::kOverA" = list(k = 1000, A=0))
#'
#' # Keep only samples with at least 100000 counts across all features.
#' obj <- filter_samples(obj, list("EpOverA" = list(A = 100000)))
#' }
NULL
#> NULL

#' @rdname filter
#' @export
filter_features <- function(object, filters, exprs_values = "counts", tolerate = 0) {
  filter_functions <- construct_filters(filters)

  keep <- apply(SummarizedExperiment::assay(object, i = exprs_values), 1, n_filterfun(n = tolerate, filter_functions))
  SingleCellExperiment:::int_metadata(object)$feature_filter <- filters
  return(object[keep, ])
}

#' @rdname filter
#' @export
filter_samples <- function(object, filters, exprs_values = "counts", tolerate = 0) {
  filter_functions <- construct_filters(filters)

  keep <- apply(Matrix::t(SummarizedExperiment::assay(object, i = exprs_values)), 1, n_filterfun(n = tolerate, filter_functions))
  SingleCellExperiment:::int_metadata(object)$sample_filter <- filters
  return(object[, keep])
}

#' Basic filter functions
#'
#' Basic filter functions return functions with bindings for A. The function evaluates to \code{TRUE} if some operation on the input results in a value greater/smaller than A.
#'
#' The basic filter functions in detail:
#' \itemize{
#' \item \code{EpOverA} evaluates to \code{TRUE} if the aggregate of the arguments elements are larger than A.
#' \item \code{FanoOverA} evaluates to \code{TRUE} if the Fano factor of the argument is larger than A.
#' \item \code{MeanBelowA} evaluates to \code{TRUE} if the mean of the arguments elements are smaller than A.
#' \item \code{MeanOverA} evaluates to \code{TRUE} if the mean of the arguments elements are greater than A.
#' }
#'
#' @param A The value that should / should not be exceeded.
#' @param na.rm Whether NAs should be removed.
#'
#' @return A function with bindings for A.
#' @name filterfunctions
NULL
#> NULL

#' @rdname filterfunctions
#' @export
EpOverA <- function(A = 4, na.rm = TRUE) {
  function(x) {
    if (na.rm)
      x <- x[!is.na(x)]
    sum(x) > A
  }
}

#' @rdname filterfunctions
#' @export
FanoOverA <- function(A = 2, na.rm = TRUE) {
  function(x) {
    if (na.rm)
      x <- x[!is.na(x)]
    stats::var(x) / mean(x) > A
  }
}

#' @rdname filterfunctions
#' @export
MeanBelowA <- function(A = 4, na.rm = TRUE) {
  function(x) {
    if (na.rm)
      x <- x[!is.na(x)]
    mean(x) < A
  }
}

#' @rdname filterfunctions
#' @export
MeanOverA <- function(A = 4, na.rm = TRUE) {
  function(x) {
    if (na.rm)
      x <- x[!is.na(x)]
    mean(x) > A
  }
}

#' A filter function with bindings for parameters for \code{\link[scater]{isOutlier}}.
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

#' Convenience function to discover elements in \code{what} in a vector \code{where} using a function \code{FUN}.
#'
#' @param object A SingleCellExperiment object.
#' @param what A character vector/character with what to find
#' @param where Where to search the object
#' @param FUN Function name used to perform the search
#' @param column Column name to search in (if rowData or colData should be searched).
#'
#' @return A list with TRUE/FALSE
#'
#' @export
discover_element <- function(object, what, where = c("rownames", "colnames", "rowData", "colData"), FUN = c("grepl", "is_in"), column = NULL) {
  if(length(where) == 1) {
    x <- switch(where,
                rownames = rownames(object),
                colnames = colnames(object),
                rowData = SummarizedExperiment::rowData(object)[, column],
                colData = SummarizedExperiment::colData(object)[, column],
                where)
  } else {
    x <- where
  }

  fn <- switch(FUN,
               grepl = get("grepl", envir = asNamespace("base"), mode = "function"),
               is_in = get("is_in", envir = asNamespace("magrittr"), mode = "function"))
  args <- list(what, x)
  do.call(what = fn, args = args)
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
