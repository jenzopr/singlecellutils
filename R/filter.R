#' Conveinience function to perform gene filtering using \code{genefilter}.
#'
#' @param object A SingleCellExperiment object.
#' @param flist A named list with gene filters. Names should be callable functions, list items should be named paramters.
#' @param exprs_values String indicating which assay contains the data that should be used for filtering.
#'
#' @return A SingleCellExperiment object with features passing the filter.
#'
#' @export
filter_genes <- function(object, flist, exprs_values = "counts") {
  filter_functions <- lapply(names(flist), function(f) {
    fn <- strsplit(f, "::")[[1]]
    what <- if(length(fn)==1) {
      get(fn[[1]], envir=parent.frame(), mode="function")
    } else {
      get(fn[[2]], envir=asNamespace(fn[[1]]), mode="function")
    }
    do.call(what = what, args = flist[[f]])
  })

  keep <- genefilter::genefilter(SummarizedExperiment::assay(object, i = exprs_values), genefilter::filterfun(filter_functions))
  return(object[keep, ])
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
    var(x) / mean(x) > A
  }
}
