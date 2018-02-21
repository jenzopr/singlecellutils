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
