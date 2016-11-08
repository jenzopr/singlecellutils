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
