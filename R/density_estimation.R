#' Discretize a matrix using the bayesian blocks algorithm
#'
#' @param X Numerical matrix to be discretized column-wise.
#' @param ... Additional parameters passed to bayesian_blocks
#'
#' @return A discretized matrix
#'
#' @export
discretize.bb <- function(X, ...) {
  X <- as.matrix(X)

  ret <- apply(X, 2, function(x, ...) {
    bins <- bayesian_blocks(x, ...)
    bins[which.min(bins)] <- min(x)
    as.numeric(cut(x, breaks = bins, include.lowest = T))
  }, ...)

  ret
}

#' Wrapper for astroML.density_estimation.bayesian_blocks
#'
#' @param t A numeric
#' @param x A numeric
#' @param sigma Data error
#' @param fitness Fitness function
#' @param kwargs A named list of ot <ther parameters
#'
#' @return array containing the (N+1) bin edges
#'
#' @export
bayesian_blocks <- function(t, x=NULL, sigma=NULL, fitness="events", kwargs = list()) {
  if (!is.numeric(t)) {
    t <- as.numeric(t)
  }

  rPython::python.assign("t", t)
  rPython::python.exec("x = None")
  rPython::python.exec("sigma = None")

  if (!is.null(x)) {
    x <- as.numeric(x)
    rPython::python.assign("x", x)
  }
  if (!is.null(sigma)) {
    sigma <- as.numeric(sigma)
    rPython::python.assign("sigma", sigma)
  }

  fi <- match.arg(fitness, c("events", "regular_events", "measures"))
  rPython::python.assign("fitness", fi)

  rPython::python.exec("from astroML.density_estimation import bayesian_blocks")
  rPython::python.exec("bins = bayesian_blocks(t=t, x=x, sigma=sigma, fitness=fitness)")
  bins <- rPython::python.method.call("bins","tolist")
  bins
}
