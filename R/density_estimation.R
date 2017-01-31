#' Wrapper for astroML.density_estimation.bayesian_blocks
#'
#' @param t A numeric
#' @param x A numeric
#' @param sigma Data error
#' @param fitness Fitness function
#' @param ... Other parameters
#'
#' @return array containing the (N+1) bin edges
#'
#' @export
bayesian_blocks <- function(t, x=NULL, sigma=NULL, fitness="events", ...) {
  if (!is.numeric(t)) {
    t <- as.numeric(t)
  }

  rPython::python.assign("t", t)
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
  rPython::python.exec("bins = bayesian_blocks(t, fitness=fitness)")
  bins <- rPython::python.method.call("bins","tolist")
  bins
}
