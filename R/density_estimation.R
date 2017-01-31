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

#' Calculates the redundant contribution of x1 and x2 regarding the target variable z
#'
#' If we consider S = {X1 , X2 } and a target variable Z, the calculates the amount of information provided by each variable within S about each state of the target Z.
#'
#' @param z The target varibale
#' @param x1 The first source variable
#' @param x2 The second source variable
#'
#' @return The amount of information provided by each variable within S about each state of the target Z
#'
#' @export
redundancy <- function(z, x1, x2) {
  prob_z <- sapply(z, function(x) sum(z==x)/length(z))
  min_ispec <- sapply(z, function(x) min(ispec(x1, z, x),ispec(x2, z, x)))
  sum(prob_z * min_ispec)
}

#' Calculates the specific information which quantifies the information provided by one variable about a specific state of another variable [89, 90]
#'
#' @param x The variable to quantify the information from
#' @param y The other variable
#' @param state The specific state of y
#'
#' @return The information provided by x about a specific state of y
ispec <- function(x, y, state) {
  xy <- condprop(x,y)
  p_x_state <- xy[as.character(x),state]
  p_state <- sum(y==state)/length(y)

  yx <- condprop(y,x)
  p_state_x <- yx[state,as.character(x)]
  ispec <- p_x_state * (log(1/p_state) - log(1/p_state_x))
  ispec[!is.finite(ispec)] <- 0
  sum(ispec)
}

#' Calculates the conditional propability p(x|given) for two vectors
#'
#' @param x The first vector
#' @param given The second vector
#'
#' @return A vector with conditional propabilities.
condprop <- function(x,given) {
  t <- as.matrix(table(x,given))
  xy <- prop.table(t, 2)
  xy
  #diag(xy[as.character(x), as.character(given)])
}
