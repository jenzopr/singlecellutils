#' Orders cells along a pseudotime between two or more known states.
#' This method implements the fraction of identity approach of Treutlein et. al. (2016)
#'
#' @param data The single-cell expression matrix
#' @param states A matrix containing the expression values for known states
#'
#' @return A list with the following elements:
#' identity - a matrix with one column per state, giving a numerical "identity" of each cell with that state, as well as the error per cell,
#' order - a matrix with one column per state, giving the decreasing order of the cells identity with that state.
#'
#' @export
calcFractionOfIdentity <- function(data, states) {
  if (nrow(data) != nrow(states)) {
    stop("Input matrices did not contain the same number of genes.")
  }

  num.states <- ncol(states)
  X <- as.matrix(states)

  # Calculate scaling for states
  # from http://stackoverflow.com/questions/28381855/r-function-solve-qp-error-constraints-are-inconsistent-no-solution
  M <- chol(t(X) %*% X)
  Dmat <- t(M) %*% M
  scaling <- norm(Dmat, "2")

  # Calculate fraction of identity for each cell
  identity <- t(apply(data, 2, function(x) {
    Y <- as.matrix(x)
    C <- cbind(rep(1, num.states), diag(num.states))
    b <- c(1, rep(0, num.states))
    d <- t(Y) %*% X
    QP <- quadprog::solve.QP(Dmat / scaling, d / scaling, C, b, meq = 1, factorized = FALSE)
    e <- sum(abs(Y - X %*% QP$solution))
    return(matrix(c(QP$solution, e), nrow = 1))
  }))

  # Build up
  rownames(identity) <- colnames(data)
  colnames(identity) <- c(colnames(states), "error")

  # Add cell ordering
  ordering <- apply(identity[, -ncol(identity)], 2, order, decreasing = TRUE)

  return(list(
    identity = identity,
    order = ordering
  ))
}
