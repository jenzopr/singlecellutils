#' Adds results from a particular dimension reduction technique to the object.
#'
#' @param object A SingleCellExperiment object.
#' @param exprs_values String indicating which assay contains the data that should be used to calculate the fraction of identity.
#' @param features A character vector (of feature names), a logical vector or numeric vector (of indices) specifying the features to use for calculation. The default of NULL will use all features.
#' @param slot Determines which entry of the \code{reducedDims} slot to use for reduced embedding.
#' @param ... Additional parameters passed to \code{calcFractionOfIdentity}.
#'
#' @return A SingleCellExperiment object with modified \code{reducedDims} slot.
#'
#' @export
add_fraction_of_identity <- function(object, exprs_values = "counts", features = NULL, slot = 'fraction_of_identity', ...) {
  foi <- calcFractionOfIdentity(data = SummarizedExperiment::assay(object, i = exprs_values)[features,], ...)
  SingleCellExperiment::reducedDim(object, slot) <- foi$identity
  return(object)
}


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
