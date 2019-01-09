#' Calculates the term frequency-inverse document frequency.
#'
#' @param object A SingleCellExperiment object.
#' @param exprs_values String indicating which assay contains the data that should be used to calculate TF-IDF.
#' @param features A character vector (of feature names), a logical vector or numeric vector (of indices) specifying the features to use for calculation The default of NULL will use all features.
#' @param tf.method A character vector indicating the variant for calculating the TF.
#' @param idf.weight A character vector indicating the weighting for IDF calculation.
#'
#' @return A \code{\link[Matrix]{sparseMatrix}}
#'
#' @export
tf_idf <- function(object, exprs_values = "counts", features = NULL, tf.method = c("raw", "adjusted", "logscaled"), idf.weight = c("unary", "idf", "idf.smooth", "pidf")) {
  input <- as.matrix(SummarizedExperiment::assay(object, i = exprs_values))
  if (!is.null(features)) {
    input <- as.matrix(SummarizedExperiment::assay(object, i = exprs_values)[features,])
  }

  tf <- switch(tf.method,
               "raw" = input,
               "adjusted" = t(t(input) / Matrix::colSums(input)),
               "logscaled" = log(input + 1))

  idf <- switch(idf.weight,
                "unary" = rep(1, nrow(input)),
                "idf" = -log(Matrix::rowSums(input) / ncol(input)),
                "idf.smooth" = log(1 + ncol(input) / Matrix::rowSums(input)),
                "pidf" = log((ncol(input) - Matrix::rowSums(input)) / Matrix::rowSums(input)))
  idf <- methods::as(idf, "sparseVector")
  methods::as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% tf
}
