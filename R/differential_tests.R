#' Tests differential expression/accessibility using logistic regression
#'
#' This function implements the test proposed by Ntranos et. al. (doi: 10.1038/s41592-018-0303-9) for \code{\link[SingleCellExperiment]{SingleCellExperiment}} objects.
#' Compared to the traditional approach of using cell labels as covariate, logistic regression is carried out to predict cell labels from
#' quantification of associated features. Those features can be transcripts or accessibility peaks associated with e.g. genes or regulatory regions.
#'
#' @param object A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object.
#' @param response Character vector indicating the response (class) variable. In case of length of one, it should determines the column name of the \code{\link[SummarizedExperiment]{colData}} slot.
#' @param exprs_values String indicating which assay contains the data that should be used to perform grouping and testing.
#' @param group_by Character vector indicating the grouping variable for expression values. In case of length of one, it should determines the column name of the \code{\link[SummarizedExperiment]{rowData}} slot.
#' @param reference The level of \code{response} that should be used as reference.
#' @param features A character vector (of feature names), a logical vector or numeric vector (of indices) specifying the features to test. The default of NULL will use all features.
#' @param ... Additional parameters passed to \code{\link[stats]{p.adjust}} for multiple testing correction.
#'
#' @return A tibble with the following columns:
#'   - group - the grouping variable
#'   - pvalue - the p-value from comparison of the two LR models
#'   - p.adjusted - the multiple testing corrected \code{pvalue}.
#'
#' @export
lrde_test <- function(object, response, exprs_values = "logcounts", group_by = NULL, reference = NULL, features = NULL, ...) {
  input <- as.matrix(SummarizedExperiment::assay(object, i = exprs_values))
  if (!is.null(features)) {
    input <- as.matrix(SummarizedExperiment::assay(object, i = exprs_values)[features,])
  }

  # Check / complete input
  if(length(response) == 1) {
    response <- SummarizedExperiment::colData(object)[, response]
  }

  if(!is.null(group_by)) {
    if(length(group_by) == 1) {
      if(!is.null(features)) {
        group_by <- SummarizedExperiment::rowData(object)[features, group_by]
      } else {
        group_by <- SummarizedExperiment::rowData(object)[, group_by]
      }
    }

    data.frame(group = group_by, index = 1:length(group_by)) %>%
      dplyr::filter(!is.na(group)) %>%
      dplyr::group_by(group) %>%
      dplyr::do(pvalue = .lrde_pval(t(input[.$index, ]), response, reference)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(pvalue = unlist(pvalue)) -> result
  } else {
    result <- data.frame(pvalue = .lrde_pval(t(input), response, reference))
  }

  result %>%
    dplyr::mutate(p.adjusted = stats::p.adjust(pvalue, ...)) %>%
    dplyr::arrange(p.adjusted)
}

#' Function to perform logistic regression test
#'
#' @param data Matrix with the predictor variables in columns
#' @param response Vector with two levels indicating the response (class) variable.
#' @param response.ref The level of \code{response} that should be used as reference.
#'
#' This function needs the package \code{lmtest} to compare GLMs.
#' In case no response reference level is set, but the response contains more than two levels, a warning is reported and the first level is automatically chosen.
#'
#' @return A p-value of the \code{\link[lmtest]{lrtest}} between the full and null model.
.lrde_pval <- function(data, response, response.ref = NULL) {
  if (!requireNamespace("lmtest", quietly = TRUE)) {
    stop("Package lmtest needed for this function to work. Please install it.", call. = FALSE)
  }

  if(!is.factor(response)) {
    response <- factor(response)
  }
  if(length(levels(response)) < 2) {
    stop(paste("Response variable can't have less than two levels:", levels(response), collapse = T))
  }
  if(length(levels(response)) > 2) {
    if(is.null(response.ref)) {
      warning("No response reference set, will use the first level of the response variable.")
      response.ref <- levels(response)[1]
    }
    non_response.ref <- paste0("non-", response.ref)
    response <- factor(ifelse(response == response.ref, response.ref, non_response.ref), levels = c(response.ref, non_response.ref))
  }

  if(dim(data)[1] == 1) {
    data <- as.numeric(data)
  }

  full_model <- stats::glm(response ~ data, family = stats::binomial())
  null_model <- stats::glm(response ~ 1, family = stats::binomial())

  test <- lmtest::lrtest(full_model, null_model)
  return(test$Pr[2])
}


