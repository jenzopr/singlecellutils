#' Calculate a weight matrix using false-negative curves
#'
#' @param data Expression data matrix
#' @param condition A factor describing different conditions
#' @param means A matrix containing mean expression values, e.g. from bulk measurements
#' @param expression_cutoff A numerical value
#' @param supress.plot Whether to supress plotting (default TRUE)
#'
#' @export
calcFNWeight <- function(data, condition = NULL, means = NULL, expression_cutoff = 10, supress.plot = T) {
  if (!is.numeric(data))
    stop("Argument data should be numeric.")
  data <- as.matrix(data)
  n <- nrow(data)
  m <- ncol(data)

  if (is.null(condition)) {
    condition <- factor(rep("a", m))
  } else {
    condition <- as.factor(condition)
  }

  if (is.null(means)) {
    means <- calcMeansPerCondition(data, condition)
  } else {
    means <- as.matrix(means)
    if (ncol(means) != length(levels(condition)))
      stop("Mean matrix should have the same number of columns as conditions has levels.")
    colnames(means) <- levels(condition)
  }

  theils.t <- norm.theil.t(data, winsorize = T, treat.zeros = "eps")
  residuals <- scale(theils.t$theils.t - (exp(theils.t$i + theils.t$d * theils.t$dropin) - 1))
  theils.pval <- pnorm(residuals, lower.tail = T)
  use <- names(theils.t$theils.t[which(theils.pval < 0.10)])

  # if(!supress.plot) {
  #
  # }

  # Calculate means that serve for binning the genes
  mu <- log2(rowMeans(data[use,]) + 1)
  bins <- split(mu, lsr::quantileCut(mu, 20))
  bin_mu <- sapply(bins, function(s) {
    in_bin <- names(s)
    mean(mu[in_bin])
  })

  # Calculate detection rate per bin
  l_detection_rate <- lapply(bins, function(s) {
    in_bin <- names(s)
    e <- data[in_bin,] > expression_cutoff
    colSums(e) / nrow(e)
  })

  detectionRate <- do.call("rbind", l_detection_rate)

  # Fit a GLM model per cell, capturing the relationship between detection rate and bin mean
  parameters <- t(apply(detectionRate, 2, function(dr) {
    mod <- glm(dr ~ bin_mu, family = gaussian(link="log"))
    coefs <- mod$coefficients
    names(coefs) <- c("Intercept", "Slope")
    coefs
  }))

  # Calculate weight per cell and gene, based on the GLMs parameters and the mean expression
  weights <- t(sapply(1:n, function(i) {
    mu <- means[i, condition]
    1 / (1 + exp(-(parameters[,1] - parameters[,2] * mu)))
  }))
  return(weights)
}

#' Calculates means per condition
#'
#' @param data Expression data matrix
#' @param condition A factor describing different conditions
#'
calcMeansPerCondition <- function(data, condition) {
  condition <- as.factor(condition)
  sapply(levels(condition), function(c) {
    i <- which(condition == c)
    apply(data[,i], 1, function(g) mean(singlecellutils::winsorize(g)))
  })
}
