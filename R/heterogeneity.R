#' Adds a heterogeneity statistic to the object.
#'
#' @param object A SingleCellExperiment object.
#' @param exprs_values String indicating which assay contains the data that should be used calculate heterogeneous genes.
#' @param column Determines the column name of the \code{\link{rowData}} slot to store results to.
#' @param ... Additional parameters passed to \code{\link{heterogeneity}}
#'
#' @return A SingleCellExperiment object with modified \code{rowData} slot.
#'
#' @export
add_heterogeneity <- function(object, exprs_values = "counts", column = ".heterogeneity", ...) {
  het <- heterogeneity(SummarizedExperiment::assay(object, i = exprs_values), ...)
  object <- scater::mutate(object, !!column := het)
  return(object)
}

#' Calculates the coefficient of variation
#'
#' @param x The vector of values
#'
#' @return The coefficient of variation for x
#'
#' @export
cv2.fun <- function(x) {
  x <- x[!is.na(x)]
  if (mean(x) == 0) {
    return(0)
  }
  stats::var(x) / mean(x)^2
}

#' Calculates the drouput percentage
#'
#' @param x The vector of values
#' @param limit The detection limit
#' @return The percentage of dropouts
#'
#' @export
dropout.fun <- function(x, limit=0) {
  x <- x[!is.na(x)]
  sum(x<=limit)/length(x)
}

#' Calculates the gini index
#'
#' @param x The vector of values
#' @return The gini index for x
#'
#' @export
gini.fun <- function(x) {
  x <- x[!is.na(x)]
  #lawstat::gini.index(x)$statistic
  ox = x[order(x)]
  n = length(x)
  dsum = drop(crossprod(2 * 1:n - n - 1, ox))
  dsum / (mean(x) * n * (n - 1))
}

#' Calculates the inequality of a numeric vector
#'
#' @param x The numeric vector of values
#' @return The inequality measure for x
#'
#' @export
inequality.fun <- function(x) {
  x <- x[!is.na(x)]
  Y <- sum(x, na.rm = T)
  return(sum(x / Y * log( x/Y * length(x)), na.rm = T))
}

#' Calculates the inequality between and within groups
#'
#' @param x The numeric vector of values
#' @param g The grouping factor
#' @return The inequality measure for x given g
#' @export
bw.inequality.fun <- function(x, g) {
  sum(theils.between(x, g)) / sum(theils.within(x, g))
}

#' Calculates various heterogeneity statistics
#'
#' @param data A numerical matrix
#' @param statistic The heterogeneity measure to calculate. Can be \code{cv} for the coefficient of variation, \code{dropout} for the percentage of dropouts, \code{gini} for the gini index or \code{inequality} for Theil's T statistic.
#' @param normalization The normalization method. Can be \code{none} for no normalization, \code{bins} to normalize in bins calculated from the mean or \code{windows} to normalize within rolling windows.
#' @param order_by A numeric vector by which the calculated statistic should be ordered in case of normalization.
#' @param groups An optional grouping factor.
#' @param ... Additional arguments passed to the normalization functions.
#'
#' @return The statistics value for each observation.
#'
#' @export
heterogeneity <- function(data, statistic = c("cv", "dropout", "gini", "inequality", "bwt", "mean"), normalization = c("none", "bins", "windows"), order_by = log2(Matrix::rowMeans(data / 10, na.rm = T) + 1), groups = NULL, ...) {
  # Transform the data based on the statistic given
  stat <- match.arg(statistic)
  transformed_data <- switch(stat,
                             cv = apply(data, 1, cv2.fun),
                             dropout = apply(data, 1, dropout.fun),
                             mean = log(rowMeans(data))/log(10),
                             gini = apply(data, 1, gini.fun),
                             inequality = apply(data, 1, inequality.fun),
                             bwt = apply(data, 1, bw.inequality.fun, g = groups))

  # Normalize data
  norm <- match.arg(normalization)
  ret <- switch(norm,
                none = transformed_data,
                bins = normByBin(transformed_data, order_by = order_by, ...),
                windows = normByWindow(transformed_data, order_by = order_by, ...))
  ret
}

#' A generic function to normalize a statistic by binning another statistic
#'
#' @param stat A numeric vector.
#' @param order_by A numeric vector to bin by.
#' @param func A function to summarize the statistic within a bin.
#' @param n_bins The number of bins to create.
#' @param normalize Whether or not to normalize to the standard deviation.
#'
#' @return A vector with the difference to the summarized statistic
#'
#' @export
normByBin <- function(stat, order_by = stat, func = mean, n_bins = 30, normalize = TRUE) {
  if (is.null(names(stat))) {
    names(stat) <- 1:length(stat)
  }

  bins <- split(stat, lsr::quantileCut(order_by, n_bins))
  bin_summary <- sapply(bins, function(s) {
    in_bin <- names(s)
    func(stat[in_bin])
  })

  if (normalize) {
    bin_sd <- sapply(bins, function(s) {
      in_bin <- names(s)
      stats::sd(stat[in_bin])
    })
  } else {
    bin_sd <- rep(1, n_bins)
  }

  ret_l <- sapply(1:n_bins, function(i) {
    in_bin <- names(bins[[i]])
    (stat[in_bin] - bin_summary[i]) / bin_sd[i]
  })
  unlist(ret_l)
}

#' A generic function to normalize a statistic by a summary statistic across windows
#'
#' @param stat A numeric vector.
#' @param order_by A numeric vector to order stat by.
#' @param func A function to summarize the statistic within a window.
#' @param window The window size to calculate the summary statistic in.
#' @param normalize Whether or not to normalize to the standard deviation.
#'
#' @return A vector with the difference to the summarized statistic
#'
#' @export
normByWindow <- function(stat, order_by = stat, func = stats::median, window = 100, normalize = TRUE) {
  if (is.null(names(stat))) {
    names(stat) <- 1:length(stat)
  }

  # Assign names to norm_to
  names(order_by) <- names(stat)

  # Sort the means
  s <- names(sort(order_by))

  # Calculate summary statistic in a window of x genes
  r <- zoo::rollapply(data = stat[s], width = window, FUN = func, fill = "extend")
  r2 <- zoo::rollapply(data = stat[s], width = window, FUN = stats::sd, fill = "extend")
  names(r) <- s
  names(r2) <- s

  if (normalize) {
    # Return standardized distance to the summary statistic
    ret <- (stat - r[names(stat)]) / r2[names(stat)]
  } else {
    ret <- stat - r[names(stat)]
  }

  ret[is.na(ret)] <- 0
  ret
}

#' Calculate the Theil's T statistic for each gene and normalize it to the dropout rate. Also calculates a p-value for the residual deviation from a normal distribution.
#'
#' @param data A numerical matrix or data.frame.
#' @param use.quantile The quantile that has to be exceeded for fitting.
#' @param ... other parameters passed to theils.t
#'
#' @return A list with the following components:
#' \itemize{
#'  \item d The linear models coeffictient,
#'  \item i The linear models intercept,
#'  \item theils.t A named vector containing Theils T statistic,
#'  \item dropin A named vector with 1 - dropout values,
#'  \item pval The p-value for the residuals deviation from the normal distribution.
#'  }
#'
#' @export
norm.theil.t <- function(data, use.quantile = 0.05, ...) {
  dropin <- 1 - apply(data, 1, function(x) sum(x == 0) / length(x))
  theil <- apply(data, 1, inequality.fun)#, ...)

  # Fit model
  use <- theil > stats::quantile(theil, use.quantile)
  fit <- stats::lm(log(theil[use] + 1)~dropin[use])
  i <- fit$coefficients[1]
  d <- fit$coefficients[2]

  theoretical <- exp(i + d * dropin) - 1
  res <- scale(theil - theoretical)
  pval <- stats::pnorm(res, lower.tail = F)

  # Order by pval and assemble return
  o <- order(pval)
  return(list(
    d = d,
    i = i,
    theils.t = theil[o],
    dropin = dropin[o],
    pval = pval[o]
  ))
}

#' Winsorizes a vector to either a fraction or a absolute number
#'
#' @param x A numeric vector
#' @param fraction The fraction to calculate the limits
#' @param absolute The absolute number to calculate the limits
#' @param two.sided If winsorization should be done two-sided
#'
#' @return The winsorized vector
#'
#' @export
winsorize <- function(x, fraction = 0.05, absolute = NULL, two.sided = TRUE) {
    if (!is.null(absolute)) {
        fraction <- absolute / length(x)
    }
    if (length(fraction) != 1 || fraction < 0 || fraction > 0.5) {
        stop("bad value for 'fraction'")
    }
    lim <- stats::quantile(x, probs = c(fraction, 1 - fraction))
    if (two.sided) {
        x[x < lim[1]] <- lim[1]
    }
    x[x > lim[2]] <- lim[2]
    x
}

#' Calculate Theil's T-Statistic for a vector
#'
#' @param x A numeric vector
#' @param winsorize If winsorization should be applied to x
#' @param treat.zeros How zeros should be treated. One of c("exclude","epsilon")
#' @param epsilon In case of treat.zeros equals to "epsilon", a small number is added to observations of x equal to zero.
#' @param ... Parameters passed to winsorize, if applicable.
#'
#' @return Theil's T statistic for x
#'
#' @export
theils.t <- function(x, winsorize = FALSE, treat.zeros = "exclude", epsilon = 0, ...) {
  z <- match.arg(treat.zeros, c("exclude", "epsilon"))
  if (winsorize) {
    x <- winsorize(x, ...)
  }
  x <- switch (z,
    exclude = x[!(x == 0)],
    epsilon = x + epsilon
  )
  m <- mean(x, na.rm = T)
  n <- length(x)
  te <- x / m * log(x / m)
  1 / n * sum(te, na.rm = T)
}


#' Calculates the Theils index between groups.
#'
#' @param x A numeric vector
#' @param groups A grouping factor.
#'
#' @return The Theils between.
#' @export
theils.between <- function(x, groups) {
  P <- length(groups)
  Y <- sum(x)
  bgte <- sapply(unique(groups), function(g) {
    ind <- (groups == g)
    y <- sum(x[ind])
    p <- sum(ind)
    return(y / Y * log( (y / Y) / (p / P)))
  })

  #T <- sum(bgte, na.rm = T)
  T <- bgte
  return(T)
}

#' Calculates the Theils index within groups.
#'
#' @param x A numeric vector
#' @param groups A grouping factor.
#'
#' @return The Theils within.
#' @export
theils.within <- function(x, groups) {
  ug <- unique(groups)
  Y <- sum(x)
  Yg <- sapply(ug, function(g) sum(x[groups == g]))
  wgte <- sapply(1:length(ug), function(gi) {
    ind <- (groups == ug[gi])
    sum((x[ind] / Yg[gi]) * log((x[ind] / Yg[gi]) / (1 / sum(ind))), na.rm = T)
  })

  #T <- sum(Yg/Y * wgte)
  T <- Yg / Y * wgte
  return(T)
}

#' Highly variable genes by Brennecke et. al.
#'
#' @param data A normalized count table.
#' @param min.cv2 Minimum coefficient of variation for genes used in fit.
#' @param mean.quantile The quantile that should be used to define a minimum mean for genes used in fit.
#' @param padj.method Method used in adjusting p-values.
#'
#' @return A list with the fit coefficients, the ordered variance and corresponding adjusted p-values.
#'
#' @export
hvg <- function(data, min.cv2 = .3, mean.quantile = .95, padj.method = "fdr") {
  means <- rowMeans(data)
  vars <- apply(data, 1, stats::var)
  cv2 <- vars / means^2

  min.mean <- unname(stats::quantile(means[which(cv2 > min.cv2)], mean.quantile))
  use <- means >= min.mean
  fit <- statmod::glmgam.fit(cbind( a0 = 1, a1tilde = 1 / means[use] ), cv2[use])
  a0 <- unname( fit$coefficients["a0"] )
  a1 <- unname( fit$coefficients["a1tilde"])
  afit <- a1 / means + a0
  varFitRatio <- vars / (afit * means^2)
  varorder <- order(varFitRatio, decreasing = T)
  df <- ncol(data) - 1
  pval <- stats::pchisq(varFitRatio * df, df = df, lower.tail = F)
  adj.pval <- stats::p.adjust(pval, method = padj.method)
  list(a0 = a0, a1 = a1, order = varorder, pval = adj.pval[varorder], df = df)
}

