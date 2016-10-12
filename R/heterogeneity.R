#' Calculate the distance to the median
#'
#' @param data A numerical matrix or data.frame.
#' @param window The window size to calculate the summary statistic in.
#' @param normalize Whether or not to normalize to the standard deviation.
#' @param func The function to calculate the summary statistic within windows.
#'
#' @return A vector with the distances to the median.
#'
#' @export
dm <- function(data, ...) {
    ret <- normByWindow(data = data, ..., stat_fun = cv2.fun)
    ret
}

#' Calculate the Theil's T statistic for each gene and normalize it to the dropout rate. Also calculates a p-value for the residual deviation from a normal distribution.
#'
#' @param data A numerical matrix or data.frame.
#' @param use.qantile The quantile that has to be exceeded for fitting.
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
norm.theil.t <- function(data, use.quantile = 0.2, ...) {
  dropin <- 1 - apply(data, 1, function(x) sum(x == 0) / length(x))
  theil <- apply(data, 1, theils.t, ...)

  # Fit model
  use <- theil > quantile(theil, use.quantile)
  fit <- lm(log(theil[use] + 1)~dropin[use])
  i <- fit$coefficients[1]
  d <- fit$coefficients[2]

  theoretical <- exp(i + d * dropin) - 1
  res <- scale(theil - theoretical)
  pval <- pnorm(res, lower.tail = F)

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

#' Calculates the gini index
#'
#' @param data A numerical matrix or data.frame
#' @param window The window size to calculate the summary statistic in.
#' @param normalize Whether or not to normalize to the standard deviation.
#' @param func The function to calculate the summary statistic within windows.
#'
#' @return A vector with gini indices.
#'
#' @export
norm.gini <- function(data, winsorize = FALSE, ...) {
    if (winsorize) {
        gini.fun <- function(x) lawstat::gini.index(winsorize(x))$statistic
    } else {
        gini.fun <- function(x) lawstat::gini.index(x)$statistic
    }
    ret <- normByWindow(data = data, ..., stat_fun = gini.fun)
    ret
}

#' A generic function to normalize a statistic by a summary statistic across windows
#'
#' @param data A numerical matrix or data.frame.
#' @param stat_fun A function to calculate the statistic.
#' @param norm_to A numeric vector to normalize to.
#' @param func A function to summarize the statistic within a window.
#' @param window The window size to calculate the summary statistic in.
#' @param normalize Whether or not to normalize to the standard deviation.
#'
#' @return A vector with the summarized statistic
normByWindow <- function(data, stat_fun = NULL, norm_to = log2(rowMeans(data / 10) + 1), func = median, window = 100, normalize = TRUE) {
    if (is.null(rownames(data))) {
        rownames(data) <- 1:nrow(data)
    }

    # Assign names to norm_to
    names(norm_to) <- rownames(data)

    # Calculate the statistic with the supplied function
    stat <- apply(data, 1, stat_fun)

    # Sort the means
    s <- names(sort(norm_to))

    # Calculate summary statistic in a window of x genes
    r <- zoo::rollapply(data = stat[s], width = window, FUN = func, fill = "extend")
    r2 <- zoo::rollapply(data = stat[s], width = window, FUN = sd, fill = "extend")
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
    lim <- quantile(x, probs = c(fraction, 1 - fraction))
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
  m <- mean(x)
  n <- length(x)
  te <- x / m * log(x / m)
  1 / n * sum(te, na.rm = T)
}

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

cv2.fun <- function(x) {
  x <- x[!is.na(x)]
  if (mean(x) == 0) {
    return(0)
  }
  log2(var(x) / mean(x)^2 + 1)
}

#' Highly variable genes by Brennecke et. al.
#'
#' @param data A normalized count table.
#' @param min.cv2 Minimum coefficient of variation for genes used in fit.
#' @param mean.quantile The quantile that should be used to define a minimum mean for genes used in fit.
#' @param padj.method Method used in adjusting p-values.
#' @param windsorize Whether windsorization should be applied to data
#'
#' @return A list with the fit coefficients, the ordered variance and corresponding adjusted p-values.
#'
#' @export
hvg <- function(data, min.cv2 = .3, mean.quantile = .95, padj.method = "fdr") {
  means <- rowMeans(data)
  vars <- apply(data, 1, var)
  cv2 <- vars / means^2

  min.mean <- unname(quantile(means[which(cv2 > min.cv2)], mean.quantile))
  use <- means >= min.mean
  fit <- statmod::glmgam.fit(cbind( a0 = 1, a1tilde = 1 / means[use] ), cv2[use])
  a0 <- unname( fit$coefficients["a0"] )
  a1 <- unname( fit$coefficients["a1tilde"])
  afit <- a1 / means + a0
  varFitRatio <- vars / (afit * means^2)
  varorder <- order(varFitRatio, decreasing = T)
  df <- ncol(data) - 1
  pval <- pchisq(varFitRatio * df, df = df, lower.tail = F)
  adj.pval <- p.adjust(pval, method = padj.method)
  list(a0 = a0, a1 = a1, order = varorder, pval = adj.pval[varorder], df = df)
}
