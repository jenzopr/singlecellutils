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

#' Calculate the normalized Theil's T statistic
#'
#' @param data A numerical matrix or data.frame.
#' @param window The window size to calculate the summary statistic in.
#' @param normalize Whether or not to normalize to the standard deviation.
#' @param func The function to calculate the summary statistic within windows.
#'
#' @return A vector with the normalize Theil's T statistic.
#'
#' @export
norm.theil.t <- function(data, ...) {
  dropout <- apply(data,1,function(x) sum(x==0)/length(x))
  ret <- normByWindow(data = data, norm_to = dropout, ..., stat_fun = theils.t)
  ret
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
#' @param func A function to summarize the statistic within a window.
#' @param window The window size to calculate the summary statistic in.
#' @param normalize Whether or not to normalize to the standard deviation.
#'
#' @return A vector with the summarized statistic
normByWindow <- function(data, stat_fun = NULL, norm_to = log2(rowMeans(data/10) + 1), func = median, window = 100, normalize = TRUE) {
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
        ret <- (stat - r[names(stat)])/r2[names(stat)]
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
        fraction <- absolute/length(x)
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

cv2.fun <- function(x) {
  x <- x[!is.na(x)]
  if (mean(x) == 0) {
    return(0)
  }
  log2(var(x)/mean(x)^2 + 1)
}

theils.t <- function(x) {
  n <- length(x)
  m <- mean(x) # maybe winsorize
  te <- x/m * log(x/m)
  te[is.nan(te)] <- 0
  1/n*sum(te)
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
  cv2 <- vars/means^2

  min.mean <- unname(quantile(means[which(cv2 > min.cv2)], mean.quantile))
  use <- means >= min.mean
  fit <- statmod::glmgam.fit(cbind( a0 = 1, a1tilde = 1/means[use] ), cv2[use])
  a0 <- unname( fit$coefficients["a0"] )
  a1 <- unname( fit$coefficients["a1tilde"])
  afit <- a1/means+a0
  varFitRatio <- vars/(afit*means^2)
  varorder <- order(varFitRatio,decreasing=T)
  df <- ncol(data) - 1
  pval <- pchisq(varFitRatio*df,df=df,lower.tail=F)
  adj.pval <- p.adjust(pval,method = padj.method)
  list(a0 = a0, a1 = a1, order = varorder, pval = adj.pval[varorder], df = df)
}
