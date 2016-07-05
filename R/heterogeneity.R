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
    cv2.fun = function(x) {
      x <- x[!is.na(x)]
      if(mean(x) == 0) {
        return(0)
      }
      log2((sd(x)/mean(x))^2+1)
    }
    ret <- normByWindow(data = data, ..., stat_fun = cv2.fun)
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
  if(winsorize) {
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
normByWindow <- function(data, stat_fun = NULL, func = median, window = 100, normalize = TRUE) {
  if (is.null(rownames(data))) {
    rownames(data) <- 1:nrow(data)
  }

  # Calculate log2 mean of the rows
  m <- log2(rowMeans(data/10) + 1)
  names(m) <- rownames(data)

  # Calculate the statistic with the supplied function
  stat <- apply(data, 1, stat_fun)

  # Sort the means
  s <- names(sort(m))

  # Calculate summary statistic in a window of x genes
  r <- zoo::rollapply(data = stat, width = window, FUN = func, fill = "extend")
  r2 <- zoo::rollapply(data = stat, width = window, FUN = sd, fill = "extend")
  names(r) <- s
  names(r2) <- s

  if(normalize) {
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
winsorize <- function (x, fraction = 0.05, absolute = NULL, two.sided = TRUE) {
  if(!is.null(absolute)) {
    fraction <- absolute/length(x)
  }
  if(length(fraction) != 1 || fraction < 0 || fraction > 0.5) {
    stop("bad value for 'fraction'")
  }
  lim <- quantile(x, probs=c(fraction, 1-fraction))
  if(two.sided) {
    x[ x < lim[1] ] <- lim[1]
  }
  x[ x > lim[2] ] <- lim[2]
  x
}
