#' Calculate the distance to the median
#'
#' @param data A numerical matrix or data.frame.
#' @param window The window size to calculate the summary statistic in.
#' @param normalize Whether or not to normalize to the standard deviation.
#' @param func The function to calculate the summary statistic within windows.
#'
#' @return A vector with the distances to the median.
#'
dm = function(data, window = 100, normalize = T, func = median) {
    if (is.null(rownames(expr))) {
        rownames(expr) <- 1:nrow(expr)
    }

    # Calculate log2 mean of the rows
    m <- log2(rowMeans(expr/10) + 1)
    names(m) <- rownames(expr)

    # Calculate the squared cv
    cv2 <- apply(expr, 1, calc_cv2)

    # Sort the means
    s <- names(sort(m))

    # Calculate median squared cv in a window of x genes
    r <- zoo::rollapply(data = log2(cv2[s] + 1), width = window, FUN = func, fill = "extend")
    r2 <- zoo::rollapply(data = log2(cv2[s] + 1), width = window, FUN = sd, fill = "extend")
    names(r) <- s
    names(r2) <- s

    if(normalize) {
      # Return standardized distance to the median
      ret <- (log2(cv2 + 1) - r[names(cv2)])/r2[names(cv2)]
    } else {
      ret <- log10(cv2 + 1) - r[names(cv2)]
    }

    ret[is.na(ret)] <- 0
    ret
}


