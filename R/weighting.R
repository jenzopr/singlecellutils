#' Calculates weights/probabilities based on false-negative curves.
#'
#' @param data Expression data matrix
#' @param means A vector containing mean expression values, e.g. from a bulk measurement
#' @param nq Number of quantiles used for binning the means
#' @param expression_cutoff A numerical value
#' @param supress.plot Whether to supress plotting (default TRUE)
#'
#' @return A matrix with probabilites that a gene i in cell j is not expressed, given that it it not detected P(not expressed|not detected). The values are only meaningful if the expression of gene i in cell j is below the detection limit.
#'
#' @export
calcFNWeight <- function(data, means = NULL, nq = 30, expression_cutoff = 10, supress.plot = T) {
  if (!is.numeric(data))
    stop("Argument data should be numeric.")
  data <- as.matrix(data)
  n <- nrow(data)
  m <- ncol(data)

  if (is.null(means)) {
    nz_data <- data
    nz_data[nz_data < 1] <- NA
    means <- rowMeans(log2(nz_data +1), na.rm = T)
  } else {
    means <- as.matrix(means)
  }

  # Either use all constitutively expressed genes or housekeeping genes here
  theils.t <- norm.theil.t(data, winsorize = T, treat.zeros = "eps")
  residuals <- scale(theils.t$theils.t - (exp(theils.t$i + theils.t$d * theils.t$dropin) - 1))
  theils.pval <- pnorm(residuals, lower.tail = T)
  use <- names(theils.t$theils.t[which(theils.pval < 0.2)])

  # Calculate means that serve for binning the genes
  # Use bulk means here?
  mu <- means[use]
  bins <- split(mu, lsr::quantileCut(mu, nq))
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

  dropoutRate <- 1 - do.call("rbind", l_detection_rate)

  # Fit a GLM model per cell, capturing the relationship between detection rate and bin mean
  parameters <- t(apply(dropoutRate, 2, function(do) {
    mod <- glm(do ~ bin_mu, family = binomial(link="logit"))
    coefs <- mod$coefficients
    names(coefs) <- c("Intercept", "Slope")
    coefs
  }))

  # Calculate P(not detected | expressed)
  p_nd_e <- t(sapply(means, sigmoid, a = parameters[,1], b = parameters[,2]))

  # Plot
  if (!supress.plot) {
    opar <- par()
    par(mfrow=c(4,4))
    for(i in 1:m) {
      plot(bin_mu, dropoutRate[,i], pch=16)
      lines(bin_mu, sigmoid(bin_mu, a = parameters[i,1], b = parameters[i,2]), type="l", col="red")
    }
    par(opar)
  }

  # Calculare priors P(detection) and P(expression)
  p_d <- rowSums(data > expression_cutoff)/m
  p_e <- p_d / rowMeans(1 - p_nd_e)

  # Calculate P(expressed | not detected)
  p_e_nd <- (p_nd_e * p_e) / (1 - p_d)

  # Calculate P(not expressed | not detected)
  p_ne_nd <- 1 - p_e_nd

  # Clip to [0,1]
  p_ne_nd[p_ne_nd < 0] <- 0

  return(p_ne_nd)
}

#' Returns the value of the sigmoid function
#' 1 - (1 / (1 + exp(a + (b*x)))
#'
#' @param x Parameter x of the function
#' @param a Parameter a of the function
#' @param b Parameter b of the function
#'
#' @return The evaluated function value.
#'
sigmoid <- function(x, a, b) {
  1-(1 / (1 + exp(a + (x * b))))
}

#' Calculates means per condition
#'
#' @param data Expression data matrix
#' @param condition A factor describing different conditions
#'
calcMeansPerCondition <- function(data, condition, log = T) {
  condition <- as.factor(condition)
  m <- sapply(levels(condition), function(c) {
    i <- which(condition == c)
    apply(data[,i], 1, function(g) mean(singlecellutils::winsorize(g)))
  })
  if (log) {
    return(log2(m+1))
  } else {
    return(m)
  }
}
