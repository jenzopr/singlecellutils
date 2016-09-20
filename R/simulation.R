#' Simulates an single-cell RNA-seq experiment
#'
#' @param genes Number of genes in the experiment
#' @param cells Number of cells in the experiment
#' @param beta.shape1,beta.shape2 Shape parameter for a beta distribution, from which burst size and burst frequency are estimated
#' @param logis.loc Location parameter of the bivariate logistic distribution
#' @param logis.scale Scale parameter of the bivariate logistic distribution
#'
#' @return An FPM expression matrix
#'
#' @export
simulatescSeq <- function(genes=5500, cells=150, beta.shape1 = 3, beta.shape2 = 5, logis.loc = 1.2, logis.scale = 1.5) {
  burst.size <- 10*rbeta(genes, beta.shape1, beta.shape2)
  burst.freq <- 10*rbeta(genes, beta.shape1, beta.shape2)
  expression.scale <- (burst.size*burst.freq)^2

  expr_mat <- sapply(1:n, function(x) {
    mean <- expression.scale[x] * rbeta(1, burst.size[x], 10-burst.freq[x])
    base <- rpois(N, lambda = mean)
    rate <- calc_dropout(burst.size[x],burst.freq[x], location = logis.loc, scale = logis.scale)
    i <- sample(1:N, size=ceiling(N*rate))
    base[i] = 0
    return(base)
  })

  return(expr_mat)
}

#' Caclulates drop-outs based on burst size and frequency
#'
#' @param size The genes burst size
#' @param frequency The genes burst frequency
#' @param location The location parameter of the bivariate logistic distribution
#' @param scale The scale parameter of the bivariate logistic distribution
#'
#' @return The gene-specific drop-out rate
calc_dropout <- function(size, frequency, location = 1.2, scale = 1.5) {
  scale.density <- VGAM::dbilogis(location, location, loc1 = location, loc2 = location, scale1 = scale, scale2 = scale)
  VGAM::dbilogis(size, frequency, loc1 = location, loc2 = location, scale1 = scale, scale2 = scale)/scale.density
}




