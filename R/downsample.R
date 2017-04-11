#' Thinning of a count matrix
#' @param y The count matrix
#' @param current.lib.sized The current library sizes
#' @param target.lib.sizes The target library sizes
#'
#' @return A list with two components:
#'   counts: The thinned count matrix
#'   lib.sizes: The new library sizes
#' @export
thin.counts = function(y, current.lib.sizes = colSums(y),
                       target.lib.sizes = min(current.lib.sizes)) {
  n = dim(y)[1];
  ncols = dim(y)[2];

  if (any (target.lib.sizes > current.lib.sizes)) {
    stop("Error: taget.lib.sizes > current.lib.sizes!");
  }

  if (length(target.lib.sizes) == 1) {
    target.lib.sizes = rep(target.lib.sizes, ncols);
  }

  y.thinned = y;

  for (k in 1:ncols) {
    keep.p =  target.lib.sizes[k]/current.lib.sizes[k]; # proportion of reads to keep
    if (keep.p == 1) {
      y.thinned[, k] <- y[,k]
    }
    else {
      y.thinned[,k] = rbinom(n, y[,k], keep.p);
    }
  }
  list(counts =y.thinned, lib.sizes = target.lib.sizes);
}

#' Downsampling of a count matrix to a target size
#'
#' @param counts The count matrix
#' @param totalReads The library size vector (in millions)
#' @param mean The mean target library size (in millions)
#' @param sd The standard deviation around the mean, used to create target library sizes
#' @param lowerDetectionLimit Number of counts considered as lower detection limit
#'
#' @return A list with the following components
#'   expressed:
downsample <- function(counts, totalReads, mean = 1, sd = 0.1, lowerDetectionLimit = 10) {
  target_size <- rnorm(ncol(counts), mean=mean, sd=sd)

  fraction <- abs(target_size) / totalReads
  fraction <- ifelse(fraction > 1, 1, fraction)

  new_total_reads <- totalReads * fraction
  target_libsizes <- colSums(counts)*fraction
  #plot(colSums(counts), target_libsizes, xlab="Total counts of core set", ylab="Downsampled counts of core set", pch=20, main=paste("Downsampled to mean", mean, "and sd",sd,sep=" "))

  ds_counts <- thin.counts(counts, current.lib.sizes = colSums(counts), target.lib.sizes = target_libsizes)
  ds_counts$counts[is.na(ds_counts$counts)] <- 0

  expressed <- ds_counts$counts > lowerDetectionLimit

  cell.order <- order(new_total_reads, decreasing = F)

  r <- rollapply(cell.order, width=20, function(ci) { mean(colSums(expressed[,ci]))}, fill="extend")

  r <- r / length(selected_features)

  list(expressed = r, fraction = fraction[cell.order], total = totalReads[cell.order], downsampled = new_total_reads[cell.order])
  # plot(1:length(cell.order), r, xlab="Window starting at cell", ylab="Percentage of mean detected features", pch=20, main=paste("Downsampled to mean", mean, "and sd",sd,sep=" "))
  #
  # m <- rollapply(cell.order, width=20, function(ci) {mean(new_total_reads[ci])}, fill="extend")
  # plot(m, r, xlab="Mean total sequenced reads", ylab="Percent of detected genes", pch=20, main=paste("Downsampled to mean", mean, "and sd",sd,sep=" "))
}
