% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dimension_reduction.R
\name{umap}
\alias{umap}
\title{Conveinience function to perform UMAP via reticulate.}
\usage{
umap(object, exprs_values, features = NULL, scale = T,
  n_neighbors = 5L, min_dist = 0.1, metric = "correlation",
  use_dimred = NULL, seed = NULL)
}
\arguments{
\item{object}{A SingleCellExperiment object.}

\item{exprs_values}{String indicating which assay contains the data that should be used to perform UMAP.}

\item{features}{A character vector (of feature names), a logical vector or numeric vector (of indices) specifying the features to use for UMAP. The default of NULL will use all features.}

\item{scale}{Logical scalar indicating whether data should be scaled before UMAP.}

\item{n_neighbors}{This determines the number of neighboring points used in local approximations of manifold structure. Larger values will result in more global structure being preserved at the loss of detailed local structure. In general this parameter should often be in the range 5 to 50, with a choice of 10 to 15 being a sensible default.}

\item{min_dist}{This controls how tightly the embedding is allowed compress points together. Larger values ensure embedded points are more evenly distributed, while smaller values allow the algorithm to optimise more accurately with regard to local structure. Sensible values are in the range 0.001 to 0.5, with 0.1 being a reasonable default.}

\item{metric}{This determines the choice of metric used to measure distance in the input space.}

\item{use_dimred}{String indicating which slot in \code{reducedDim} contains the data that should be used to perform UMAP (if \code{exprs_values = NULL}).}

\item{seed}{A numeric seed to initialize the random number generator.}
}
\value{
A matrix with the two-dimensional embedding.
}
\description{
Performs Uniform Manifold Approximation and Projection (UMAP) on a SingleCellExperiment object.
The algorithm is explained in detail in \href{https://arxiv.org/abs/1802.03426}{McInnes et. al.} and requires the python packages \href{https://pypi.org/project/umap-learn/}{umap-learn} and \href{https://pypi.org/project/numpy/}{numpy} to be installed.
}
