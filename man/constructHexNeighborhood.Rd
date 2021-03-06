% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dimension_reduction.R
\name{constructHexNeighborhood}
\alias{constructHexNeighborhood}
\title{Construct the neighborhood of a hexagonal r*q matrix}
\usage{
constructHexNeighborhood(r, q, even.layout = T)
}
\arguments{
\item{r}{The number of rows}

\item{q}{The number of columns}

\item{even.layout}{Whether the layout is even}
}
\value{
A list of length n*q with their neighbors on the grid.
}
\description{
Construct the neighborhood of a hexagonal r*q matrix
}
