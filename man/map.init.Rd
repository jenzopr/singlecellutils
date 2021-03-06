% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dimension_reduction.R
\name{map.init}
\alias{map.init}
\title{Initializes a SOM map with the first two PCA components in a linear fashion}
\usage{
map.init(data, resolution = 1)
}
\arguments{
\item{data}{The training data of the SOM}

\item{resolution}{A size factor increasing or decreasing the map resolution}
}
\value{
A list with initialization parameters
}
\description{
Initializes a SOM map with the first two PCA components in a linear fashion
}
