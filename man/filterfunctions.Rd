% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter.R
\name{filterfunctions}
\alias{filterfunctions}
\alias{EpOverA}
\alias{FanoOverA}
\alias{MeanBelowA}
\alias{MeanOverA}
\title{Basic filter functions}
\usage{
EpOverA(A = 4, na.rm = TRUE)

FanoOverA(A = 2, na.rm = TRUE)

MeanBelowA(A = 4, na.rm = TRUE)

MeanOverA(A = 4, na.rm = TRUE)
}
\arguments{
\item{A}{The value that should / should not be exceeded.}

\item{na.rm}{Whether NAs should be removed.}
}
\value{
A function with bindings for A.
}
\description{
Basic filter functions return functions with bindings for A. The function evaluates to \code{TRUE} if some operation on the input results in a value greater/smaller than A.
}
\details{
The basic filter functions in detail:
\itemize{
\item \code{EpOverA} evaluates to \code{TRUE} if the aggregate of the arguments elements are larger than A.
\item \code{FanoOverA} evaluates to \code{TRUE} if the Fano factor of the argument is larger than A.
\item \code{MeanBelowA} evaluates to \code{TRUE} if the mean of the arguments elements are smaller than A.
\item \code{MeanOverA} evaluates to \code{TRUE} if the mean of the arguments elements are greater than A.
}
}
