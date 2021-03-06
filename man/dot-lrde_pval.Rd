% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/differential_tests.R
\name{.lrde_pval}
\alias{.lrde_pval}
\title{Function to perform logistic regression test}
\usage{
.lrde_pval(data, response, response.ref = NULL)
}
\arguments{
\item{data}{Matrix with the predictor variables in columns}

\item{response}{Vector with two levels indicating the response (class) variable.}

\item{response.ref}{The level of \code{response} that should be used as reference.

This function needs the package \code{lmtest} to compare GLMs.
In case no response reference level is set, but the response contains more than two levels, a warning is reported and the first level is automatically chosen.}
}
\value{
A p-value of the \code{\link[lmtest]{lrtest}} between the full and null model.
}
\description{
Function to perform logistic regression test
}
