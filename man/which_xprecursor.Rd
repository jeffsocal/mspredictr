% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extendr-wrappers.R
\name{which_xprecursor}
\alias{which_xprecursor}
\title{Returns the boolean index of values that are in proximity to mz}
\usage{
which_xprecursor(f, mz)
}
\arguments{
\item{f}{A vector of numerical floats}

\item{mz}{The the value to look for}
}
\description{
Returns the boolean index of values that are in proximity to mz
}
\examples{
which_xprecursor(c(123.45, 234.56, 345.67), 233.61)
}
