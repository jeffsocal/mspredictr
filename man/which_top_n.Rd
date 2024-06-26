% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extendr-wrappers.R
\name{which_top_n}
\alias{which_top_n}
\title{Returns the boolean index of the top N largest values}
\usage{
which_top_n(f, n)
}
\arguments{
\item{f}{A vector of numerical floats}

\item{n}{The number of top values to keep}
}
\description{
Returns the boolean index of the top N largest values
}
\examples{
which_top_n(c(123.45, 234.56, 345.67), 2)
}
