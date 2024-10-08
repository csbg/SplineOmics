% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_limma_splines.R
\name{remove_intercept}
\alias{remove_intercept}
\title{Remove intercept from a formula}
\usage{
remove_intercept(formula)
}
\arguments{
\item{formula}{A formula object. The formula can include an intercept (`1`)
and other terms. If a `1` is found, it is replaced with `0`.}
}
\value{
A modified formula with the intercept removed. The first standalone
  occurrence of `1` will be replaced by `0`.
}
\description{
This function modifies a given formula by replacing the first occurrence
of a standalone intercept (`1`) with `0`. It works even if the `1` is
preceded by a tilde (`~`), ensuring that the intercept is removed while
leaving other parts of the formula intact.
}
