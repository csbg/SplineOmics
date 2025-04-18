% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_limma_splines.R
\name{build_contrasts_for_pair}
\alias{build_contrasts_for_pair}
\title{Construct contrast names for a pair of condition levels}
\usage{
build_contrasts_for_pair(condition, level_pair, design_matrix)
}
\arguments{
\item{condition}{A string giving the name of the condition factor.}

\item{level_pair}{A character vector of length 2 specifying the two levels
to be compared.}

\item{design_matrix}{The design matrix used for model fitting, whose column
names determine which condition level is modeled explicitly.}
}
\value{
A named list with two elements:
  - `condition`: name of the coefficient for the condition-only effect
  - `interaction`: vector of coefficient names for the interaction terms
}
\description{
Builds the coefficient names used for extracting the condition-only and
condition-time interaction effects between two condition levels. Assumes
a design matrix where one condition level is used as the reference.
}
