% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_gsea.R
\name{check_genes}
\alias{check_genes}
\title{Check Valid Gene IDs}
\usage{
check_genes(genes, max_index_overall = NA)
}
\arguments{
\item{genes}{A character vector containing gene IDs.}

\item{max_index_overall}{An integer, specifying the highest index of all
features across all levels.}
}
\value{
None. This function stops execution and provides
an error message if the vector does not meet the criteria,
including the first offending element and its index.
}
\description{
This function checks whether a character vector `genes`
contains only valid gene IDs. Each gene ID must consist
solely of alphabetic letters and numbers.
}
