% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster_hits.R
\name{check_between_level_pattern}
\alias{check_between_level_pattern}
\title{Check for Between-Level Patterns in Top Tables}
\usage{
check_between_level_pattern(top_tables)
}
\arguments{
\item{top_tables}{A list where each element is itself a list containing
named elements.}
}
\value{
A list with two elements:
\describe{
  \item{between_levels}{A logical value indicating whether any element names
  match the between-level pattern.}
  \item{index_with_pattern}{The index of the first element in `top_tables`
  where all names match the between-level pattern, or NA if no match is
  found.}
}
}
\description{
This function checks if any of the elements within a list of top tables
contain element names that match the specified between-level pattern.
}
\details{
The function iterates over each element in `top_tables`. For each element
that
is a list, it checks if all names within that inner list match the pattern
`".+_vs_.+"`. If a match is found, the function sets `between_levels` to TRUE
and records the index of the matching element. The search stops at the first
match.
}
