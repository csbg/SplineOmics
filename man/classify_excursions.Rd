% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_pvc.R
\name{classify_excursions}
\alias{classify_excursions}
\title{Classify Peaks, Valleys, and Cliffs from Compound Contrast P-values}
\usage{
classify_excursions(data, meta, pvc_pvals, alpha = 0.05, symmetry_ratio = 0.3)
}
\arguments{
\item{data}{A numeric matrix or data frame of expression values, with 
rows as features and columns matching \code{meta$Time}.}

\item{meta}{A data frame containing metadata with at least a column 
\code{Time} indicating timepoint assignment.}

\item{alpha}{Significance threshold for p-value filtering. Only 
timepoints with \code{p < alpha} are considered for classification.}

\item{symmetry_ratio}{Minimum ratio of the smaller to larger slope 
required to consider a pattern symmetric enough to be a peak or valley 
(default: 0.3). Values below this threshold are classified as cliffs.}

\item{peak_valley_pvals}{A matrix of adjusted p-values from a compound 
contrast test, with rows matching \code{data} and columns corresponding 
to internal timepoints (i.e., excluding first and last timepoints).}
}
\value{
A data frame with the same number of rows as \code{data} and 
one column per timepoint. Each entry is a character label: \code{"p"}, 
\code{"v"}, \code{"t"}, \code{"b"}, or \code{"."} for no signal.
}
\description{
Assigns a label to each internal timepoint (T₂ to Tₙ₋₁) indicating 
whether it is a peak (\code{"p"}), valley (\code{"v"}), top of a 
cliff (\code{"t"}), or bottom of a cliff (\code{"b"}). Timepoints 
that are not statistically significant or cannot be classified are 
labeled as \code{"."}.
}
\details{
The classification is based on compound contrast p-values from the 
\code{peak_valley_test()} function, as well as the relative symmetry 
of changes between adjacent timepoints.
}
