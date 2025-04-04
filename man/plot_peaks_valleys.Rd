% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/detect_excursions_playground.R
\name{plot_peaks_valleys}
\alias{plot_peaks_valleys}
\title{Plot Excursions in Time-Series Omics Data}
\usage{
plot_peaks_valleys(results, data, meta, meta_replicates_column)
}
\arguments{
\item{results}{A list returned from `detect_excursions()`, containing 
`results_df` (excursion matrix) and `pairwise_pvals` (one-tailed p-values).}

\item{data}{A numeric matrix, where rows correspond to features (e.g., genes, 
proteins, metabolites) and columns correspond to samples.}

\item{meta}{A data frame containing metadata for the samples. Must include
a column named `"Time"` that specifies the timepoint for each sample.}

\item{meta_replicates_column}{A character string specifying the column name 
in `meta` that indicates biological replicates.}
}
\value{
A named list of ggplot objects, where each element corresponds to a 
feature with detected excursions. Each plot displays the expression levels 
across timepoints, with replicates distinguished by different shapes.
}
\description{
This function generates scatter plots for features that exhibit excursions
in time-series omics data. Excursion points are highlighted in red, while 
normal points remain grey. Significance stars are added to indicate the 
statistical significance of adjacent pairwise comparisons.
}
\details{
- The function first extracts features with at least one excursion.
- Each feature's expression is plotted across time using `ggplot2`.
- Replicates are displayed with distinct shapes, while excursion points 
  are highlighted in red.
- Significance stars (`*`, `**`, `***`, `****`) are placed between 
  adjacent timepoints if their pairwise comparison is significant.
- The significance stars are positioned above horizontal lines, with 
  right-side comparisons always placed slightly higher for clarity.
}
