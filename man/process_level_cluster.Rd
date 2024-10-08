% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster_hits.R
\name{process_level_cluster}
\alias{process_level_cluster}
\title{Process Level Cluster}
\usage{
process_level_cluster(
  top_table,
  cluster_size,
  level,
  meta,
  condition,
  spline_params,
  mode
)
}
\arguments{
\item{top_table}{A dataframe containing the top table results from limma.}

\item{cluster_size}{The size of clusters to generate.}

\item{level}{The level within the condition to process.}

\item{meta}{A dataframe containing metadata.}

\item{condition}{A character string specifying the condition.}

\item{spline_params}{A list of spline parameters for the analysis.}

\item{mode}{A character string specifying the mode
('isolated' or 'integrated').}
}
\value{
A list containing the clustering results, including curve values and
        the design matrix.
}
\description{
Processes clustering for a specific level within a condition using the
provided top table and spline parameters.
}
\seealso{
\code{\link{get_curve_values}}, \code{\link{normalize_curves}},
\code{\link{hierarchical_clustering}}
}
