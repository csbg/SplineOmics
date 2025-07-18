% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster_hits.R
\name{cluster_hits}
\alias{cluster_hits}
\title{cluster_hits()}
\usage{
cluster_hits(
  splineomics,
  nr_clusters,
  adj_pthresholds = c(0.05),
  adj_pthresh_avrg_diff_conditions = 0,
  adj_pthresh_interaction_condition_time = 0,
  genes = NULL,
  plot_info = list(y_axis_label = "Value", time_unit = "min", treatment_labels = NA,
    treatment_timepoints = NA),
  plot_options = list(cluster_heatmap_columns = FALSE, meta_replicate_column = NULL),
  raw_data = NULL,
  report_dir = here::here(),
  report = TRUE,
  max_hit_number = 25
)
}
\arguments{
\item{splineomics}{An S3 object of class `SplineOmics` that contains all the
necessary data and parameters for the analysis, including:
\itemize{
  \item \code{data}: The original expression dataset used for differential
  expression analysis.
  \item \code{meta}: A dataframe containing metadata corresponding to the
  \code{data}, must include a 'Time' column and any columns specified by
  \code{conditions}.
  \item \code{design}: A character of length 1 representing the limma
  design formula.
  \item \code{condition}: Character of length 1 specifying the column name
  in \code{meta} used to define groups for analysis.
  \item \code{spline_params}: A list of spline parameters for the analysis.
  \item \code{meta_batch_column}: A character string specifying the column
  name in the metadata used for batch effect removal.
  \item \code{meta_batch2_column}: A character string specifying the second
  column name in the metadata used for batch effect removal.
  \item \code{limma_splines_result}: A list of data frames, each representing
   a top table from differential expression analysis, containing at least
   'adj.P.Val' and expression data columns.
}}

\item{nr_clusters}{A list whose length matches `top_tables`; each element is
a numeric vector of positive integers (e.g. `1:1`, `2:8`) giving the
candidate number(s) of clusters for the corresponding condition level.}

\item{adj_pthresholds}{Numeric vector of p-value thresholds for filtering
hits in each top table.}

\item{adj_pthresh_avrg_diff_conditions}{p-value threshold for the results
from the average difference of the condition limma result. Per default 0 (
turned off).}

\item{adj_pthresh_interaction_condition_time}{p-value threshold for the
results from the interaction of condition and time limma result. Per default
0 (turned off).}

\item{genes}{A character vector containing the gene names of the features to
be analyzed. The entries should be in the same order as they appear in data.}

\item{plot_info}{List containing the elements y_axis_label (string),
time_unit (string), treatment_labels (character vector),
treatment_timepoints (integer vector). All can also be NA.
This list is used to add this info to the spline plots.
time_unit is used to label the x-axis, and treatment_labels
and -timepoints are used to create vertical dashed lines,
indicating the positions of the treatments (such as
feeding, temperature shift, etc.).}

\item{plot_options}{List with specific fields (cluster_heatmap_columns =
Bool) that allow for customization of plotting behavior.}

\item{raw_data}{Optional. Data matrix with the raw (unimputed) data, still 
containing NA values. When provided, it highlights the datapoints in the 
spline plots that originally where NA and that were imputed.}

\item{report_dir}{Character string specifying the directory path where the
HTML report and any other output files should be saved.}

\item{report}{Boolean TRUE or FALSE value specifing if a report should be
generated.}

\item{max_hit_number}{Maximum number of hits which are plotted within each
cluster. This can be used to limit the computation time and size of
the HTML report in the case of many hits. Default is 100.}
}
\value{
A list where each element corresponds to a group factor and contains
the clustering results,
        including `clustered_hits` data frame, hierarchical clustering object
        `hc`, `curve_values`
        data frame with normalized spline curves, and `top_table` with
        cluster assignments.
}
\description{
Performs clustering on hits from top tables generated by differential
expression analysis.
This function filters hits based on adjusted p-value thresholds, extracts
spline coefficients for
significant features, normalizes these coefficients, and applies hierarchical
clustering. The results,
including clustering assignments and normalized spline curves, are saved in a
specified directory and
compiled into an HTML report.
}
\seealso{
\code{\link[limma]{topTable}}, \code{\link[stats]{hclust}}
}
