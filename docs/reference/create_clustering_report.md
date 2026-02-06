# Generate an HTML clustering report from precomputed results

Creates an HTML report and associated visualizations for clustering
results obtained from a prior SplineOmics analysis. This function
performs no statistical modeling or clustering itself; instead, it
consumes a precomputed *report payload* and is solely responsible for
report generation, plotting, and writing outputs to disk.

## Usage

``` r
create_clustering_report(
  report_payload,
  plot_info = list(y_axis_label = "Value", time_unit = "min", treatment_labels = NA,
    treatment_timepoints = NA),
  plot_options = list(cluster_heatmap_columns = FALSE, meta_replicate_column = NULL,
    condition_colours = NULL),
  verbose = FALSE,
  max_hit_number = 25,
  raw_data = NULL,
  report_dir = tempdir()
)
```

## Arguments

- report_payload:

  `list`: A structured report payload object as returned by upstream
  analysis and clustering functions (e.g.,
  [`cluster_hits()`](https://csbg.github.io/SplineOmics/reference/cluster_hits.md)).
  The payload must contain all data and metadata required to generate
  the report, including (non-exhaustive):

  - `splineomics`: A `SplineOmics` object containing the original input
    data, metadata, design, spline parameters, and analysis settings.

  - `all_levels_clustering`: Clustering assignments for all condition
    levels and analysis categories.

  - `genes`: [`character()`](https://rdrr.io/r/base/character.html)
    Vector of gene or feature names aligned to the rows of the input
    data.

  - `predicted_timecurves`: Predicted spline trajectories on a common
    time grid.

  - `category_2_and_3_hits`: Data describing significant features from
    average-difference and interaction analyses.

  - Adjusted p-value thresholds and additional metadata required for
    report annotation.

- plot_info:

  `list`: List with optional elements used to annotate spline plots:

  - `y_axis_label`: `character(1)` Label for the y-axis.

  - `time_unit`: `character(1)` Unit used in the x-axis label.

  - `treatment_labels`: `list(character(1))` Named list of labels
    describing experimental interventions.

  - `treatment_timepoints`: `list(numeric(1))` Named list of time points
    at which interventions occur.

  If any treatment list is supplied, both must be present and must share
  identical names corresponding to levels of `meta[[condition]]`.
  Vertical dashed lines are drawn at the specified time points in the
  spline plots and annotated with the provided labels.

- plot_options:

  `list`: Named list controlling optional plot customization. Supported
  entries include:

  - `cluster_heatmap_columns` `logical(1)` (default = `FALSE`): Whether
    to cluster columns in the heatmap.

        \item \code{meta_replicate_column} `character(1)`
          (default = \code{NULL}):
          Column in \code{meta} encoding replicate information. When supplied,
          spline plot data points are colored by replicate.

        \item \code{condition_colours} `list`
          (default = \code{NULL}):
          Optional named list mapping condition levels to colours. Names must
          correspond to values in the condition column of \code{meta}, and
          values must be valid R colour specifications accepted by
          \code{ggplot2}. When provided, these colours override the default
          colours for the corresponding conditions.

- verbose:

  `logical(1)`: Logical flag controlling the verbosity of status and
  progress messages during report generation.

- max_hit_number:

  `integer(1)`: Maximum number of hits plotted per cluster. This
  parameter can be used to limit report size and computation time when
  many significant features are present.

- raw_data:

  `matrix`: Optional matrix containing the raw (unimputed) data with
  missing values still present. When supplied, data points that were
  originally missing and subsequently imputed are highlighted in the
  spline plots.

- report_dir:

  `character(1)`: Path to the directory where the HTML report and all
  associated output files will be written. The directory is created if
  it does not already exist.

## Value

A named list of plots corresponding to the visualizations included in
the generated HTML report. This includes:

- Cluster-specific spline plots for time effects.

- Heatmaps summarizing cluster assignments.

- Plots for average-difference and interaction (condition–time)
  clustering results, when available.

The returned plots are invisibly written to disk as part of the report
generation process.

## Details

This strict separation between computation and reporting ensures that
analytical results can be inspected, tested, and reused independently of
visualization or file-system side effects.
