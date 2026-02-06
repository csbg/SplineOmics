# Create PVC HTML report from statistical PVC results

Generates per-feature PVC plots and writes an HTML report to disk from
the statistical output of
[`find_pvc()`](https://csbg.github.io/SplineOmics/reference/find_pvc.md).
This function does not recompute any statistics; it consumes
`pvc_results` (including its attributes) and the original input data
stored in `splineomics`.

## Usage

``` r
create_pvc_report(
  splineomics,
  pvc_results,
  plot_info = list(y_axis_label = "Value", time_unit = "min", treatment_labels = NA,
    treatment_timepoints = NA),
  verbose = FALSE,
  report_dir = tempdir()
)
```

## Arguments

- splineomics:

  `list`: Preprocessed time-series input. Must include at least `data`,
  `meta`, `annotation`, `condition`, `meta_batch_column`,
  `meta_batch2_column` (optional), `report_info`, and
  `feature_name_columns`.

- pvc_results:

  `list`: Output from
  [`find_pvc()`](https://csbg.github.io/SplineOmics/reference/find_pvc.md),
  a named list by condition level. Each level must contain `alpha` and
  `pvc_adj_pvals`. Attributes `padjust_method` and `support` are used
  for report header settings.

- plot_info:

  `list`: Plot annotation options passed to `plot_pvc()`. See
  [`find_pvc()`](https://csbg.github.io/SplineOmics/reference/find_pvc.md)
  for the expected structure. Treatment annotations are applied per
  condition level if `treatment_labels` and `treatment_timepoints` are
  provided.

- verbose:

  Boolean flag indicating if messages should be shown.

- report_dir:

  `character(1)`: Output directory for the HTML report and any
  associated files.

## Value

The input `pvc_results` with an additional `plots` element added under
each condition level. The report is written to `report_dir`.

## Details

For each condition level, the function subsets `splineomics$data` and
`splineomics$meta`, generates per-feature plots via `plot_pvc()`, then
calls `generate_report_html()` with the plots, metadata, and report
header settings derived from `pvc_results` and `splineomics`.
