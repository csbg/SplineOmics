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

## Examples

``` r
## Runnable example using packaged extdata (small subset for speed)
data_in <- readRDS(
  xzfile(
    system.file(
      "extdata",
      "phosphoproteomics_data.rds.xz",
      package = "SplineOmics"
    )
  )
)

meta <- read.csv(
  system.file(
    "extdata",
    "phosphoproteomics_meta.csv",
    package = "SplineOmics"
  ),
  stringsAsFactors = FALSE
)

first_na_col <- which(is.na(data_in[1, ]))[1]
annotation_all <- data_in[, (first_na_col + 1):ncol(data_in), drop = FALSE]

## Extract a small feature subset to keep examples fast
data_mat <- SplineOmics::extract_data(
  data = data_in,
  feature_name_columns = c("T..Gene"),
  use_row_index = TRUE,
  top_row = 1,
  bottom_row = 50,
  right_col = 36,
  left_col = 1
)

## IMPORTANT: subset annotation to match the selected features
keep_ids <- rownames(data_mat)
annotation <- annotation_all[
    annotation_all$T..Gene %in% keep_ids, , drop = FALSE
    ]
annotation <- annotation[
    match(keep_ids, annotation$T..Gene), , drop = FALSE
    ]

splineomics <- SplineOmics::create_splineomics(
  data = data_mat,
  meta = meta,
  annotation = annotation,
  condition = "Phase",
  meta_batch_column = "Reactor"
)

pvc_results <- SplineOmics::find_pvc(
  splineomics = splineomics,
  alphas = 0.025,
  padjust_method = "BH",
  support = 1
)
#> design matrix of interest not specified. Assuming a one-group experiment.
#> Warning: Partial NA coefficients for 2 probe(s)
#> Warning: Partial NA coefficients for 12 probe(s)
#> design matrix of interest not specified. Assuming a one-group experiment.
#> Warning: Partial NA coefficients for 2 probe(s)
#> Warning: Partial NA coefficients for 9 probe(s)

report_dir <- file.path(tempdir(), "pvc_report_example")
dir.create(report_dir, showWarnings = FALSE, recursive = TRUE)

out <- SplineOmics::create_pvc_report(
  splineomics = splineomics,
  pvc_results = pvc_results,
  verbose = FALSE,
  report_dir = report_dir
)

names(out)
#> [1] "Exponential" "Stationary" 
```
