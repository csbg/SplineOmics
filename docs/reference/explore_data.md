# Automatically perform exploratory data analysis (EDA)

This function automatically generates exploratory data analysis (EDA)
plots from the provided data. These include density plots, boxplots, PCA
plots, MDS plots, variance explained plots, violin plots, mean
correlation with time, first lag autocorrelation, lag1 differences, and
coefficient of variation. The function returns all EDA plots as a list
and, by default, creates an HTML report containing the plots, saving it
to the specified report directory.

## Usage

``` r
explore_data(splineomics, report_dir = tempdir(), report = TRUE)
```

## Arguments

- splineomics:

  `SplineOmics`: A named list containing all required inputs for the
  splineomics workflow. Must contain the following elements:

  - `data`: `matrix` The data matrix with the values. The columns are
    the samples (timepoint + replicate combo) and the rows are the
    features (e.g. genes or proteins).

  - `meta`: `data.frame` A dataframe containing metadata corresponding
    to the `data`. Must include a 'Time' column and any columns
    specified by `conditions`. In general, the columns of meta
    correspond to different metadata types, and each row corresponds to
    a column of data (metadata for that sample).

  - `annotation`: `data.frame` A dataframe that maps the rows of `data`
    to annotation info, such as gene names or database identifiers.

  - `report_info`: `list` A named list describing the experiment. Must
    include the following fields (`character(1)`): -
    `"omics_data_type"` - `"data_description"` -
    `"data_collection_date"` - `"analyst_name"` - `"contact_info"` -
    `"project_name"`

    May also include the following optional fields (`character(1)`): -
    `"method_description"` - `"results_summary"` - `"conclusions"`

  - `condition`: `character(1)` Character vector of length 1 specifying
    the column name in `meta` used to define groups for analysis.

  - `meta_batch_column`: `character(1)` Character vector of length 1
    specifying the column name in `meta` that contains the info for the
    batch effect.

  - `meta_batch2_column`: `character(1)` Character vector of length 1
    specifying the column name in `meta` that contains the info for the
    second batch effect.

- report_dir:

  `character(1)`: A non-empty string specifying the report directory.
  The output HTML reports will be placed there.

- report:

  `logical(1)`: A Boolean TRUE or FALSE value, specifying if a report
  should be generated or not.

## Value

A list of ggplot objects representing various exploratory plots.

## Examples

``` r
# Temporary output dir for the HTML report
report_dir <- file.path(tempdir(), "explore_data_demo")
if (!dir.exists(report_dir)) dir.create(report_dir, recursive = TRUE)

## --- Load example inputs from inst/extdata ---
data <- readRDS(xzfile(system.file(
    "extdata", "proteomics_data.rds.xz",
    package = "SplineOmics"
)))

meta <- read.csv(
    system.file("extdata", "proteomics_meta.csv", package = "SplineOmics"),
    stringsAsFactors = FALSE
)

# Derive the annotation table from the data frame (as in the vignette)
first_na_col <- which(is.na(data[1, ]))[1]
annotation <- data |>
    dplyr::select((first_na_col + 1):ncol(data)) |>
    dplyr::slice(-c(1:3))

data <- SplineOmics::extract_data(
    # The dataframe with the numbers on the left and info on the right.
    data = data,
    # Use this annotation column for the feature names.
    feature_name_columns = c("Gene_name"),
    top_row = 4,
    bottom_row = 1165,
    right_col = 37,
    left_col = 2
)

# Minimal report metadata
report_info <- list(
    omics_data_type = "PTX",
    data_description = "Demo proteomics dataset from extdata",
    data_collection_date = "2024",
    analyst_name = "Example Analyst",
    contact_info = "analyst@example.org",
    project_name = "DemoProject"
)

# Build the SplineOmics input object
splineomics <- SplineOmics::create_splineomics(
    data               = data,
    meta               = meta,
    annotation         = annotation,
    report_info        = report_info,
    condition          = "Phase", # column in `meta` defining groups
    meta_batch_column  = "Reactor" # optional: for batch-effect removal
)

# Run EDA and write the HTML report to the temp dir
plots <- SplineOmics::explore_data(
    splineomics = splineomics,
    report_dir  = report_dir,
    report      = TRUE
)
#> Making density plots...
#> Making violin plots...
#> Making PCA plots...
#> Making MDS plots...
#> Making correlation heatmaps...
#> Subsampled to top 1000 most variable features (after filtering rows with > 20% missing) for correlation heatmap.
#> Making mean correlation with time plots...
#> Making lag1 differences plots...
#> Making first lag auto-correlation with time plots...
#> Making cv plots...
#> Making density plots...
#> Making violin plots...
#> Making PCA plots...
#> Making MDS plots...
#> Making correlation heatmaps...
#> Subsampled to top 1000 most variable features (after filtering rows with > 20% missing) for correlation heatmap.
#> Making mean correlation with time plots...
#> Making lag1 differences plots...
#> Making first lag auto-correlation with time plots...
#> Making cv plots...
#> 
#>  Info Exploratory data analysis completed successfully.
#>  Your HTML reports are located in the directory:  /tmp/RtmpyRGyv7/explore_data_demo .
#>  Please note that due to embedded files, the reports might be flagged as
#>  harmful by other software. Rest assured that they provide no harm.

# Inspect what was written
list.files(report_dir, recursive = TRUE)
#> [1] "explore_batch_corrected_data_PTX_05_02_2026-12_38_55.html"
#> [2] "explore_data_PTX_05_02_2026-12_38_55.html"                

# `plots` is a named list of ggplot objects (e.g., plots$raw_data$pca, etc.)
# print(plots$raw_data$pca)
```
