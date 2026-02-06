# Generate HTML report with p-value histograms of all the features

Generates an HTML report based on the results of a limma analysis with
splines. The report includes various plots and sections summarizing the
analysis results for time effects, average differences between
conditions, and interaction effects between condition and time.

## Usage

``` r
create_limma_report(
  splineomics,
  adj_pthresh = 0.05,
  report_dir = tempdir(),
  verbose = FALSE
)
```

## Arguments

- splineomics:

  `SplineOmics`: An S3 object of class `SplineOmics` that contains all
  the necessary data and parameters for the analysis, including:

  - `limma_splines_result`: `list` A list containing top tables from
    differential expression analysis for the three different limma
    results.

  - `meta`: `data.frame` A data frame with sample metadata. Must contain
    a column `"Time"`.

  - `condition`: `character(1)` A character string specifying the column
    name in the metadata (`meta`) that defines groups for analysis. This
    column contains levels such as `"exponential"` and `"stationary"`
    for phases, or `"drug"` and `"no_drug"` for treatments.

  - `annotation`: `data.frame` A data frame containing feature
    information, such as gene and protein names, associated with the
    expression data.

  - `report_info`: `list` A list containing metadata about the analysis
    for reporting purposes.

- adj_pthresh:

  `numeric(1)`: A numeric value specifying the adjusted p-value
  threshold for significance. Default is 0.05. Must be \> 0 and \< 1.

- report_dir:

  `character(1)`: A string specifying the directory where the report
  should be saved. Default is the current working directory.

- verbose:

  `logical(1)`: Boolean flag controlling the display of messages.

## Value

A list of plots included in the generated HTML report.

## Examples

``` r
set.seed(1)

# --- Toy data: 4 features x 6 samples ---
toy_data <- matrix(
    rnorm(4 * 6),
    nrow = 4, ncol = 6,
    dimnames = list(paste0("feat", 1:4), paste0("S", 1:6))
)

# --- Metadata with required columns (Time, Condition) ---
toy_meta <- data.frame(
    SampleID = colnames(toy_data),
    Time = c(0, 0, 1, 1, 2, 2),
    Condition = factor(c("Ctrl", "Ctrl", "Ctrl", "Trt", "Trt", "Trt"),
        levels = c("Ctrl", "Trt")
    ),
    row.names = colnames(toy_data),
    check.names = FALSE
)

# --- Minimal annotation (feature-level) ---
toy_anno <- data.frame(
    feature_id = rownames(toy_data),
    gene = paste0("G", 1:4),
    row.names = rownames(toy_data),
    check.names = FALSE
)

# --- Helper to fabricate limma-like topTables ---
make_tt <- function(n = 4) {
    p <- runif(n)
    ap <- p.adjust(p, method = "BH")
    data.frame(
        logFC = rnorm(n),
        AveExpr = rnorm(n, 5),
        t = rnorm(n),
        P.Value = p,
        adj.P.Val = ap,
        B = rnorm(n),
        row.names = paste0("feat", 1:n),
        check.names = FALSE
    )
}

# Structure expected by create_limma_report():
# time_effect: list of per-condition tables
# avrg_diff_conditions: list of per-contrast tables
# interaction_condition_time: list of per-contrast tables
toy_limma_res <- list(
   time_effect = list(
       "Condition_Ctrl" = make_tt(),
       "Condition_Trt"  = make_tt()
   ),
   avrg_diff_conditions = list(
       "avrg_diff_Ctrl_vs_Trt" = make_tt()
   ),
   interaction_condition_time = list(
       "time_interaction_Ctrl_vs_Trt" = make_tt()
   )
)

# Build SplineOmics object
so <- create_splineomics(
    data = toy_data,
    meta = toy_meta,
    condition = "Condition",
    annotation = toy_anno,
    design = "1 ~ Condition + Time",
    spline_params = list(spline_type = c("n", "n"), dof = c(2L, 2L)),
    report_info = list(
        omics_data_type       = "RNA-seq (toy)",
        data_description      = "Simulated expression matrix (4x6)",
        data_collection_date  = "2025-10-07",
        analyst_name          = "Analyst A",
        contact_info          = "analyst@example.org",
        project_name          = "SplineOmics Demo"
    )
)

# Attach limma results to the object
so <- update_splineomics(so, limma_splines_result = toy_limma_res)

# --- Generate the HTML report into a temporary directory ---
out_plots <- create_limma_report(
    splineomics = so,
    adj_pthresh = 0.05,
    report_dir  = tempdir(),
    verbose     = FALSE
)
```
