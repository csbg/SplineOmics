# Create an object containing variables often used by SplineOmics functions

Creates a SplineOmics object containing variables that are commonly used
across multiple functions in the package. This object is then passed as
an argument to the other functions of this package.

## Usage

``` r
create_splineomics(
  data,
  meta,
  condition,
  rna_seq_data = NULL,
  annotation = NULL,
  report_info = NULL,
  meta_batch_column = NULL,
  meta_batch2_column = NULL,
  feature_name_columns = NULL,
  design = NULL,
  use_array_weights = FALSE,
  dream_params = NULL,
  mode = "isolated",
  spline_params = NULL,
  padjust_method = "BH",
  bp_cfg = NULL
)
```

## Arguments

- data:

  `matrix`: The actual omics data. If the `rna_seq_data` argument is
  used, still provide this argument. In that case, input the data matrix
  here (for example, the `$E` part of a `voom` object). Assign your
  feature names as row headers; otherwise, numeric indices will be used.

- meta:

  `data.frame`: Metadata associated with the omics data.

- condition:

  `character(1)`: Condition variable describing the experimental groups.

- rna_seq_data:

  `list` \| `NULL`: An object containing the preprocessed RNA-seq data,
  such as the output from
  [`limma::voom`](https://rdrr.io/pkg/limma/man/voom.html) or a similar
  pipeline. This argument is not validated directly; input checks rely
  on [`limma::lmFit()`](https://rdrr.io/pkg/limma/man/lmFit.html).

- annotation:

  `data.frame` \| `NULL`: Feature annotations (optional) providing
  descriptive information about each feature in `data`.

- report_info:

  `list`: Named list describing the experiment. Must include the
  following fields (all `character(1)`):

  - `"omics_data_type"`

  - `"data_description"`

  - `"data_collection_date"`

  - `"analyst_name"`

  - `"contact_info"`

  - `"project_name"`

  Optional fields (all `character(1)`):

  - `"method_description"`

  - `"results_summary"`

  - `"conclusions"`

- meta_batch_column:

  `character(1)` \| `NULL`: Column name in `meta` specifying batch
  information (optional).

- meta_batch2_column:

  `character(1)` \| `NULL`: Column name in `meta` specifying secondary
  batch information (optional).

- feature_name_columns:

  [`character()`](https://rdrr.io/r/base/character.html): Vector of
  column names from `annotation` that describe the features. Used in the
  HTML report to define how feature names displayed above each spline
  plot were created. Use the same vector that was used to create the row
  headers for the data matrix.

- design:

  `matrix` \| `NULL`: Design matrix or similar object (optional).

- use_array_weights:

  `logical(1)`: Boolean flag indicating whether to use the robust
  fitting strategy to handle heteroskedasticity. If `NULL`, this is
  determined automatically via the Levene test: if at least 10% of
  features are significant, the robust strategy is enabled. For RNA-seq
  data, this uses
  [`limma::voomWithQualityWeights()`](https://rdrr.io/pkg/limma/man/voomWithQualityWeights.html),
  otherwise
  [`limma::arrayWeights()`](https://rdrr.io/pkg/limma/man/arrayWeights.html)
  with `robust = TRUE` in
  [`limma::eBayes()`](https://rdrr.io/pkg/limma/man/ebayes.html). These
  approaches down-weight samples with higher variance, improving
  validity of statistical inference.

- dream_params:

  `list` \| `NULL`: Optional named list controlling mixed-model fitting.
  When not `NULL`, may include:

  - `dof` `integer(1)` Degrees of freedom for the DREAM `topTable`.

  - `KenwardRoger` `logical(1)` Whether to use the Kenward-Roger
    correction.

  Random effects are specified directly in the design formula, not here.

- mode:

  `character(1)`: Either `"isolated"` or `"integrated"`. Determines
  whether conditions are analysed independently (`"isolated"`) or
  jointly (`"integrated"`). The integrated mode fits a single model
  across all levels.

- spline_params:

  `list` \| `NULL`: Parameters for spline functions. Must contain:

  - `spline_type`: `character(1)` `"n"` for natural cubic or `"b"` for
    B-splines.

  - `dof`: `integer(1)` Degrees of freedom. If set to `0`, `SplineOmics`
    automatically determines the best value using leave-one-out
    cross-validation.

  - `degree`: `integer(1)` Degree of the spline (B-splines only).

- padjust_method:

  `character(1)`: Method for p-value adjustment. One of `"none"`,
  `"BH"`, `"BY"`, `"holm"`, `"bonferroni"`, `"hochberg"`, or `"hommel"`.
  Defaults to `"BH"` (Benjamini–Hochberg).

- bp_cfg:

  [`numeric()`](https://rdrr.io/r/base/numeric.html) \| `NULL`: Named
  numeric vector specifying parallelization settings, with expected
  names `"n_cores"` and `"blas_threads"`. Controls the number of R
  worker processes (`n_cores`) and BLAS/OpenBLAS threads per process
  (`blas_threads`). If `bp_cfg` is `NULL` or missing, both default to
  `1`, disabling parallelization and avoiding thread oversubscription.

## Value

A SplineOmics object.

## Examples

``` r
set.seed(1)

# 6 samples, 4 features
toy_data <- matrix(
    rnorm(4 * 6, mean = 0, sd = 1),
    nrow = 4, ncol = 6,
    dimnames = list(
        paste0("gene", 1:4),
        paste0("S", 1:6)
    )
)

# Sample metadata
toy_meta <- data.frame(
    SampleID = colnames(toy_data),
    Time = c(0, 0, 1, 1, 2, 2),
    Condition = factor(c("Ctrl", "Ctrl", "Ctrl", "Trt", "Trt", "Trt"),
        levels = c("Ctrl", "Trt")
    ),
    Batch = factor(c("B1", "B1", "B1", "B2", "B2", "B2")),
    stringsAsFactors = FALSE,
    row.names = colnames(toy_data)
)

# Condition vector (must align with samples)
cond <- toy_meta$Condition

# Minimal annotation (feature-level info)
toy_anno <- data.frame(
    feature_id = rownames(toy_data),
    symbol = c("G1", "G2", "G3", "G4"),
    stringsAsFactors = FALSE,
    row.names = rownames(toy_data)
)

# Spline parameters (natural splines with df = 3)
toy_spline <- list(spline_type = "n", dof = 3)

# Parallel config (single-threaded for examples)
toy_bp <- c(n_cores = 1, blas_threads = 1)

# Dream params example (optional)
toy_dream <- list(dof = 3L, KenwardRoger = FALSE)

# Simple design matrix (intercept + condition + time)
toy_design <- stats::model.matrix(~ Condition + Time, data = toy_meta)

# Required report fields
toy_report <- list(
    omics_data_type = "RNA-seq (toy)",
    data_description = "Simulated expression matrix (4x6)",
    data_collection_date = "2025-10-07",
    analyst_name = "Analyst A",
    contact_info = "analyst@example.org",
    project_name = "SplineOmics Demo",
    method_description = "Toy example to construct a SplineOmics object"
)

so <- create_splineomics(
    data                 = toy_data,
    meta                 = toy_meta,
    condition            = cond,
    rna_seq_data         = NULL, # not used in this toy
    annotation           = toy_anno,
    report_info          = toy_report,
    meta_batch_column    = "Batch",
    meta_batch2_column   = NULL,
    feature_name_columns = c("feature_id", "symbol"),
    design               = toy_design,
    use_array_weights    = FALSE,
    dream_params         = toy_dream,
    mode                 = "isolated",
    spline_params        = toy_spline,
    padjust_method       = "BH",
    bp_cfg               = toy_bp
)

class(so)
#> [1] "SplineOmics"
str(so, max.level = 1)
#> List of 16
#>  $ data                : num [1:4, 1:6] -0.626 0.184 -0.836 1.595 0.33 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ rna_seq_data        : NULL
#>  $ meta                :'data.frame':    6 obs. of  4 variables:
#>  $ condition           : Factor w/ 2 levels "Ctrl","Trt": 1 1 1 2 2 2
#>  $ annotation          :'data.frame':    4 obs. of  2 variables:
#>  $ report_info         :List of 7
#>  $ meta_batch_column   : chr "Batch"
#>  $ meta_batch2_column  : NULL
#>  $ feature_name_columns: chr [1:2] "feature_id" "symbol"
#>  $ design              : num [1:6, 1:3] 1 1 1 1 1 1 0 0 0 1 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   ..- attr(*, "assign")= int [1:3] 0 1 2
#>   ..- attr(*, "contrasts")=List of 1
#>  $ use_array_weights   : logi FALSE
#>  $ dream_params        :List of 2
#>  $ mode                : chr "isolated"
#>  $ spline_params       :List of 2
#>  $ padjust_method      : chr "BH"
#>  $ bp_cfg              : Named num [1:2] 1 1
#>   ..- attr(*, "names")= chr [1:2] "n_cores" "blas_threads"
#>  - attr(*, "class")= chr "SplineOmics"
```
