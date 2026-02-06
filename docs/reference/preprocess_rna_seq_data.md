# Perform essential preprocessing steps for raw RNA-seq counts

The `preprocess_rna_seq_data()` function performs essential
preprocessing steps for raw RNA-seq counts. This includes creating a
`DGEList` object, normalizing the counts using the default TMM (Trimmed
Mean of M-values) normalization via the
[`edgeR::calcNormFactors`](https://rdrr.io/pkg/edgeR/man/calcNormFactors.html)
function, and applying the `voom` transformation from the `limma`
package to obtain log-transformed counts per million (logCPM) with
associated precision weights. If you require a different normalization
method, you can supply your own custom normalization function.

## Usage

``` r
preprocess_rna_seq_data(splineomics, normalize_func = NULL, verbose = FALSE)
```

## Arguments

- splineomics:

  `SplineOmics`: An S3 object of class `SplineOmics` that must contain
  the following elements:

  - `data`: `matrix` A matrix of the omics dataset, with feature names
    optionally as row headers (genes as rows, samples as columns).

  - `meta`: `data.frame` A dataframe containing metadata corresponding
    to the `data`. The dataframe must include a 'Time' column and a
    column specified by the `condition`.

  - `design`: `character(1)` A character string representing the design
    formula for the limma analysis (e.g., `'~ 1 + Phase*X + Reactor'`).

  - `spline_params`: `list` A list of spline parameters used in the
    analysis. This can include:

    - `spline_type`: `character(1)` A character string specifying the
      type of spline. Must be either `'n'` for natural cubic splines or
      `'b'` for B-splines.

    - `dof`: `integer(1)` An integer specifying the degrees of freedom.
      Required for both natural cubic splines and B-splines.

    - `degree`: `integer(1)` An integer specifying the degree of the
      spline (for B-splines only).

  - `use_array_weights`: `logical(1)` \| `NULL` Boolean flag indicating
    if the robust fit strategy to deal with heteroscedasticity should be
    used or not. If set to NULL, then this is handled implicitly based
    on the result of the Levene test. If this test is significant for at
    least 10% of the features, then the robust strategy is used. The
    robust strategy uses the function voomWithQualityWeights for RNA-seq
    data instead of the normal voom function. For other, non-count-based
    data, the function limma::arrayWeights is used instead, combined
    with setting the robust argument to TRUE for the limma::eBayes
    function. In summary, the strategy employed by those functions is to
    downweight samples with higher variance. This can be necessary,
    because linear models have the assumption of homoscedasticity, which
    means that the variance is (approx.) the same across all datapoints
    where the linear model is fitted. If this is violated, then the
    resulting p-values cannot be trusted (common statistical wisdom).

  - `bp_cfg`: [`numeric()`](https://rdrr.io/r/base/numeric.html) \|
    `NULL` A named numeric vector specifying the parallelization
    configuration, with expected names `"n_cores"` and `"blas_threads"`.

    This controls how many **R worker processes** (`n_cores`) and how
    many **BLAS/OpenBLAS threads per process** (`blas_threads`) should
    be used during parallel computation.

    If `bp_cfg` is `NULL`, missing, or any of its required fields is
    `NA`, both `n_cores` and `blas_threads` default to `1`. This
    effectively disables parallelization and avoids oversubscription of
    CPU threads.

- normalize_func:

  `function` \| `NULL`: An optional normalization function. If provided,
  this function will be used to normalize the `DGEList` object. If not
  provided, TMM normalization (via
  [`edgeR::calcNormFactors`](https://rdrr.io/pkg/edgeR/man/calcNormFactors.html))
  will be used by default. Must take as input the y of: y \<-
  edgeR::DGEList(counts = raw_counts) and output the y with the
  normalized counts.

- verbose:

  `logical(1)`: Boolean flag controlling the display of messages.

## Value

The updaed `splineomics` object, now containing the `voom` object, which
includes the log2-counts per million (logCPM) matrix and
observation-specific weights. Additionally, the splineparams are updated
with the identified optimal dof based on LOOCV, when dof = 0L (for
auto-dof)

## See also

edgeR::DGEList, edgeR::calcNormFactors

## Examples

``` r
if (requireNamespace("edgeR", quietly = TRUE) &&
    requireNamespace("limma", quietly = TRUE)) {
    # Toy raw counts: 4 genes x 6 samples (2 conditions x 3 timepoints)
    set.seed(1)
    counts <- matrix(
        rpois(24, lambda = 20),
        nrow = 4, ncol = 6,
        dimnames = list(paste0("g", 1:4), paste0("s", 1:6))
    )

    meta <- data.frame(
        Time = c(0, 1, 2, 0, 1, 2),
        condition = rep(c("WT", "KO"), each = 3),
        row.names = colnames(counts)
    )

    # Simple fixed-effects design (no random effects → uses limma::voom)
    design_str <- "~ 0 + Time + condition"

    # Spline params: natural cubic with fixed dof (avoid auto-dof path)
    sp <- list(
        spline_type = "n", dof = 3L, degree = NA_integer_
    )

    # Build minimal SplineOmics object
    so <- list(
        data = counts,
        meta = meta,
        design = design_str,
        condition = "condition",
        spline_params = sp,
        use_array_weights = FALSE, # skip homoscedasticity auto-check
        bp_cfg = list(n_cores = 1L, blas_threads = 1L)
    )
    class(so) <- "SplineOmics"

    # Run preprocessing (TMM + voom)
    so2 <- preprocess_rna_seq_data(so)

    # Inspect the voom results (logCPM matrix dimensions)
    if (!is.null(so2$rna_seq_data)) {
        dim(so2$rna_seq_data$E)
    }
}
#> Coefficients not estimable: X3 
#> Warning: Partial NA coefficients for 4 probe(s)
#> [1] 4 6
```
