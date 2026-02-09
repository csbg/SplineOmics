# Perform a limma spline analysis to identify significant features

This is the core function, which performs a limma spline analysis to
identify significant time-dependent changes in features (e.g., proteins)
within an omics time-series dataset. It evaluates features within each
condition level and between levels by comparing average differences and
interactions between time and condition.

## Usage

``` r
run_limma_splines(splineomics, verbose = FALSE)
```

## Arguments

- splineomics:

  `SplineOmics`: An S3 object of class `SplineOmics` that contains the
  following elements:

  - `data`: `matrix` The matrix of the omics dataset, with the feature
    names optionally as row headers.

  - `rna_seq_data`: `list` An optional object containing the
    preprocessed RNA-seq data, such as the output from
    [`limma::voom`](https://rdrr.io/pkg/limma/man/voom.html) or a
    similar preprocessing pipeline. This must only be provided when the
    input is RNA-seq data.

  - `meta`: `data.frame` A dataframe containing metadata corresponding
    to the `data`, must include a 'Time' column and the column specified
    by `condition`. The columns of meta describe the meta info, such as
    the time and condition, and each row corresponds to a column in
    data, and therefore, contains the meta info for that data column. It
    is important that meta and data are matched in this way.

  - `padjust_method`: `character(1)` Statistical method that is used for
    multiple hypothesis correction. Supported methods include all that
    are included in the p.adjust() function in R: "holm", "hochberg",
    "hommel", "bonferroni", "BH" or "fdr", "BY", or "none" for no
    correction. Default of this package is "BH".

  - `design`: `character(1)` A character string representing the limma
    design formula, such as "~ 1 + Phase\*Time + Reactor" for an
    integrated design, or "~ 1 + Time + Reactor" for an isolated design.

  - `dream_params`: `list | NULL` A named list or NULL. When not NULL,
    it can contain the following named elements: - `dof`: `integer(1)`
    An integer greater than 1, specifying the degrees of freedom for the
    dream topTable. If set to 0, then the best dof is automatically
    found with the help of leave-one-out-crossvalidation (loocv). The
    dof with the lowest error on the loocv is chosen. - `KenwardRoger`:
    `logical(1)` A boolean indicating whether to use the Kenward-Roger
    approximation for mixed models. Note that random effects are now
    directly specified in the design formula and not in `dream_params`.

  - `mode`: `character(1)` Specifies how the design formula is
    constructed: either `"isolated"` or `"integrated"`.

    - `"isolated"`: Each level is analyzed independently, using only the
      subset of data corresponding to that level. The design formula
      does not include the condition variable, since only one condition
      is present in each subset.

    - `"integrated"`: All levels are analyzed together in a single
      model, using the full dataset. The design formula includes the
      condition variable (and optionally interaction terms with it) so
      that results are estimated jointly across all levels.

    - `condition`: `character(1)` A character string specifying the
      column name in `meta` used to define groups for analysis. The
      condition column contains the levels of the experiment (such as
      control and treatment).

    - `spline_params`: `list` A list of spline parameters used in the
      analysis, including:

      - `spline_type`: `character(1)`The type of spline (e.g., `"n"` for
        natural splines or `"b"` for B-splines).

      - `dof`: `integer(1)` Degrees of freedom for the spline. If set to
        `0`, `SplineOmics` automatically determines the optimal number
        of degrees of freedom using leave-one-out cross-validation and
        selects the value that yields the best predictive performance.

      - `degree`: `integer(1)` Polynomial degree of the spline basis
        (for B-splines only).

    - `use_array_weights`: `logical(1)` Boolean flag indicating if the
      robust fit strategy to deal with heteroscedasticity should be used
      or not. If set to NULL, then this is handeled implicitly based on
      the result of the Levene test. If this test is significant for at
      least 10% of the features, then the robust strategy is used. The
      robust strategy uses the function voomWithQualityWeights for
      RNA-seq data instead of the normal voom function. For other,
      non-count-based data, the function limma::arrayWeights is used
      instead, combined with setting the robust argument to TRUE for the
      limma::eBayes function. In summary, the strategy employed by those
      functions is to downweights samples with higher variance. This can
      be neccessary, because linear models have the assumption of
      homoscedasticity, which means that the variance is (approx.) the
      same across all datapoints where the linear model is fitted. If
      this is violated, then the resulting p-values cannot be trusted
      (common statistical wisdom).

    - `bp_cfg`: [`numeric()`](https://rdrr.io/r/base/numeric.html) A
      named numeric vector specifying the parallelization configuration,
      with expected names `"n_cores"`: `integer(1)` and
      `"blas_threads"`: `integer(1)`.

    This controls how many R worker processes (`n_cores`) and how many
    BLAS/OpenBLAS threads per process (`blas_threads`) should be used
    during parallel computation.

    If `bp_cfg` is `NULL`, missing, or any of its required fields is
    `NA`, both `n_cores` and `blas_threads` default to `1`. This
    effectively disables parallelization and avoids oversubscription of
    CPU threads.

- verbose:

  `logical(1)`: Boolean flag controlling the display of messages.

## Value

The SplineOmics object, updated with a list with three elements: -
`time_effect`: `list` A list of top tables for each level with the time
effect. - `avrg_diff_conditions`: `list` A list of top tables for each
comparison between the levels. The comparison is the average difference
of the values. - `interaction_condition_time`: `list` A list of top
tables for each comparison between levels. The comparison is the
interaction between the condition and the time.

## Examples

``` r
# Toy data: 4 features x 6 samples (two conditions, three time points)
toy_data <- matrix(
    c(
        3, 5, 8, 12, 17, 23, # f1
        23, 17, 13, 9, 6, 4, # f2
        5, 3, 2, 2, 3, 5, # f3
        1, 4, 9, 8, 4, 1, # f4
        10, 10, 10, 10, 10, 10, # f5
        2, 2, 2, 9, 12, 15, # f6
        4, 5, 7, 10, 14, 19, # f7
        12, 11, 9, 8, 9, 12 # f8
    ),
    nrow = 8, ncol = 6, byrow = TRUE,
    dimnames = list(paste0("f", 1:8), paste0("s", 1:6))
)

toy_meta <- data.frame(
    Time = c(0, 1, 2, 0, 1, 2),
    condition = rep(c("WT", "KO"), each = 3),
    Replicate = rep(c("R1", "R2"), each = 3),
    row.names = colnames(toy_data),
    stringsAsFactors = FALSE
)

toy_annot <- data.frame(
    feature_nr = 1:8,
    gene = c("G1", "G2", "G3", "G4"),
    stringsAsFactors = FALSE
)

# Stub limma "top tables" with minimal required fields
# (feature_nr + adj.P.Val)
tt_wt <- data.frame(feature_nr = 1:4, adj.P.Val = c(0.01, 0.20, 0.04, 0.60))
tt_ko <- data.frame(feature_nr = 1:4, adj.P.Val = c(0.50, 0.03, 0.70, 0.02))
tt_c2 <- data.frame(feature_nr = 1:4, adj.P.Val = c(0.04, 0.70, 0.80, 0.90))
tt_c3 <- data.frame(feature_nr = 1:4, adj.P.Val = c(0.20, 0.90, 0.03, 0.80))

design_str <- "~ 1 + Time*condition"

# Minimal spline parameters required by spline machinery
spline_params <- list(
    spline_type = "n", # natural cubic splines
    dof = 1L # degrees of freedom for the spline basis
)

toy_splineomics <- list(
    data = toy_data,
    meta = toy_meta,
    annotation = toy_annot,
    report_info = list(
        omics_data_type = "RNA-seq",
        data_description = "toy example",
        data_collection_date = "2025-01-01",
        analyst_name = "Example",
        contact_info = "example@example.org",
        project_name = "ToyProject"
    ),
    design = design_str,
    mode = "integrated",
    condition = "condition",
    spline_params = spline_params,
    meta_batch_column = NULL,
    meta_batch2_column = NULL,
    limma_splines_result = list(
        time_effect                  = list(WT = tt_wt, KO = tt_ko),
        avrg_diff_conditions         = tt_c2,
        interaction_condition_time   = tt_c3
    ),
    feature_name_columns = "gene"
)
class(toy_splineomics) <- "SplineOmics"

toy_splineomics <- run_limma_splines(toy_splineomics)
#> Running Levene's test both feature wise and sample wise to implicitly
#> • decide whether to use the limma array weights or not.
#> No random effects: fitting model with lmFit()...
#> Running feature wise Levene's test...
#> 
#> ------------------------------------------------------------
#> Fraction of features violating homoscedasticity
#>     (p < 0.050): 0.00% (0/8 features)
#> No violating features found.
#> ------------------------------------------------------------
#> Running Levene's test across samples to detect
#> • inter-sample variance differences...
#> Levene's test p-value (sample-level): 0.1425
#> ✅ No strong evidence of inter-sample variance differences.
#> ------------------------------------------------------------
#> ✅ No strong evidence for heteroscedasticity.
#> • Proceeding WITHOUT using robust strategy.
#> ------------------------------------------------------------
```
