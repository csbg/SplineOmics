# screen_limma_hyperparams()

This function screens through various combinations of hyperparameters
for limma analysis, including designs, modes, and degrees of freedom. It
validates inputs, generates results for all combinations, and plots the
outcomes. Finally, it may also be involved in generating an HTML report
as part of a larger analysis workflow.

## Usage

``` r
screen_limma_hyperparams(
  splineomics,
  datas,
  datas_descr,
  metas,
  designs,
  modes,
  spline_test_configs,
  report_dir = here::here(),
  adj_pthresholds = c(0.05),
  rna_seq_datas = NULL,
  time_unit = "min",
  padjust_method = "BH"
)
```

## Arguments

- splineomics:

  An S3 object of class \`SplineOmics\` that contains all the necessary
  data and parameters for the analysis, including:

  - `condition`: A string specifying the column name of the meta
    dataframe, that contains the levels that separate the experiment
    ('treatment' can be a condition, and 'drug' and 'no drug' can be the
    levels of such a condition).

  - `report_info`:

  - `meta_batch_column`: A character string specifying the meta batch
    column.

  - `meta_batch2_column`: A character string specifying the second meta
    batch column (the limma function removeBatchEffect supports a
    maximum of two batch columns.)

- datas:

  A list of matrices containing the datasets to be analyzed.

- datas_descr:

  A description object for the data.

- metas:

  A list of data frames containing metadata for each dataset in
  \`datas\`.

- designs:

  A character vector of design formulas for the limma analysis.

- modes:

  A character vector that must have the same length as 'designs'. For
  each design formula, you must specify either 'isolated' or
  'integrated'. Isolated means limma determines the results for each
  level using only the data from that level. Integrated means limma
  determines the results for all levels using the full dataset (from all
  levels).

- spline_test_configs:

  A configuration object for spline tests.

- report_dir:

  A non-empty string specifying the report directory.

- adj_pthresholds:

  A numeric vector of p-value thresholds for significance determination.

- rna_seq_datas:

  A list of RNA-seq data objects, such as the voom object derived from
  the limma::voom function.

- time_unit:

  A character string specifying the time unit label for plots.

- padjust_method:

  A character string specifying the method for p-value adjustment.

## Value

Returns a list of plots generated from the limma analysis results. Each
element in the list corresponds to a different combination of
hyperparameters.
