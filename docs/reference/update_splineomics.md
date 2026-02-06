# Update the variables in a SplineOmics object

Updates a SplineOmics object by modifying existing fields or adding new
ones.

## Usage

``` r
update_splineomics(splineomics, ...)
```

## Arguments

- splineomics:

  `SplineOmics`: A SplineOmics object to be updated.

- ...:

  Named arguments with new values for fields to be updated or added.

## Value

The updated SplineOmics object.

## Examples

``` r
set.seed(1)
toy_data <- matrix(rnorm(12),
    nrow = 3,
    dimnames = list(paste0("gene", 1:3), paste0("S", 1:4))
)
toy_meta <- data.frame(
    SampleID = colnames(toy_data),
    Condition = c("Ctrl", "Ctrl", "Trt", "Trt"),
    stringsAsFactors = FALSE,
    row.names = colnames(toy_data)
)

so <- create_splineomics(
    data = toy_data,
    meta = toy_meta,
    condition = toy_meta$Condition
)

# Update the mode and add a new design matrix
new_design <- model.matrix(~Condition, data = toy_meta)
so_updated <- update_splineomics(so,
    mode = "integrated",
    design = new_design
)

str(so_updated, max.level = 1)
#> List of 16
#>  $ data                : num [1:3, 1:4] -0.626 0.184 -0.836 1.595 0.33 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ rna_seq_data        : NULL
#>  $ meta                :'data.frame':    4 obs. of  2 variables:
#>  $ condition           : chr [1:4] "Ctrl" "Ctrl" "Trt" "Trt"
#>  $ annotation          : NULL
#>  $ report_info         : NULL
#>  $ meta_batch_column   : NULL
#>  $ meta_batch2_column  : NULL
#>  $ feature_name_columns: NULL
#>  $ design              : num [1:4, 1:2] 1 1 1 1 0 0 1 1
#>   ..- attr(*, "dimnames")=List of 2
#>   ..- attr(*, "assign")= int [1:2] 0 1
#>   ..- attr(*, "contrasts")=List of 1
#>  $ use_array_weights   : logi FALSE
#>  $ dream_params        : NULL
#>  $ mode                : chr "integrated"
#>  $ spline_params       : NULL
#>  $ padjust_method      : chr "BH"
#>  $ bp_cfg              : NULL
#>  - attr(*, "class")= chr "SplineOmics"
```
