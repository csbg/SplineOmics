# Print function for SplineOmics objects

This function provides a summary print of the SplineOmics object,
showing relevant information such as the number of features, samples,
metadata, RNA-seq data, annotation, and spline parameters.

## Usage

``` r
# S3 method for class 'SplineOmics'
print(x, ...)
```

## Arguments

- x:

  `SplineOmics`: A SplineOmics object created by the
  `create_splineomics` function.

- ...:

  Additional arguments passed to or from other methods.

## Value

The function does not return a value. It prints a summary of the
SplineOmics object.

## Details

This function is automatically called when a SplineOmics object is
printed. It provides a concise overview of the object's contents and
attributes, including the dimensions of the data, available metadata,
and other relevant information such as annotations and spline
parameters.

## Examples

``` r
# Example: create and print a SplineOmics object
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
    condition = toy_meta$Condition,
    spline_params = list(spline_type = "n", dof = 3),
    padjust_method = "BH"
)

# The print method is automatically called:
so
#> data:SplineOmics Object
#> -------------------
#> Number of features (rows): 3 
#> Number of samples (columns): 4 
#> Meta data columns: 2 
#> First few meta columns:
#>    SampleID Condition
#> S1       S1      Ctrl
#> S2       S2      Ctrl
#> S3       S3       Trt
#> Condition: Ctrl Ctrl Trt Trt 
#> No RNA-seq data provided.
#> No annotation provided.
#> Spline parameters are set:
#> $spline_type
#> [1] "n"
#> 
#> $dof
#> [1] 3
#> 
#> P-value adjustment method: BH 

# Or explicitly:
print(so)
#> data:SplineOmics Object
#> -------------------
#> Number of features (rows): 3 
#> Number of samples (columns): 4 
#> Meta data columns: 2 
#> First few meta columns:
#>    SampleID Condition
#> S1       S1      Ctrl
#> S2       S2      Ctrl
#> S3       S3       Trt
#> Condition: Ctrl Ctrl Trt Trt 
#> No RNA-seq data provided.
#> No annotation provided.
#> Spline parameters are set:
#> $spline_type
#> [1] "n"
#> 
#> $dof
#> [1] 3
#> 
#> P-value adjustment method: BH 
```
