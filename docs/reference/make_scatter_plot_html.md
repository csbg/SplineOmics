# make_scatter_plot_html()

This function is used to make scatter plots for the raw data of all the
features. It generates an HTML report in the fashion of the other
functions of the SplineOmics package which contains all the scatter
plots.

## Usage

``` r
make_scatter_plot_html(
  data,
  meta,
  output_file = "scatter_plots_base64.html",
  meta_replicate_column = NULL
)
```

## Arguments

- data:

  A matrix with features as rows and samples as columns. Row names
  should be feature names.

- meta:

  A data frame with the meta information. Must contain a numeric column
  "Time".

- output_file:

  The name of the HTML output file.

- meta_replicate_column:

  Column name of the column in meta that contains the info about the
  replicates, such as reactor.

## See also

[`render`](https://pkgs.rstudio.com/rmarkdown/reference/render.html)

## Examples

``` r
if (FALSE) { # \dontrun{
# Example Data
data <- matrix(rnorm(50), nrow = 5)
meta <- data.frame(Time = seq(1, 10, length.out = 10))

# Generate HTML report (only if you want to test it)
make_scatter_plot_html(data, meta, "scatter_report.html")
} # }
```
