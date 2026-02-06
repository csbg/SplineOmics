# make_scatter_plots_html()

This function is used to make scatter plots for the raw data of all the
features. It generates an HTML report in the fashion of the other
functions of the SplineOmics package which contains all the scatter
plots.

## Usage

``` r
make_scatter_plots_html(
  data,
  meta,
  output_file = "scatter_report",
  meta_replicate_column = NULL,
  features_per_file = 500
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

- features_per_file:

  Integer specifying how many features (for example proteins) should be
  plotted per file. Per default 500.

## Value

Invisibly returns `NULL`. This function is called for its side effects:
it renders one or more HTML reports to the **current working
directory**, named `paste0(output_file, "_chunk_", i, ".html")` for
chunk index `i`, and prints progress messages during generation.

## Examples

``` r
# Create toy data: 3 features × 6 samples
set.seed(123)
toy_data <- matrix(
  rnorm(18, mean = 5, sd = 2),
  nrow = 3,
  ncol = 6,
  dimnames = list(
    c("geneA", "geneB", "geneC"),
    paste0("sample", 1:6)
  )
)

# Meta data: must include a numeric Time column
toy_meta <- data.frame(
  Time = rep(c(0, 1, 2), each = 2),
  Replicate = rep(c("R1", "R2"), times = 3),
  row.names = colnames(toy_data)
)

# Write HTML reports into a temporary directory
old_wd <- setwd(tempdir())
make_scatter_plots_html(
  data = toy_data,
  meta = toy_meta,
  output_file = "scatter_demo",
  meta_replicate_column = "Replicate",
  features_per_file = 2
)
#> Total features: 3
#> Generating 2 HTML reports in chunks of 2
#> Warning: The `size` argument of `element_line()` is deprecated as of ggplot2 3.4.0.
#> ℹ Please use the `linewidth` argument instead.
#> ℹ The deprecated feature was likely used in the SplineOmics package.
#>   Please report the issue at <https://github.com/csbg/SplineOmics/issues>.
#> 
#> 
#> processing file: file235a167767d.Rmd
#> 1/1
#> output file: file235a167767d.knit.md
#> /usr/lib/rstudio/resources/app/bin/quarto/bin/tools/x86_64/pandoc +RTS -K512m -RTS file235a167767d.knit.md --to html4 --from markdown+autolink_bare_uris+tex_math_single_backslash --output scatter_demo_chunk_1.html --lua-filter /home/thomas/.cache/R/renv/library/SplineOmics-2e6a5dbb/linux-ubuntu-jammy/R-4.5/x86_64-pc-linux-gnu/rmarkdown/rmarkdown/lua/pagebreak.lua --lua-filter /home/thomas/.cache/R/renv/library/SplineOmics-2e6a5dbb/linux-ubuntu-jammy/R-4.5/x86_64-pc-linux-gnu/rmarkdown/rmarkdown/lua/latex-div.lua --lua-filter /home/thomas/.cache/R/renv/library/SplineOmics-2e6a5dbb/linux-ubuntu-jammy/R-4.5/x86_64-pc-linux-gnu/rmarkdown/rmarkdown/lua/table-classes.lua --embed-resources --standalone --variable bs3=TRUE --section-divs --template /home/thomas/.cache/R/renv/library/SplineOmics-2e6a5dbb/linux-ubuntu-jammy/R-4.5/x86_64-pc-linux-gnu/rmarkdown/rmd/h/default.html --no-highlight --variable highlightjs=1 --variable theme=bootstrap --mathjax --variable 'mathjax-url=https://mathjax.rstudio.com/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML' --include-in-header /tmp/Rtmpheq4uA/rmarkdown-str235a3ba1f2e0.html 
#> 
#> Output created: scatter_demo_chunk_1.html
#> [DONE] Rendered: scatter_demo_chunk_1.html
#> 
#> 
#> processing file: file235a11161a5.Rmd
#> 1/1
#> output file: file235a11161a5.knit.md
#> /usr/lib/rstudio/resources/app/bin/quarto/bin/tools/x86_64/pandoc +RTS -K512m -RTS file235a11161a5.knit.md --to html4 --from markdown+autolink_bare_uris+tex_math_single_backslash --output scatter_demo_chunk_2.html --lua-filter /home/thomas/.cache/R/renv/library/SplineOmics-2e6a5dbb/linux-ubuntu-jammy/R-4.5/x86_64-pc-linux-gnu/rmarkdown/rmarkdown/lua/pagebreak.lua --lua-filter /home/thomas/.cache/R/renv/library/SplineOmics-2e6a5dbb/linux-ubuntu-jammy/R-4.5/x86_64-pc-linux-gnu/rmarkdown/rmarkdown/lua/latex-div.lua --lua-filter /home/thomas/.cache/R/renv/library/SplineOmics-2e6a5dbb/linux-ubuntu-jammy/R-4.5/x86_64-pc-linux-gnu/rmarkdown/rmarkdown/lua/table-classes.lua --embed-resources --standalone --variable bs3=TRUE --section-divs --template /home/thomas/.cache/R/renv/library/SplineOmics-2e6a5dbb/linux-ubuntu-jammy/R-4.5/x86_64-pc-linux-gnu/rmarkdown/rmd/h/default.html --no-highlight --variable highlightjs=1 --variable theme=bootstrap --mathjax --variable 'mathjax-url=https://mathjax.rstudio.com/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML' --include-in-header /tmp/Rtmpheq4uA/rmarkdown-str235a3c5e7fc2.html 
#> 
#> Output created: scatter_demo_chunk_2.html
#> [DONE] Rendered: scatter_demo_chunk_2.html
#> All reports generated successfully.
setwd(old_wd)

# Inspect generated HTML files in tempdir():
list.files(tempdir(), pattern = "scatter_demo_chunk_.*html$")
#> [1] "scatter_demo_chunk_1.html" "scatter_demo_chunk_2.html"
```
