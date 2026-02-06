# Download Enrichr databases

Download gene sets from one or more Enrichr libraries and return them as
a table with one row per gene-set membership.

If `output_dir` is provided, the same table is also written to disk as a
tab-separated file (`.tsv`). If `filename` is `NULL`, a timestamped name
is generated to ensure uniqueness.

## Usage

``` r
download_enrichr_databases(gene_set_lib, output_dir = NULL, filename = NULL)
```

## Arguments

- gene_set_lib:

  [`character()`](https://rdrr.io/r/base/character.html). A character
  vector of Enrichr library names to download, e.g.
  `c("WikiPathways_2019_Human", "NCI-Nature_2016")`.

- output_dir:

  `character(1)` or `NULL`. Output directory for writing a `.tsv` file.
  If `NULL`, no file is written.

- filename:

  `character(1)` or `NULL`. Output file name (not a path). If `NULL` and
  `output_dir` is not `NULL`, a default name of the form
  `enrichr_databases_YYYYmmdd-HHMMSS.tsv` is used. Due to commas in some
  terms, `.tsv` is recommended.

## Value

A `data.frame` with three columns:

- DB:

  Enrichr library name.

- Geneset:

  Gene set or pathway term within that library.

- Gene:

  Gene symbol contained in the gene set.

If a requested library cannot be downloaded, it may be omitted from the
result. The function errors if no gene sets can be retrieved.

## Examples

``` r
if (interactive()) {
  libs <- c("WikiPathways_2019_Human")
  out <- download_enrichr_databases(
    gene_set_lib = libs,
    output_dir = tempdir(),
    filename = "enrichr_demo.tsv"
  )
  head(out)
}
```
