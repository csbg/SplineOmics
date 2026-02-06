# Perform over-representation analysis with the results from cluster_hits()

This function performs over-representation analysis (ORA) on clustered
feature sets using the R package **clusterProfiler**. For each grouping
column in the provided cluster table, enrichment is evaluated separately
for all clusters and specified databases.

The function computes and returns the raw enrichment results only. No
plots are generated and no files are written. Visualization and HTML
report generation are intentionally handled by a separate reporting
function (e.g.,
[`create_ora_report()`](https://csbg.github.io/SplineOmics/reference/create_ora_report.md)),
enabling a strict separation between statistical analysis and reporting.

## Usage

``` r
run_ora(
  cluster_table,
  databases,
  clusterProfiler_params = NA,
  mapping_cfg = list(method = "none", from_species = NULL, to_species = NULL),
  enrichGO_cfg = NULL,
  universe = NULL,
  verbose = FALSE
)
```

## Arguments

- cluster_table:

  `tibble` A tibble with one row per feature and at least one column
  named `gene`. The `gene` column must contain gene identifiers as
  `character(1)` (empty strings and `NA` are ignored).

  All remaining columns (except `feature_nr` and `feature_name`) are
  treated as *grouping columns*. Each grouping column defines a set of
  categories (clusters) via its distinct non-`NA` values. For a given
  grouping column, all rows that share the same value form one cluster,
  and overrepresentation analysis is performed separately for each
  cluster. The data type of grouping columns is not restricted (e.g.,
  integer, character, factor); values are compared by equality after
  being coerced to character.

  A value of `NA` in a grouping column indicates that the feature does
  not belong to any cluster for that grouping column and is excluded
  from the corresponding analysis.

  The columns `feature_nr` and `feature_name`, if present, are ignored.
  These columns are typically added by upstream functions in the
  SplineOmics pipeline.

- databases:

  `data.frame`: A `data.frame` that defines the gene set collections to
  be tested in the overrepresentation analysis. Must contain exactly
  three columns:

  - `DB`: `character(1)` The database identifier (e.g., KEGG, GO_BP,
    Reactome).

  - `Geneset`: `character(1)` The name of the gene set or pathway within
    the database.

  - `Gene`: `character(1)` A gene identifier belonging to the gene set
    (e.g., gene symbol, Ensembl ID).

  Each row corresponds to one `(database, geneset, gene)` association.
  The same gene may appear in multiple gene sets.

- clusterProfiler_params:

  `list` \| `NULL`: A named list of arguments passed directly to the
  corresponding functions in the **clusterProfiler** package. Typical
  entries include `pvalueCutoff`, `pAdjustMethod`, `minGSSize`,
  `maxGSSize`, and `qvalueCutoff`. The names must match the argument
  names in clusterProfiler; see the clusterProfiler documentation for
  details. If `NULL` (default), the standard clusterProfiler defaults
  are used.

- mapping_cfg:

  `list` \| `NULL`: A named list that controls the optional behavior of
  automatically mapping gene symbols across species. This is useful when
  your input gene symbols (e.g., from CHO cells) do not match the
  species used by the enrichment databases (e.g., human or mouse). By
  default, no mapping is performed and gene symbols are used as-is. If
  mapping is desired, this list must contain the following three
  elements:

  method

  :   `character(1)`: Mapping method to use. One of `none` (default; no
      mapping), `gprofiler` (online, via the g:Profiler API), or
      `orthogene` (offline, if installed).

  from_species

  :   `character(1)`: Source species code, e.g. `cgriseus` for CHO. Must
      match the expected format for the selected tool.

  to_species

  :   `character(1)`: Target species code, e.g. `hsapiens` for human.
      This must be the species used in your ORA database and must also
      match the expected format for the selected tool.

- enrichGO_cfg:

  `list` \| `NULL`: A named list specifying the configuration for
  running GO enrichment with Bioconductor's
  [`enrichGO`](https://rdrr.io/pkg/clusterProfiler/man/enrichGO.html).
  This is only needed when you want to perform GO Biological Process
  (BP), Molecular Function (MF), or Cellular Component (CC) enrichment
  using Bioconductor's organism databases (e.g., `org.Mm.eg.db` for
  mouse).

  The list must be named according to the GO ontology, e.g., `"GO_BP"`,
  `"GO_MF"`, `"GO_CC"`. Each entry must provide:

  - `OrgDb`: `character(1)` The organism database, e.g., `org.Mm.eg.db`.

  - `keyType`: `character(1)` The gene identifier type, e.g.,
    `"SYMBOL"`.

  - `ontology`: `character(1)` One of `"BP"`, `"MF"`, or `"CC"`.

  If `enrichGO_cfg` is `NULL` (default), no Bioconductor-based GO
  enrichment is performed. All enrichment runs through
  [`enricher`](https://rdrr.io/pkg/clusterProfiler/man/enricher.html)
  with the provided TERM2GENE mappings.

- universe:

  [`character()`](https://rdrr.io/r/base/character.html) \| `NULL`:
  Enrichment background data. This is a parameter of clusterProfiler;
  for details, please check the documentation of the clusterProfiler R
  package.

- verbose:

  `logical(1)`: Boolean flag controlling the display of messages.

## Value

A named list with two elements:

- `all_results`:

  A nested, named list containing the raw over-representation analysis
  results for each grouping column in `cluster_table`. Each top-level
  element corresponds to one grouping column and contains the field
  `ora_results`, a nested list of enrichment results organized by
  cluster and database.

- `report_payload`:

  A structured list containing all information required to generate an
  HTML ORA report without recomputing any enrichment results. This
  payload is intended to be passed to
  [`create_ora_report()`](https://csbg.github.io/SplineOmics/reference/create_ora_report.md).

## Examples

``` r
{
    set.seed(1)

    # toy cluster table (two "conditions")
    toy_genes <- paste0("G", 1:8)
    cluster_table <- tibble::tibble(
        gene          = toy_genes,
        cluster_condA = c(1, 1, 2, 2, NA, NA, 1, 2),
        cluster_condB = c(NA, 1, NA, 2, 1, 2, 1, NA)
    )

    # toy TERM2GENE database
    databases <- data.frame(
        DB = rep("ToyDB", 6),
        Geneset = c(rep("SetA", 3), rep("SetB", 3)),
        Gene = c("G1", "G2", "G7", "G3", "G4", "G6"),
        stringsAsFactors = FALSE
    )

    # permissive params for tiny example 
    clusterProfiler_params <- list(
        pvalueCutoff = 1,
        minGSSize    = 1,
        maxGSSize    = 500
    )

    # run ORA
    res <- run_ora(
        cluster_table            = cluster_table,
        databases                = databases,
        clusterProfiler_params   = clusterProfiler_params,
        verbose                  = TRUE
    )

    # see sections and files written
    names(res)
}
#> 
#> 
#>  Running clusterProfiler for column: cluster_condA
#> 
#> Cluster: cluster_1
#> Database: ToyDB
#> Foreground genes:3
#> Foreground genes overlapping with database: 3 (100%)
#> 
#> Cluster: cluster_2
#> Database: ToyDB
#> Foreground genes:3
#> Foreground genes overlapping with database: 2 (66.7%)
#> 
#> 
#>  Running clusterProfiler for column: cluster_condB
#> 
#> Cluster: cluster_1
#> Database: ToyDB
#> Foreground genes:3
#> Foreground genes overlapping with database: 2 (66.7%)
#> 
#> Cluster: cluster_2
#> Database: ToyDB
#> Foreground genes:2
#> Foreground genes overlapping with database: 2 (100%)
#> [1] "all_results"    "report_payload"
```
