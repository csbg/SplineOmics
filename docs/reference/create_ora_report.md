# Generate an HTML report for over-representation analysis results

This function generates an HTML report for over-representation analysis
(ORA) results produced by
[`run_ora()`](https://csbg.github.io/SplineOmics/reference/run_ora.md).
It consumes a structured `report_payload` object and creates dotplots
summarizing enrichment results for each grouping column in the cluster
table.

The function is responsible for all visualization and reporting side
effects. Statistical analysis is not recomputed. Calling this function
implies that an HTML report should be generated; therefore, a valid
output directory must be provided.

## Usage

``` r
create_ora_report(
  report_payload,
  report_info = NULL,
  cluster_hits_report_name = NULL,
  verbose = FALSE,
  report_dir = tempdir()
)
```

## Arguments

- report_payload:

  `list`: A structured ORA report payload as returned by
  [`run_ora()`](https://csbg.github.io/SplineOmics/reference/run_ora.md).
  It must contain the raw enrichment results and all metadata required
  for report generation.

- report_info:

  `list` \| `NULL`: Optional list with experiment metadata used to
  annotate the HTML report (e.g., omics data type, project name,
  analyst, and data description). If `NULL`, a minimal report is
  generated.

- cluster_hits_report_name:

  `character(1)` \| `NULL`: Optional name of the
  [`cluster_hits()`](https://csbg.github.io/SplineOmics/reference/cluster_hits.md)
  report that produced the clustering results used for ORA. When
  provided, it is displayed in the report to document provenance.

- verbose:

  `logical(1)`: Logical flag controlling the display of progress and
  status messages.

- report_dir:

  `character(1)`: Directory where the HTML report and all associated
  output files are written. The directory is created if it does not
  already exist.

## Value

A named list of `ggplot` objects containing the ORA dotplots generated
for the report. Each element is named after the corresponding grouping
column in the cluster table, allowing easy identification and reuse of
the plots outside the HTML report.

## See also

[`run_ora()`](https://csbg.github.io/SplineOmics/reference/run_ora.md),
[`cluster_hits()`](https://csbg.github.io/SplineOmics/reference/cluster_hits.md),
[`clusterProfiler::enricher()`](https://rdrr.io/pkg/clusterProfiler/man/enricher.html)

## Examples

``` r
# Toy example: build a minimal payload and generate a report in tempdir().
# This runs only if the internal report helpers are available.

tmp <- tempdir()

# Minimal cluster table with two cluster labels
cluster_table <- data.frame(
  gene = c("GeneA", "GeneB", "GeneC", "GeneD"),
  cluster = c("1", "1", "2", "2"),
  stringsAsFactors = FALSE
)

# Minimal "databases" descriptor as expected by the report
databases <- data.frame(
  DB = c("GO_BP"),
  stringsAsFactors = FALSE
)

# A tiny ORA result table resembling clusterProfiler output
toy_enrich <- data.frame(
  ID = c("GO:0008150", "GO:0009987"),
  Description = c("biological_process", "cellular process"),
  p.adjust = c(0.01, 0.02),
  Count = c(2L, 1L),
  stringsAsFactors = FALSE
)

# Minimal nested structure:
# all_results[[section]][[cluster]][[db]] -> data.frame
all_results <- list(
  cluster = list(
    `1` = list(GO_BP = toy_enrich),
    `2` = list(GO_BP = toy_enrich)
  )
)

report_payload <- list(
  all_results = all_results,
  databases = databases,
  clusterProfiler_params = list(pAdjustMethod = "BH", pvalueCutoff = 0.05),
  cluster_table = cluster_table,
  universe = unique(cluster_table$gene)
)

# Only run if the report backend is present (avoids failures on minimal
# check environments where reporting deps might not be installed).
if (exists("generate_report_html", mode = "function")) {
  plots <- suppressMessages(
    create_ora_report(
      report_payload = report_payload,
      report_info = list(project = "Toy ORA"),
      cluster_hits_report_name = "toy_cluster_hits",
      verbose = FALSE,
      report_dir = file.path(tmp, "splineomics_ora_toy_report")
    )
  )

  # Inspect returned ggplot objects (one per section)
  names(plots)
  plots[[1]]
}
```
