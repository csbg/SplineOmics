# Extract gene set annotations from Bioconductor organism databases

This function extracts gene-to-ontology mappings from a specified
Bioconductor organism annotation package (e.g., `org.Hs.eg.db`,
`org.Mm.eg.db`) and saves the gene sets to a `.tsv` file in a
standardized format. The output includes mappings for Gene Ontology (GO)
Biological Process (BP), Molecular Function (MF), Cellular Component
(CC), and KEGG pathways. The resulting file can be used directly with
enrichment functions such as
[`clusterProfiler::enricher()`](https://rdrr.io/pkg/clusterProfiler/man/enricher.html)
with `TERM2GENE`.

## Usage

``` r
extract_gene_sets(
  organism_db = "org.Hs.eg.db",
  output_dir = tempdir(),
  filename = NULL
)
```

## Arguments

- organism_db:

  `character(1)`: A string specifying the Bioconductor organism
  annotation database to use (e.g., `"org.Hs.eg.db"` for human or
  `"org.Mm.eg.db"` for mouse).

- output_dir:

  `character(1)`: A string specifying the output directory where the
  `.tsv` file will be saved.

- filename:

  `character(1)` \| `NULL`: An optional string specifying the filename
  for the output file. If `NULL` (default), a filename is generated
  automatically with a timestamp.

## Value

A `data.frame` of gene set annotations with three columns:

- DB:

  Ontology/database source, e.g. `"GO_BP"`, `"GO_MF"`, `"GO_CC"`, or
  `"KEGG"` (if available).

- Geneset:

  Ontology term ID or pathway ID (e.g. GO ID, KEGG ID).

- Gene:

  Gene symbol (`SYMBOL`).

## Details

The TSV has three columns:

- DB:

  Ontology/database source, e.g., `"GO_BP"`, `"GO_MF"`, `"GO_CC"`, or
  `"KEGG"` (if available).

- Geneset:

  Ontology term ID or pathway ID (e.g., GO ID, KEGG ID).

- Gene:

  Gene symbol (`SYMBOL`).

Note: Some `org.*.eg.db` packages no longer include KEGG mappings; in
such cases the KEGG mappings may be absent.

In addition to returning the `data.frame`, the function also writes the
same table to disk as a `.tsv` file in the specified `output_dir`.

## Examples

``` r
# Minimal real example (runs only if org package is installed)
tmp <- tempdir()
if (requireNamespace("org.Mm.eg.db", quietly = TRUE) &&
    requireNamespace("AnnotationDbi", quietly = TRUE)) {
    gs <- extract_gene_sets(
        organism_db = "org.Mm.eg.db",
        output_dir  = tmp,
        filename    = "mm_genesets.tsv"
    )
    head(gs)
    # The file path:
    file.path(tmp, "mm_genesets.tsv")
}
#> 
#> 'select()' returned 1:many mapping between keys and columns
#> 'select()' returned 1:many mapping between keys and columns
#> 
#> Gene set extraction complete complete! The file has been saved as: /tmp/RtmpyRGyv7/mm_genesets.tsv
#> [1] "/tmp/RtmpyRGyv7/mm_genesets.tsv"

# If the organism package is not installed, you can still see the TSV format:
tiny <- data.frame(
    DB = c("GO_BP", "GO_MF"),
    Geneset = c("GO:0008150", "GO:0003674"),
    Gene = c("Trp53", "Egfr"),
    stringsAsFactors = FALSE
)
utils::write.table(
    tiny,
    file = file.path(tmp, "example_genesets.tsv"),
    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
)
```
