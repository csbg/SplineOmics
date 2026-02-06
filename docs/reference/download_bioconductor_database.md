# Download Gene Set Annotations from Bioconductor Organism Databases

This function extracts gene-to-ontology mappings from a specified
Bioconductor organism annotation package (e.g., \`org.Hs.eg.db\`,
\`org.Mm.eg.db\`) and saves the gene sets to a \`.tsv\` file in a
standardized format. The output includes mappings for Gene Ontology (GO)
Biological Process (BP), Molecular Function (MF), Cellular Component
(CC), and KEGG pathways. The resulting file can be used directly with
enrichment functions such as \`clusterProfiler::enricher()\` with
\`TERM2GENE\`.

## Usage

``` r
download_bioconductor_database(
  organism_db = "org.Hs.eg.db",
  output_dir = here::here(),
  filename = NULL
)
```

## Arguments

- organism_db:

  A string specifying the Bioconductor organism annotation database to
  use (e.g., \`"org.Hs.eg.db"\` for human or \`"org.Mm.eg.db"\` for
  mouse).

- output_dir:

  A string specifying the output directory where the \`.tsv\` file will
  be saved. Defaults to the current project directory as defined by
  \`here::here()\`.

- filename:

  An optional string specifying the filename for the output file. If
  \`NULL\` (default), a filename is generated automatically with a
  timestamp.

## Value

A \`data.frame\` of gene set annotations with three columns:

- DB:

  Ontology/database source, e.g. \`"GO_BP"\`, \`"GO_MF"\`, \`"GO_CC"\`,
  or \`"KEGG"\` (if available).

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
such cases the KEGG section will be empty.

In addition to returning the \`data.frame\`, the function also writes
the same table to disk as a \`.tsv\` file in the specified
\`output_dir\`.
