# Gene-centric multi-omics clustering across data modalities

Performs gene-centric clustering of multi-omics feature representations
that are supplied as modality-specific matrices. The function harmonizes
genes across modalities, constructs a gene-centric joint feature matrix,
computes a UMAP neighborhood graph using uwot, and performs spectral
clustering on the resulting UMAP fuzzy graph (`fgraph`).

The input matrices are assumed to be precomputed gene-level trajectories
(e.g., spline values or coefficients) or many-to-one feature-level
matrices (e.g., phospho sites) that are summarized into gene-level
pattern signature vectors inside this function.

This function performs statistical computation and returns the UMAP fit
for downstream visualization. Clustering is performed on the UMAP graph
rather than on the two-dimensional embedding coordinates.

## Usage

``` r
cluster_genes_multiomics(
  data,
  meta,
  k,
  gene_mode = c("intersection", "union"),
  n_neighbors = 15L,
  min_graph_purity = 0,
  verbose = FALSE
)
```

## Arguments

- data:

  Named list defining the multi-condition, multi-modality input data.
  The outer list corresponds to experimental conditions (e.g. control,
  and treatment). Each element of the outer list must itself be a named
  list of numeric matrices, one per data modality.

  For each condition, the inner list entries represent modalities and
  must share the same modality names across conditions. Each modality
  matrix must have row names and can be of two types:

  One-to-one (gene-level) modality

  :   Rows represent genes directly. Row names must be gene identifiers
      in the format: `<gene_id>`. Columns represent the
      modality-specific representation used for clustering, such as
      spline values at fixed time points, spline coefficients, or other
      numeric features.

  Many-to-one (feature-level) modality

  :   Rows represent features that map to genes (e.g., phospho sites,
      probes). Row names must follow the pattern
      `<gene_id>_<feature_id>` where the gene identifier precedes the
      first underscore. Such modalities are aggregated into gene-level
      pattern signature vectors using an internal feature clustering
      step prior to integration.

  The \<gene_id\> parts of the rownames of all modalities across all
  conditions should match, otherwise, the gene-centric clustering is not
  possible! Modality matrices may differ in their number of columns both
  within and across conditions. After aggregation (for many-to-one
  modalities) and normalization, all condition- and modality-specific
  gene-level matrices are aligned according to `gene_mode` and
  concatenated to form the gene-centric feature matrix used to construct
  the UMAP graph and clustering.

  When multiple conditions are supplied, many-to-one modalities are
  processed by pooling all feature-level observations across conditions
  prior to constructing gene-level pattern signatures. This means that
  the global archetype basis used for feature aggregation is learned
  jointly across conditions.

  As a consequence, gene representations from different conditions are
  expressed in a shared feature space, enabling direct comparison of
  gene trajectories across conditions. Running the function separately
  for each condition would instead learn independent archetype bases and
  may therefore produce different gene-level representations for
  many-to-one modalities.

- meta:

  Data frame with one row per modality, providing modality-level
  parameters required for gene-centric representation building. The data
  frame must contain at least the following columns:

  `modality`

  :   Character scalar giving the modality identifier. Each value must
      match a modality name present in `data[[condition]]`. The order of
      rows defines the modality ordering used throughout downstream
      processing.

  `many_to_one_k`

  :   Integer or `NA`. If `NA`, the modality is treated as one-to-one
      (gene-level) and passed through unchanged. If an integer, the
      modality is treated as many-to-one (feature-level) and collapsed
      to gene-level signatures using `many_to_one_k` global archetypes.

  `modality_w`

  :   Numeric, non-negative scalar giving the relative weight of the
      modality in the joint feature space. Weights are normalized
      internally to sum to one across modalities before being applied.

  Additional columns may be present but are ignored.

- k:

  Integer. Number of gene clusters for spectral clustering.

- gene_mode:

  Character string specifying how to harmonize genes across modalities
  prior to constructing the joint feature matrix:

  `"intersection"`

  :   Retain only genes present in all modalities.

  `"union"`

  :   Retain genes present in any modality. Gene vectors are constructed
      using available modalities; missing modality blocks are handled
      internally.

- n_neighbors:

  Integer. Size of the local neighborhood used by UMAP to construct the
  k-nearest-neighbor graph. Larger values emphasize broader structure,
  smaller values emphasize local structure.

- min_graph_purity:

  Numeric scalar in \\\[0, 1\]\\ controlling post-clustering filtering
  based on graph-based assignment purity. For each gene, purity is
  defined as the fraction of its total graph connectivity (`fgraph`)
  that lies within its assigned cluster. Genes with
  `purity < min_graph_purity` are reassigned to cluster `0` ("other").
  Default `0` disables filtering.

- verbose:

  Logical scalar indicating whether progress messages should be emitted
  via [`rlang::inform()`](https://rlang.r-lib.org/reference/abort.html)
  by internal helpers.

## Value

A named list with the following elements:

- `cluster_table`:

  A tibble with one row per gene and columns `gene` and `cluster`.

- `centroid_info`:

  A tibble with one row per `(condition, modality, cluster)` that
  summarizes cluster centroids, block coverage, and within-cluster
  coherence statistics. The centroid is stored as a list-column and its
  interpretation depends on modality type: for one-to-one modalities the
  centroid is a mean normalized trajectory/feature vector (typically
  row-wise z-scored), whereas for many-to-one modalities the centroid is
  a raw archetype mixture vector (signature fractions over '
  `many_to_one_k` global archetypes; not Hellinger-transformed) and
  therefore sums to 1 when defined. Coherence is reported via
  `qc_method` and `mean_qc`/`sd_qc`, computed in the normalized
  clustering space (Pearson R\\^2\\ for one-to-one, Bhattacharyya
  coefficient in Hellinger space for many-to-one). The columns
  `centroid_type` and `centroid_feature_names` provide machine-readable
  metadata for plotting the centroid vector.

- `many_to_one_clustering_qc`:

  A tibble (or `NULL`) with QC diagnostics for many-to-one feature
  clustering steps, if any are present.

- `many_to_one_archetypes`:

  A named list (or `NULL`) with the global archetypes learned for each
  many-to-one modality. Archetypes are fit once per modality on pooled
  feature-level trajectories across all supplied conditions. Each entry
  is keyed by modality name and contains: `archetypes` (a `k x p`
  numeric matrix of archetype patterns), `k` (the number of archetypes),
  and `feature_names` (column names for the `p` features / time points).

- `umap_fit`:

  The object returned by
  [`uwot::umap()`](https://jlmelville.github.io/uwot/reference/umap.html),
  including `$embedding` (genes x `n_components`) and `$fgraph` (the
  UMAP fuzzy graph) when requested via `ret_extra`.

## Examples

``` r
set.seed(1)
genes <- paste0("gene", 1:8)

rna <- matrix(
  rnorm(length(genes) * 5),
  nrow = length(genes),
  dimnames = list(genes, NULL)
)

prot <- matrix(
  rnorm(length(genes) * 3),
  nrow = length(genes),
  dimnames = list(genes, NULL)
)

data <- list(
  condition1 = list(rna = rna, protein = prot)
)

meta <- data.frame(
  modality = c("rna", "protein"),
  many_to_one_k = c(NA_real_, NA_real_),
  modality_w = c(1, 1),
  stringsAsFactors = FALSE
)

res <- cluster_genes_multiomics(
  data = data,
  meta = meta,
  k = 3L,
  gene_mode = "intersection",
  n_neighbors = 5L
)

res$cluster_table
#> # A tibble: 8 × 2
#>   gene  cluster
#>   <chr>   <int>
#> 1 gene1       1
#> 2 gene2       3
#> 3 gene3       2
#> 4 gene4       2
#> 5 gene5       3
#> 6 gene6       3
#> 7 gene7       1
#> 8 gene8       1
```
