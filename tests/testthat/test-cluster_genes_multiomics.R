test_that("cluster_genes_multiomics snapshot with mixed layers", {
    set.seed(123)
    
    # genes
    genes <- paste0("gene", 1:6)
    
    # gene-level layer: rna trajectories for two blocks
    rna_block1 <- matrix(
        rnorm(length(genes) * 5),
        nrow = length(genes),
        ncol = 5,
        dimnames = list(genes, NULL)
    )
    rna_block2 <- matrix(
        rnorm(length(genes) * 5),
        nrow = length(genes),
        ncol = 5,
        dimnames = list(genes, NULL)
    )
    
    # many-to-one layer: 3 sites per gene, rownames "<gene>_<site>"
    sites <- paste0("s", 1:3)
    feat_ids <- paste0(
        rep(genes, each = length(sites)),
        "_",
        rep(sites, times = length(genes))
    )
    
    site_block1 <- matrix(
        rnorm(length(feat_ids) * 5),
        nrow = length(feat_ids),
        ncol = 5,
        dimnames = list(feat_ids, NULL)
    )
    site_block2 <- matrix(
        rnorm(length(feat_ids) * 5),
        nrow = length(feat_ids),
        ncol = 5,
        dimnames = list(feat_ids, NULL)
    )
    
    # construct blocks: two blocks, two layers each (rna + many-to-one site)
    blocks <- list(
        time_Ctrl = list(
            rna  = rna_block1,
            site = site_block1
        ),
        time_Treat = list(
            rna  = rna_block2,
            site = site_block2
        )
    )
    
    # block-level metadata
    block_meta <- data.frame(
        block           = c("time_Ctrl", "time_Treat"),
        block_k         = c(2L, 2L),
        result_category = c(1, 1),
        cond1           = c("Ctrl", "Treat"),
        cond2           = c(NA_character_, NA_character_),
        stringsAsFactors = FALSE
    )
    
    # layer-level metadata:
    # - rna: one-to-one (layer_k = NA)
    # - site: many-to-one (layer_k > 0)
    modality_meta <- data.frame(
        block   = c("time_Ctrl", "time_Ctrl",
                    "time_Treat", "time_Treat"),
        layer   = c("rna", "site", "rna", "site"),
        layer_k = c(NA_real_, 3, NA_real_, 3),
        layer_w = c(1, 1, 1, 1),
        stringsAsFactors = FALSE
    )
    
    clustering_results <- cluster_genes_multiomics(
        blocks     = blocks,
        block_meta = block_meta,
        modality_meta = modality_meta,
        gene_mode  = "union",
        verbose = FALSE
    )

    expect_snapshot_value(
        clustering_results$cluster_table,
        style = "deparse"
        )
    expect_snapshot_value(
        clustering_results$centroid_info,
        style = "deparse"
        )
})

