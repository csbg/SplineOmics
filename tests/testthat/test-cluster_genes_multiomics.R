testthat::test_that(
    "cluster_genes_multiomics works.",
    {
        set.seed(123)
        
        genes <- paste0("gene", 1:6)
        
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
        
        blocks <- list(
            time_Ctrl = list(
                rna = rna_block1,
                site = site_block1
            ),
            time_Treat = list(
                rna = rna_block2,
                site = site_block2
            )
        )
        
        block_clusters <- list(
            time_Ctrl = 2L,
            time_Treat = 2L
        )
        
        modality_meta <- data.frame(
            block = c(
                "time_Ctrl", "time_Ctrl",
                "time_Treat", "time_Treat"
            ),
            layer = c("rna", "site", "rna", "site"),
            layer_k = c(NA_real_, 3, NA_real_, 3),
            layer_w = c(1, 1, 1, 1),
            stringsAsFactors = FALSE
        )
        
        clustering_results <- cluster_genes_multiomics(
            blocks = blocks,
            block_clusters = block_clusters,
            modality_meta = modality_meta,
            gene_mode = "union",
            verbose = FALSE
        )
        
        testthat::expect_snapshot_value(
            clustering_results$cluster_table,
            style = "deparse"
        )
        testthat::expect_snapshot_value(
            clustering_results$centroid_info,
            style = "deparse"
        )
    }
)
