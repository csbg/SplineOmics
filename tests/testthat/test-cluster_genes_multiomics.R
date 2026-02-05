testthat::test_that(
    "cluster_genes_multiomics works.",
    {
        set.seed(123)
        
        genes_all <- paste0("gene", 1:12)
        genes_ctrl <- genes_all
        genes_treat <- genes_all[-c(3, 9)]
        
        rna_ctrl <- matrix(
            rnorm(length(genes_ctrl) * 6),
            nrow     = length(genes_ctrl),
            ncol     = 6,
            dimnames = list(genes_ctrl, paste0("t", 1:6))
        )
        
        rna_treat <- matrix(
            rnorm(length(genes_treat) * 6),
            nrow     = length(genes_treat),
            ncol     = 6,
            dimnames = list(genes_treat, paste0("t", 1:6))
        )
        
        make_m2o <- function(
        gene_ids,
        n_feat = 3,
        n_cols = 4
        ) {
            rn <- unlist(
                lapply(
                    gene_ids,
                    function(g) {
                        paste0(g, "_site", seq_len(n_feat))
                    }
                )
            )
            
            matrix(
                rnorm(length(rn) * n_cols),
                nrow     = length(rn),
                ncol     = n_cols,
                dimnames = list(rn, paste0("s", 1:n_cols))
            )
        }
        
        phospho_ctrl <- make_m2o(
            gene_ids = genes_ctrl,
            n_feat   = 3,
            n_cols   = 4
        )
        
        phospho_treat <- make_m2o(
            gene_ids = genes_treat,
            n_feat   = 3,
            n_cols   = 4
        )
        
        data <- list(
            Ctrl = list(
                rna     = rna_ctrl,
                phospho = phospho_ctrl
            ),
            Treat = list(
                rna     = rna_treat,
                phospho = phospho_treat
            )
        )
        
        meta <- data.frame(
            modality       = c("rna", "phospho"),
            many_to_one_k  = c(NA_real_, 3),
            modality_w     = c(1, 1),
            stringsAsFactors = FALSE
        )
        
        k <- 3L
        
        testthat::skip_if_not_installed("uwot")
        testthat::skip_if_not_installed("Matrix")
        testthat::skip_if_not_installed("RSpectra")
        
        res <- cluster_genes_multiomics(
            data         = data,
            meta         = meta,
            k            = k,
            gene_mode    = "union",
            n_neighbors  = 3L,
            verbose      = FALSE
        )
        
        testthat::expect_true(is.list(res))
        
        testthat::expect_true(
            all(
                c(
                    "cluster_table",
                    "centroid_info",
                    "many_to_one_clustering_qc",
                    "umap_fit"
                ) %in% names(res)
            )
        )
        
        testthat::expect_s3_class(
            res$cluster_table,
            "tbl_df"
        )
        
        testthat::expect_true(
            all(
                c("gene", "cluster") %in% colnames(res$cluster_table)
            )
        )
        
        testthat::expect_equal(
            nrow(res$cluster_table),
            length(unique(res$cluster_table$gene))
        )
        
        testthat::expect_true(
            all(res$cluster_table$cluster %in% seq_len(k))
        )
        
        testthat::expect_true(is.list(res$umap_fit))
        
        testthat::expect_true(
            !is.null(res$umap_fit$embedding)
        )
        
        testthat::expect_equal(
            ncol(res$umap_fit$embedding),
            2L
        )
        
        testthat::expect_snapshot_value(
            res$cluster_table,
            style = "deparse"
        )
        
        testthat::expect_snapshot_value(
            res$centroid_info,
            style = "deparse"
        )
    }
)
