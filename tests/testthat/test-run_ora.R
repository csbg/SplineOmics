testthat::test_that("run_ora example runs without error and has stable structure", {
    testthat::skip_if_not_installed("clusterProfiler")
    testthat::skip_if_not_installed("tibble")
    
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
    
    # run ORA: must not error
    expect_no_error({
        res <- run_ora(
            cluster_table            = cluster_table,
            databases                = databases,
            clusterProfiler_params   = clusterProfiler_params
        )
    })
    
    # snapshot on the high-level structure
    expect_snapshot(names(res))
})
