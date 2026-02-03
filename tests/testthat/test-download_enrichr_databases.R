testthat::test_that("download_enrichr_databases returns a tidy table", {
    testthat::local_mocked_bindings(
        enrichr_get_genesets = function(databases) {
            list(
                WikiPathways_2019_Human = list(
                    "Pathway A" = c("TP53", "BRCA1,extra"),
                    "Pathway B" = c("EGFR")
                )
            )
        }
    )
    
    out <- download_enrichr_databases(
        gene_set_lib = "WikiPathways_2019_Human",
        output_dir = NULL
    )
    
    testthat::expect_s3_class(out, "data.frame")
    testthat::expect_true(all(c("DB", "Geneset", "Gene") %in% names(out)))
    testthat::expect_true(nrow(out) == 3L)
    testthat::expect_true(all(out$DB == "WikiPathways_2019_Human"))
    testthat::expect_true(any(out$Geneset == "Pathway A"))
    testthat::expect_false(any(grepl(",", out$Gene)))
    testthat::expect_true(any(out$Gene == "BRCA1"))
})

testthat::test_that("download_enrichr_databases does not write by default", {
    testthat::local_mocked_bindings(
        enrichr_get_genesets = function(databases) {
            list(
                X = list(
                    Y = c("A")
                )
            )
        }
    )
    
    td <- withr::local_tempdir()
    
    out <- download_enrichr_databases(
        gene_set_lib = "X",
        output_dir = NULL,
        filename = "should_not_exist.tsv"
    )
    
    testthat::expect_true(nrow(out) == 1L)
    testthat::expect_false(file.exists(file.path(td, "should_not_exist.tsv")))
})

testthat::test_that("download_enrichr_databases writes a tsv if requested", {
    testthat::local_mocked_bindings(
        enrichr_get_genesets = function(databases) {
            list(
                Lib1 = list(
                    Set1 = c("G1", "G2")
                )
            )
        }
    )
    
    td <- withr::local_tempdir()
    fn <- "enrichr_test.tsv"
    fp <- file.path(td, fn)
    
    out <- download_enrichr_databases(
        gene_set_lib = "Lib1",
        output_dir = td,
        filename = fn
    )
    
    testthat::expect_true(file.exists(fp))
    disk <- utils::read.table(
        fp,
        sep = "\t",
        header = TRUE,
        stringsAsFactors = FALSE,
        quote = "",
        comment.char = ""
    )
    
    testthat::expect_true(nrow(out) == 2L)
    testthat::expect_true(nrow(disk) == 2L)
    testthat::expect_true(all(names(disk) == c("DB", "Geneset", "Gene")))
})

testthat::test_that("download_enrichr_databases errors if helper errors", {
    testthat::local_mocked_bindings(
        enrichr_get_genesets = function(databases) {
            stop("network down")
        }
    )
    
    testthat::expect_error(
        download_enrichr_databases(
            gene_set_lib = "Lib1",
            output_dir = NULL
        ),
        "Failed to download Enrichr gene set libraries"
    )
})

testthat::test_that("download_enrichr_databases errors if nothing retrieved", {
    testthat::local_mocked_bindings(
        enrichr_get_genesets = function(databases) {
            list()
        }
    )
    
    testthat::expect_error(
        download_enrichr_databases(
            gene_set_lib = "Lib1",
            output_dir = NULL
        ),
        "No Enrichr gene sets were retrieved"
    )
})
