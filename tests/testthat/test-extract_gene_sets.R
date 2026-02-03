testthat::test_that(
    "extract_gene_sets errors if organism db is missing",
    {
        testthat::expect_error(
            extract_gene_sets(
                organism_db = "org.This.Does.Not.Exist",
                output_dir = tempdir(),
                filename = "x.tsv"
            ),
            "Organism database 'org.This.Does.Not.Exist' is not installed\\."
        )
    }
)

testthat::test_that(
    "extract_gene_sets runs and writes TSV when deps are available",
    {
        if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
            testthat::skip(
                "AnnotationDbi not installed (Suggested); skipping test."
            )
        }
        
        org_candidates <- c(
            "org.Hs.eg.db",
            "org.Mm.eg.db"
        )
        
        org_available <- org_candidates[
            vapply(
                org_candidates,
                requireNamespace,
                logical(1),
                quietly = TRUE
            )
        ]
        
        if (length(org_available) == 0L) {
            testthat::skip(
                "No org.*.eg.db package installed; skipping test."
            )
        }
        
        outdir <- withr::local_tempdir()
        filename <- "genesets_test.tsv"
        
        result <- suppressMessages(
            extract_gene_sets(
                organism_db = org_available[[1]],
                output_dir = outdir,
                filename = filename
            )
        )
        
        testthat::expect_s3_class(result, "data.frame")
        
        testthat::expect_true(
            all(
                c("DB", "Geneset", "Gene") %in% colnames(result)
            )
        )
        
        testthat::expect_true(nrow(result) > 0L)
        testthat::expect_false(anyNA(result$DB))
        testthat::expect_false(anyNA(result$Geneset))
        testthat::expect_false(anyNA(result$Gene))
        
        testthat::expect_true(
            all(
                result$DB %in% c("GO_BP", "GO_MF", "GO_CC", "KEGG")
            )
        )
        
        testthat::expect_identical(
            nrow(result),
            nrow(
                unique(
                    result[, c("DB", "Geneset", "Gene"), drop = FALSE]
                )
            )
        )
        
        filepath <- file.path(outdir, filename)
        
        testthat::expect_true(file.exists(filepath))
        testthat::expect_gt(file.info(filepath)$size, 0)
        
        roundtrip <- utils::read.table(
            filepath,
            sep = "\t",
            header = TRUE,
            quote = "",
            stringsAsFactors = FALSE,
            check.names = FALSE
        )
        
        testthat::expect_true(
            all(
                c("DB", "Geneset", "Gene") %in% colnames(roundtrip)
            )
        )
        
        key_result <- paste(
            result$DB,
            result$Geneset,
            result$Gene,
            sep = "\r"
        )
        
        key_roundtrip <- paste(
            roundtrip$DB,
            roundtrip$Geneset,
            roundtrip$Gene,
            sep = "\r"
        )
        
        testthat::expect_setequal(
            key_result,
            key_roundtrip
        )
    }
)
