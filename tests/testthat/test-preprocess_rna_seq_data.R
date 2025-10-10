library(testthat)
library(SplineOmics)

test_that("splineomics object structure is as expected", {
    # Load data as in your example
    data <- readRDS(xzfile(system.file(
        "extdata",
        "rna_seq_data.rds.xz",
        package = "SplineOmics"
    )))
    data <- data[rowSums(data) > 0, ]

    meta <- readr::read_csv(
        system.file(
            "extdata",
            "rna_seq_meta.csv",
            package = "SplineOmics"
        ),
        show_col_types = FALSE
    )

    report_info <- list(
        omics_data_type = "RNA",
        data_description = "RNA-seq data of CHO cells",
        data_collection_date = "December 2024",
        analyst_name = "Thomas Rauter",
        contact_info = "thomas.rauter@plus.ac.at",
        project_name = "DGTX"
    )

    splineomics <- SplineOmics::create_splineomics(
        data = data,
        meta = meta,
        report_info = report_info,
        condition = "Condition",
        meta_batch_column = "Plate"
    )

    splineomics <- SplineOmics::update_splineomics(
        splineomics = splineomics,
        data = data,
        design = "~ 1 + Condition*Time + Plate + (1|Reactor)",
        use_array_weights = FALSE,
        mode = "integrated",
        spline_params = list(
            spline_type = c("n"),
            dof = c(3L)
        ),
        bp_cfg = c(
            n_cores = 1,
            blas_threads = 1
        )
    )

    splineomics <- SplineOmics::preprocess_rna_seq_data(splineomics)

    # Now test its structure:

    testthat::expect_s3_class(splineomics, "SplineOmics")
    testthat::expect_type(splineomics, "list")
    testthat::expect_length(splineomics, 16)

    testthat::expect_true(is.matrix(splineomics$data))
    testthat::expect_equal(dim(splineomics$data)[2], 136)

    testthat::expect_true(is.list(splineomics$rna_seq_data))
    testthat::expect_length(splineomics$rna_seq_data, 5)

    testthat::expect_s3_class(splineomics$meta, "spec_tbl_df")
    testthat::expect_equal(ncol(splineomics$meta), 5)
    testthat::expect_equal(nrow(splineomics$meta), 136)

    testthat::expect_type(splineomics$condition, "character")
    testthat::expect_equal(length(splineomics$condition), 1)

    testthat::expect_type(splineomics$report_info, "list")
    testthat::expect_length(splineomics$report_info, 6)

    testthat::expect_type(splineomics$meta_batch_column, "character")
    testthat::expect_equal(length(splineomics$meta_batch_column), 1)

    testthat::expect_type(splineomics$design, "character")

    testthat::expect_type(splineomics$use_array_weights, "logical")
    testthat::expect_false(splineomics$use_array_weights)

    testthat::expect_type(splineomics$mode, "character")
    testthat::expect_equal(splineomics$mode, "integrated")

    testthat::expect_type(splineomics$spline_params, "list")
    testthat::expect_length(splineomics$spline_params, 2)

    testthat::expect_type(splineomics$padjust_method, "character")
    testthat::expect_equal(splineomics$padjust_method, "BH")

    testthat::expect_type(splineomics$bp_cfg, "double")
    testthat::expect_equal(length(splineomics$bp_cfg), 2)
    testthat::expect_named(splineomics$bp_cfg, c("n_cores", "blas_threads"))
    testthat::expect_equal(as.numeric(splineomics$bp_cfg), c(1, 1))
})
