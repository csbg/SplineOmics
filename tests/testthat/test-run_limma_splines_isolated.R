test_that("run_limma_splines() works correctly", {
    ## -------------------------
    ##  set-up: reproduce object
    ## -------------------------
    data <- readRDS(xzfile(system.file(
        "extdata",
        "proteomics_data.rds.xz",
        package = "SplineOmics"
    )))

    meta <- read.csv(
        system.file(
            "extdata",
            "proteomics_meta.csv",
            package = "SplineOmics"
        ),
        stringsAsFactors = FALSE
    )

    extracted_data <- SplineOmics::extract_data(
        data = data,
        feature_name_columns = c("Gene_name"),
        top_row = 4,
        bottom_row = 1165,
        right_col = 37,
        left_col = 2
    )

    first_na_col <- which(is.na(data[1, ]))[1]
    annotation <- data |>
        dplyr::select((first_na_col + 1):ncol(data)) |>
        dplyr::slice(-c(1:3))

    report_info <- list(
        omics_data_type = "PTX",
        data_description = "Test data",
        data_collection_date = "2024",
        analyst_name = "Thomas Rauter",
        contact_info = "thomas.rauter@plus.ac.at",
        project_name = "DGTX"
    )

    splineomics <- SplineOmics::create_splineomics(
        data = extracted_data,
        meta = meta,
        annotation = annotation,
        report_info = report_info,
        condition = "Phase",
        meta_batch_column = "Reactor"
    )

    splineomics <- SplineOmics::update_splineomics(
        splineomics = splineomics,
        design = "~ 1 + Time + Reactor",
        mode = "isolated",
        use_array_weights = FALSE,
        spline_params = list(
            spline_type = c("n", "n"),
            dof = c(2L, 0L)
        )
    )

    # Run the function under test
    splineomics <- SplineOmics::run_limma_splines(splineomics)

    ## ---------------------------
    ##  0. Structural sanity
    ## ---------------------------
    expect_true(
        inherits(splineomics, "SplineOmics"),
        info = "Should be S3 class 'SplineOmics'"
    )
    expect_type(splineomics, "list")
    expect_length(splineomics, 18)

    # Check added element
    expect_true(
        "limma_splines_result" %in% names(splineomics),
        info = "Must contain 'limma_splines_result'"
    )
    expect_type(
        splineomics$limma_splines_result,
        "list"
    )
    expect_named(
        splineomics$limma_splines_result,
        "time_effect",
        ignore.order = TRUE
    )

    ## -------------------------------------------------
    ##  1. Number of significant features (= adj.P < .05)
    ## -------------------------------------------------
    exp_tbl <- splineomics$limma_splines_result$time_effect$Phase_Exponential
    stat_tbl <- splineomics$limma_splines_result$time_effect$Phase_Stationary

    expect_equal(
        sum(exp_tbl$adj.P.Val < 0.05, na.rm = TRUE),
        32 # ← hard-coded baseline
    )
    expect_equal(
        sum(stat_tbl$adj.P.Val < 0.05, na.rm = TRUE),
        1
    )

    testthat::expect_snapshot_value(
        list(Phase_Exponential = exp_tbl, Phase_Stationary = stat_tbl),
        style = "deparse"
    )
})
