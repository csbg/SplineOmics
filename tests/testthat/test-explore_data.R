test_that("explore_data returns correctly structured nested list", {
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

    # Run extract_data() to get the numeric matrix
    extracted_data <- SplineOmics::extract_data(
        data = data,
        feature_name_columns = c("Gene_name"),
        top_row = 4,
        bottom_row = 1165,
        right_col = 37,
        left_col = 2
    )

    # Extract annotation from original raw data (not from extract_data())
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

    report_dir <- tempdir()

    plots <- SplineOmics::explore_data(
        splineomics = splineomics,
        report_dir = report_dir,
        report = FALSE
    )

    expect_type(plots, "list")
    expect_length(plots, 2)

    expect_true(
        all(vapply(plots, function(x) is.list(x) && length(x) == 27, logical(1))),
        info = "Each element of the output should be a list of length 27"
    )
})
