test_that("extract_data returns correct numeric matrix", {
    data <- readRDS(xzfile(system.file(
        "extdata",
        "proteomics_data.rds.xz",
        package = "SplineOmics"
    )))

    result <- SplineOmics::extract_data(
        data = data,
        feature_name_columns = c("Gene_name"),
        top_row = 4,
        bottom_row = 1165,
        right_col = 37,
        left_col = 2
    )

    expect_true(
        is.matrix(result),
        info = "Result should be a matrix"
    )
    expect_true(
        is.numeric(result),
        info = "Matrix should be numeric (int or double only)"
    )
    expect_equal(
        dim(result),
        c(1162, 36),
        info = "Matrix should have 1162 rows and 36 columns"
    )
})
