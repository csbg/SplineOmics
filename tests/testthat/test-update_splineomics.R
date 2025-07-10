test_that("update_splineomics returns valid SplineOmics object", {
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
    bottom_row = 4165,
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
  
  updated_splineomics <- SplineOmics::update_splineomics(
    splineomics = splineomics,
    design = "~ 1 + Time + Reactor",
    mode = "isolated",
    spline_params = list(
      spline_type = c("n", "n"),
      dof = c(2L, 2L)
    )
  )
  
  expect_true(
    inherits(updated_splineomics, "SplineOmics"),
    info = "Should be S3 class 'SplineOmics'"
    )
  expect_type(
    updated_splineomics,
    "list"
    )
  expect_length(
    updated_splineomics,
    16
    )
})