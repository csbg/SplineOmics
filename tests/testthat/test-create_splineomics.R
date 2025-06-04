test_that("create_splineomics returns valid SplineOmics object", {
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
    data = data,
    meta = meta,
    annotation = annotation,
    report_info = report_info,
    condition = "Phase",
    meta_batch_column = "Reactor"
  )
  
  expect_true(
    inherits(splineomics, "SplineOmics"),
    info = "Should be S3 class 'SplineOmics'"
    )
  expect_type(splineomics, "list")
  expect_length(splineomics, 15)
})
