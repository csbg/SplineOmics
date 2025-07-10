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
  
  splineomics <- SplineOmics::update_splineomics(
    splineomics = splineomics,
    design = "~ 1 + Time + Reactor",
    mode = "isolated", 
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
  exp_tbl  <- splineomics$limma_splines_result$time_effect$Phase_Exponential
  stat_tbl <- splineomics$limma_splines_result$time_effect$Phase_Stationary
  
  expect_equal(
    sum(exp_tbl$adj.P.Val  < 0.05, na.rm = TRUE),
    183                                   # ← hard-coded baseline
  )
  expect_equal(
    sum(stat_tbl$adj.P.Val < 0.05, na.rm = TRUE),
    1
  )
  
  ## ---------------------------------------------------
  ##  2. Ten invariant cells must keep their exact value
  ## ---------------------------------------------------
  get_cell <- function(tbl, row, col) round(tbl[[col]][row], 4)   # ← 4-dp round
  
  ## ---- Phase_Exponential ----
  exp_expected <- c(
    12.4406, -0.0220, -0.0199, 13.2243, 3.1449,
    0.0188,  0.2097,  0.8524,  7.1278, 0.0062
  )
  exp_actual <- c(
    get_cell(exp_tbl, 2609, 3),
    get_cell(exp_tbl, 4069, 1),
    get_cell(exp_tbl, 2369, 1),
    get_cell(exp_tbl, 1098, 3),
    get_cell(exp_tbl, 1252, 4),
    get_cell(exp_tbl,  634, 5),
    get_cell(exp_tbl, 2097, 5),
    get_cell(exp_tbl, 3911, 5),
    get_cell(exp_tbl,  356, 4),
    get_cell(exp_tbl, 3954, 2)
  )
  expect_equal(
    exp_actual,
    exp_expected,
    tolerance = 1e-4
    )
  
  ## ---- Phase_Stationary ----
  stat_expected <- c(
    0.9248, -0.0245, 0.0415, 2.7826, 0.4730,
    3.6694, -0.0196, -0.0426, 16.4171, -0.1241
  )
  stat_actual <- c(
    get_cell(stat_tbl,  860, 6),
    get_cell(stat_tbl, 3944, 2),
    get_cell(stat_tbl,  259, 5),
    get_cell(stat_tbl,  481, 4),
    get_cell(stat_tbl, 2072, 5),
    get_cell(stat_tbl,  299, 4),
    get_cell(stat_tbl, 3983, 2),
    get_cell(stat_tbl, 2454, 2),
    get_cell(stat_tbl, 3720, 3),
    get_cell(stat_tbl,  836, 1)
  )
  expect_equal(
    stat_actual,
    stat_expected,
    tolerance = 1e-4
    )
})