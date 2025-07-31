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
    design = "~ 1 + Phase*Time + Reactor", # best design formula
    mode = "integrated", 
    use_array_weights = FALSE,
    spline_params = list(
      spline_type = c("n"), 
      dof = c(0L) 
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
  
  ## ---- limma_splines_result must exist and have 3 named sub-lists ----
  expect_true(
    "limma_splines_result" %in% names(splineomics),
    info = "Must contain 'limma_splines_result'"
  )
  expect_type(splineomics$limma_splines_result, "list")
  
  expected_subslots <- c(
    "time_effect",
    "avrg_diff_conditions",
    "interaction_condition_time"
  )
  expect_named(
    splineomics$limma_splines_result,
    expected_subslots,
    ignore.order = TRUE
  )
  
  
  ## -------------------------------------------------
  ##  1. Number of significant features (= adj.P < .05)
  ## -------------------------------------------------
  exp_tbl  <- splineomics$limma_splines_result$time_effect$Phase_Exponential
  stat_tbl <- splineomics$limma_splines_result$time_effect$Phase_Stationary
  
  expect_equal(
    sum(exp_tbl$adj.P.Val  < 0.05, na.rm = TRUE),
    61          # ← hard-coded baseline for integrated-mode
  )
  expect_equal(
    sum(stat_tbl$adj.P.Val < 0.05, na.rm = TRUE),
    5
  )
  
  ## ---------------------------------------------------
  ##  2. Ten invariant cells must keep their exact value
  ## ---------------------------------------------------
  get_cell <- function(tbl, row, col) round(tbl[[col]][row], 4)  # 4-dp rounding
  
  ## ---- Phase_Exponential ----
  exp_expected <- c(
    17.3959, 0.0198, 0.0382, 13.4048, 2.0101,
    0.0438, 0.3638, 0.9157,  4.7002, 0.0152
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
    0.7685, 0.0282, 0.0356, 2.6977, 0.4312,
    3.4851, -0.0088, 0.0954, 12.7730, 0.1551
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
  
  ## ---------------------------------------------------
  ##  3. avrg_diff_conditions: 10 invariant cells
  ## ---------------------------------------------------
  avg_tbl <- splineomics$limma_splines_result$avrg_diff_conditions[[1]]
  
  avg_expected <- c(
    0.092, -0.152, 0.292, -0.027, 12.987,
    131.696, 18.748, 0.864, 17.097, 1.115
  )
  avg_actual <- c(
    get_cell(avg_tbl, 2609, 3),
    get_cell(avg_tbl, 4069, 1),
    get_cell(avg_tbl, 2369, 1),
    get_cell(avg_tbl, 1098, 3),
    get_cell(avg_tbl, 1252, 4),
    get_cell(avg_tbl,  634, 5),
    get_cell(avg_tbl, 2097, 5),
    get_cell(avg_tbl, 3911, 5),
    get_cell(avg_tbl,  356, 4),
    get_cell(avg_tbl, 3954, 2)
  )
  expect_equal(
    avg_actual,
    avg_expected,
    tolerance = 1e-4
    )
  
  ## ---------------------------------------------------------------
  ##  4. interaction_condition_time: 10 invariant cells
  ## ---------------------------------------------------------------
  int_tbl <- splineomics$limma_splines_result$interaction_condition_time[[1]]
  
  int_expected <- c(
    0.8370, 0.0302, 0.0412, 2.6055, 0.4598,
    3.3136, -0.0087, 0.1143, 16.3951, -0.2668
  )
  int_actual <- c(
    get_cell(int_tbl,  860, 6),
    get_cell(int_tbl, 3944, 2),
    get_cell(int_tbl,  259, 5),
    get_cell(int_tbl,  481, 4),
    get_cell(int_tbl, 2072, 5),
    get_cell(int_tbl,  299, 4),
    get_cell(int_tbl, 3983, 2),
    get_cell(int_tbl, 2454, 2),
    get_cell(int_tbl, 3720, 3),
    get_cell(int_tbl,  836, 1)
  )
  expect_equal(
    int_actual,
    int_expected,
    tolerance = 1e-4
    )
})