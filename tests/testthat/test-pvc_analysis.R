library(testthat)
library(SplineOmics)

test_that("phosphoproteomics pipeline runs without errors", {
  data <- readRDS(xzfile(system.file(
    "extdata",
    "phosphoproteomics_data.rds.xz",
    package = "SplineOmics"
  )))
  
  meta <- read.csv(
    system.file(
      "extdata",
      "phosphoproteomics_meta.csv",
      package = "SplineOmics"
    ),
    stringsAsFactors = FALSE
  )
  
  first_na_col <- which(is.na(data[1, ]))[1]
  annotation <- data |> dplyr::select((first_na_col + 1):ncol(data))
  
  data <- SplineOmics::extract_data(
    data = data,
    feature_name_columns = c("T..Gene"),
    use_row_index = TRUE,
    top_row = 1,
    bottom_row = 1000,
    right_col = 36,
    left_col = 1
  )
  
  report_info <- list(
    omics_data_type = "PPTX",
    data_description = "Phosphoproteomics data of CHO cells",
    data_collection_date = "February 2024",
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
  
  plot_info <- list(
    y_axis_label = "log2 intensity",
    time_unit = "min",
    treatment_labels = list(
      Exponential = "feeding",
      Stationary = "feeding"
    ),
    treatment_timepoints = list(
      Exponential = 0,
      Stationary = 0
    )
  )
  
  warning_message <- capture_warnings({
    excursion_plots <- SplineOmics::find_pvc(
      splineomics = splineomics,
      alphas = 0.025,
      plot_info = plot_info,
      report_dir = withr::local_tempdir()
    )
  })
  
  expect_true(any(grepl("Partial NA coefficients", warning_message)))
  
  # Check overall structure
  expect_type(excursion_plots, "list")
  expect_named(excursion_plots, c("Exponential", "Stationary"))
  expect_length(excursion_plots, 2)
  
  # Check Exponential
  exp <- excursion_plots$Exponential
  expect_type(exp, "list")
  expect_named(exp, c("plots", "pvc_adj_pvals", "pvc_pattern_summary"))
  expect_length(exp$plots, 1)
  expect_type(exp$plots, "list")
  
  expect_type(exp$pvc_adj_pvals, "double")
  expect_equal(dim(exp$pvc_adj_pvals), c(1000, 4))
  expect_equal(exp$pvc_adj_pvals[1, 1], 0.99561, tolerance = 1e-5)
  expect_equal(exp$pvc_adj_pvals[1, 2], 0.94922, tolerance = 1e-5)
  expect_equal(exp$pvc_adj_pvals[1, 3], 0.966554, tolerance = 1e-5)
  expect_equal(exp$pvc_adj_pvals[1, 4], 0.9926, tolerance = 1e-5)
  
  # pvc_pattern_summary for Exponential
  exp_summary <- exp$pvc_pattern_summary
  expect_s3_class(exp_summary, "data.frame")
  expect_equal(dim(exp_summary), c(4, 6))
  expect_equal(names(exp_summary), c("-60", "15", "60", "90", "120", "240"))
  expect_equal(as.integer(unlist(exp_summary[1, ])), c(0, 0, 0, 0, 0, 0))
  expect_equal(as.integer(unlist(exp_summary[2, ])), c(0, 0, 0, 1, 0, 0))
  expect_equal(as.integer(unlist(exp_summary[3, ])), c(0, 0, 0, 0, 0, 0))
  expect_equal(as.integer(unlist(exp_summary[4, ])), c(0, 0, 0, 0, 0, 0))
  
  # Check Stationary
  stat <- excursion_plots$Stationary
  expect_type(stat, "list")
  expect_named(stat, c("plots", "pvc_adj_pvals", "pvc_pattern_summary"))
  expect_length(stat$plots, 11)
  expect_type(stat$plots, "list")
  
  expect_type(stat$pvc_adj_pvals, "double")
  expect_equal(dim(stat$pvc_adj_pvals), c(1000, 4))
  expect_equal(stat$pvc_adj_pvals[1, 1], 0.621858, tolerance = 1e-5)
  expect_equal(stat$pvc_adj_pvals[1, 2], 0.025854, tolerance = 1e-5)
  expect_equal(stat$pvc_adj_pvals[1, 3], 0.357782, tolerance = 1e-5)
  expect_equal(stat$pvc_adj_pvals[1, 4], 0.99784, tolerance = 1e-5)
  
  # pvc_pattern_summary for Stationary
  stat_summary <- stat$pvc_pattern_summary
  expect_s3_class(stat_summary, "data.frame")
  expect_equal(dim(stat_summary), c(4, 6))
  expect_equal(names(stat_summary), c("-60", "15", "60", "90", "120", "240"))
  expect_equal(as.integer(unlist(stat_summary[1, ])), c(0, 0, 8, 0, 0, 0))
  expect_equal(as.integer(unlist(stat_summary[2, ])), c(0, 1, 2, 0, 0, 0))
  expect_equal(as.integer(unlist(stat_summary[3, ]))[1:5], c(0, 1, 0, 0, 0))
  expect_equal(as.integer(unlist(stat_summary[4, ])), c(0, 0, 0, 0, 0, 0))
})
