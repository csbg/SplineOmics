library(testthat)
library(SplineOmics)

test_that("phosphoproteomics pipeline snapshot", {
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
      Stationary  = "feeding"
    ),
    treatment_timepoints = list(
      Exponential = 0,
      Stationary  = 0
    )
  )
  
  # run and collect warnings
  warning_message <- capture_warnings({
    excursion_plots <- SplineOmics::find_pvc(
      splineomics = splineomics,
      alphas = 0.025,
      plot_info = plot_info,
      report_dir = withr::local_tempdir()
    )
  })
  
  # light sanity checks that won’t churn snapshots
  expect_type(excursion_plots, "list")
  expect_setequal(names(excursion_plots), c("Exponential", "Stationary"))
  
  exp  <- excursion_plots$Exponential
  stat <- excursion_plots$Stationary
  
  expect_type(exp,  "list")
  expect_type(stat, "list")
  
  # prepare a compact, stable snapshot payload
  # round numerics to reduce noise across R/BLAS platforms
  round6 <- function(x) {
    if (is.numeric(x)) return(round(x, 6))
    if (is.matrix(x))  return(apply(x, 2, function(col) round(col, 6)))
    x
  }
  
  snapshot_payload <- list(
    warnings = warning_message,
    structure = list(
      exp_names  = names(exp),
      stat_names = names(stat),
      exp_plots_len  = length(exp$plots),
      stat_plots_len = length(stat$plots)
    ),
    exponential = list(
      pvc_adj_pvals_dim   = dim(exp$pvc_adj_pvals),
      pvc_adj_pvals_head5 = round6(exp$pvc_adj_pvals[seq_len(min(5, nrow(exp$pvc_adj_pvals))), , drop = FALSE]),
      pvc_pattern_summary = exp$pvc_pattern_summary
    ),
    stationary = list(
      pvc_adj_pvals_dim   = dim(stat$pvc_adj_pvals),
      pvc_adj_pvals_head5 = round6(stat$pvc_adj_pvals[seq_len(min(5, nrow(stat$pvc_adj_pvals))), , drop = FALSE]),
      pvc_pattern_summary = stat$pvc_pattern_summary
    )
  )
  
  # store as JSON-like snapshot (stable + human-readable)
  expect_snapshot_value(snapshot_payload, style = "json2")
})
