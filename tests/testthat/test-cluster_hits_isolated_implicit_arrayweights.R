test_that("cluster_hits() returns correctly structured result", {
  # Setup
  data <- readRDS(xzfile(system.file(
    "extdata", "proteomics_data.rds.xz", package = "SplineOmics"
  )))
  meta <- read.csv(system.file(
    "extdata", "proteomics_meta.csv", package = "SplineOmics"
  ), stringsAsFactors = FALSE)
  
  first_na_col <- which(is.na(data[1, ]))[1]
  annotation <- data |>
    dplyr::select((first_na_col + 1):ncol(data)) |>
    dplyr::slice(-c(1:3))
  
  extracted_data <- SplineOmics::extract_data(
    data = data,
    feature_name_columns = c("Gene_name"),
    top_row = 4,
    bottom_row = 1165,
    right_col = 37,
    left_col = 2
  )
  
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
    use_array_weights = NULL,
    spline_params = list(
      spline_type = c("n", "n"),
      dof = c(2L, 2L)
    )
  )
  
  splineomics <- SplineOmics::run_limma_splines(splineomics)
  
  # Inputs for clustering
  adj_pthresholds <- c(0.05, 0.05)
  nr_clusters <- list(6, 2)
  gene_column_name <- "Gene_symbol"
  genes <- annotation[[gene_column_name]]
  plot_info <- list(
    y_axis_label = "log2 intensity",
    time_unit = "min",
    treatment_labels = list("feeding"),
    treatment_timepoints = list(0)
  )
  plot_options <- list(
    meta_replicate_column = "Reactor",
    cluster_heatmap_columns = FALSE
  )
  raw_data <- splineomics$data
  
  clustering_results <- SplineOmics::cluster_hits(
    splineomics = splineomics,
    adj_pthresholds = adj_pthresholds,
    nr_clusters = nr_clusters,
    genes = genes,
    plot_info = plot_info,
    plot_options = plot_options,
    raw_data = raw_data,
    report_dir = withr::local_tempdir(),
    report = TRUE,
    max_hit_number = 25,
    adj_pthresh_avrg_diff_conditions = 0.05,
    adj_pthresh_interaction_condition_time = 0.05
  )

  # cluster_hits() returns a list: cluster_table + plots
  expect_type(clustering_results, "list")
  expect_named(
    clustering_results,
    c("cluster_table", "spline_results", "plots"),
    ignore.order = FALSE
  )
  expect_length(clustering_results, 3)
  
  # cluster_table: structure + minimal integrity checks
  cs <- clustering_results$cluster_table
  expect_s3_class(cs, c("tbl_df", "tbl", "data.frame"))
  
  # required metadata cols
  req_cols <- c("feature_nr", "feature_name", "gene")
  expect_true(all(req_cols %in% names(cs)))
  
  # at least two time-effect cluster columns (cat1)
  all_cluster_cols <- grep("^cluster_", names(cs), value = TRUE)
  cond_cols <- setdiff(all_cluster_cols, c("cluster_cat2", "cluster_cat3"))
  expect_gte(length(cond_cols), 2)
  
  # feature_nr is integer-like and unique
  expect_true(is.numeric(cs$feature_nr))
  expect_true(all(cs$feature_nr == as.integer(cs$feature_nr), na.rm = TRUE))
  expect_false(any(duplicated(cs$feature_nr)))
  
  # feature_name is present and non-empty
  expect_true(all(!is.na(cs$feature_name) & nzchar(cs$feature_name)))
  
  # gene column exists but may be all NA (e.g., non-gene omics)
  expect_true("gene" %in% names(cs))
  expect_true(is.character(cs$gene) | all(is.na(cs$gene)))
  
  # Optional: if cat2/cat3 exist, they are character and atomic
  opt_cols <- intersect(c("cluster_cat2", "cluster_cat3"), names(cs))
  for (cc in opt_cols) {
    expect_true(is.atomic(cs[[cc]]))
  }
  
  # plots: present and list-like (count may vary)
  expect_type(clustering_results$plots, "list")
  expect_gte(length(clustering_results$plots), 1)
  
  # Check individual elements
  expect_type(clustering_results$plots, "list")
  expect_length(clustering_results$plots, 18)
  
  testthat::expect_snapshot_value(cs, style = "json2")
})