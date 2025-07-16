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
    bottom_row = 4165,
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
    spline_params = list(
      spline_type = c("n", "n"),
      dof = c(2L, 2L)
    )
  )
  
  splineomics <- SplineOmics::run_limma_splines(splineomics)
  
  # Inputs for clustering
  adj_pthresholds <- c(0.05, 0.05)
  nr_clusters <- list(6, 2:3)
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

  # Check top-level structure
  expect_type(clustering_results, "list")
  expect_named(clustering_results,
               c("all_levels_clustering", "plots", "clustered_hits_levels"),
               ignore.order = FALSE
               )
  expect_length(clustering_results, 3)
  
  # Check individual elements
  expect_type(clustering_results$plots, "list")
  expect_length(clustering_results$plots, 10)
  
  expect_type(clustering_results$clustered_hits_levels, "list")
  expect_length(clustering_results$clustered_hits_levels, 2)
  
  expect_type(clustering_results$all_levels_clustering, "list")
  expect_length(clustering_results$all_levels_clustering, 2)
})