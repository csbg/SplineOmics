library(splinetime)

input_file_path <- here::here("data", "timeseries_proteomics_example.RData")
load(input_file_path) 



# hyperparams screen limma -----------------------------------------------------
data <- as.matrix(data)
data1 <- data
meta1 <- meta

data2 <- data   # data2 is the same as data1, but just to show the functionality
meta2 <- meta

datas <- list(data1, data2)
datas_descr <- c("full_data", "outliers_removed")
metas <- list(meta1, meta2)
designs <- c("~ 1 + Phase*X + Reactor", "~ 1 + X + Reactor")
modes <- c("integrated", "isolated")
condition <- "Phase"
DoFs <- c(2L, 3L, 4L, 5L)
feature_names <- annotation$First.Protein.Description
report_dir <- here::here("results", "hyperparams_screen_reports")
meta_batch_column = "Reactor"
pthresholds <- c(0.05, 0.1)

# Every row a combo to test.
spline_test_configs <- data.frame(spline_type = c("b", "b", "b"),
                                  degrees = c(3L, 3L, 2L),
                                  DoFs = c(NA, NA, NA),
                                  knots = I(list(c(0), c(0, 60), c(0, 60, 120))),
                                  bknots = c(NA, NA, NA))

report_info <- list(
  omics_data_type = "PTX",
  data_description = "Timeseries Proteomics",
  data_collection_date = "2024-05-15",
  analyst_name = "Thomas Rauter",
  project_name = "DGTX",
  contact_info = "thomas.rauter@plus.ac.at"
)

# Test all the different hyperparameters and their combinations.
result <- limma_hyperparams_screen(datas,
                                   datas_descr,
                                   metas,
                                   designs,
                                   modes,
                                   condition,
                                   spline_test_configs,
                                   feature_names,
                                   report_info,
                                   report_dir,
                                   pthresholds,
                                   meta_batch_column)


## Run limma splines -----------------------------------------------------------
DoFs <- c(2L, 2L)
design <- "~ 1 + Phase*X + Reactor"
spline_params = list(spline_type = c("n"),
                     DoFs = c(2L))

# Run the limma spline analysis
result <- run_limma_splines(data, 
                            meta, 
                            design, 
                            spline_params = spline_params,
                            condition,
                            feature_names, 
                            mode = "integrated")

top_tables <- result$top_tables
ttslc_factor_only <- result$ttslc_factor_only
ttslc_factor_time <- result$ttslc_factor_time


## Cluster hits ----------------------------------------------------------------
p_values <- c(0.05, 0.05)
clusters <- list(6L, 3L)
report_dir <- here::here("results", "clustering_reports")

report_info <- list(
  omics_data_type = "PTX",
  data_description = "Timeseries Proteomics",
  data_collection_date = "2024-05-15",
  analyst_name = "Thomas Rauter",
  project_name = "DGTX",
  dataset_name = "CHO Cell Study",
  limma_design = design,
  method_description = "Spline analysis for time-series data.",
  results_summary = "Identified significant changes in the proteins.",
  conclusions = "Found the hits",
  contact_info = "thomas.rauter@plus.ac.at"
)

# Cluster the results from the limma spline analysis.
clustering_results <- cluster_hits(top_tables = top_tables, 
                                   data = data, 
                                   meta = meta, 
                                   condition = condition, 
                                   spline_params = spline_params,
                                   mode = "integrated",
                                   p_values = p_values,
                                   clusters = clusters,
                                   report_info = report_info,
                                   meta_batch_column = meta_batch_column,
                                   report_dir = report_dir
)

clustering_results[[2]]$clustered_hits
