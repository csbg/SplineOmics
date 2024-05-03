library(splinetime)

input_file_path <- here::here("data", "timeseries_proteomics_example.RData")
load(input_file_path) 


# hyperparams screen limma -----------------------------------------------------
data <- as.matrix(data)
datas <- list(data)
datas_descr <- c("full_data")
metas <- list(meta, )
designs <- c("~ 1 + Phase*X + Reactor", "~ 1 + X + Reactor")
modes <- c("integrated", "isolated")
condition <- "Phase"
DoFs <- c(2L, 3L, 4L, 5L)
feature_names <- annotation$First.Protein.Description
report_dir <- here::here("results", "hyperparams_screening")
pthresholds <- c(0.05, 0.1)

spline_configs = list(spline_type = c("n", "n", "n"), 
                      degrees = c(3L, 3L, 3L),
                      knots = list(c(0), c(0, 60), c(0, 60, 120)))

result <- limma_hyperparams_screen(datas,
                                   datas_descr,
                                   metas,
                                   designs,
                                   modes,
                                   condition,
                                   spline_configs,
                                   feature_names,
                                   report_dir,
                                   pthresholds)


## Run limma splines -----------------------------------------------------------
DoFs <- c(2L, 2L)
design <- "~ 1 + Phase*X + Reactor"

spline_params = list(spline_type = c("n"),
                     DoFs = c(2L))

result <- run_limma_splines(data, 
                            meta, 
                            design, 
                            spline_params,
                            condition,
                            feature_names, 
                            "integrated")

top_tables <- result$top_tables
ttslc_factor_only <- result$ttslc_factor_only
ttslc_factor_time <- result$ttslc_factor_time


## Cluster hits ----------------------------------------------------------------
p_values <- c(0.05, 0.05)
clusters <- list(6L, 3L)
report_dir <- here::here("results", "clustering")
data <- removeBatchEffect(x = data, batch = meta$Reactor)

clustering_results <- cluster_hits(top_tables, 
                                   data, 
                                   meta, 
                                   condition, 
                                   p_values, 
                                   clusters, 
                                   report_dir)

clustering_results[[2]]$clustered_hits