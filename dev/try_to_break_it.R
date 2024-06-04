rm(list = ls(all.names = TRUE))


# Setup ------------------------------------------------------------------------
library(devtools)
devtools::load_all()

library(conflicted)
library(strict)



# Load the data ----------------------------------------------------------------
load(here::here("test_data", "timeseries_proteomics_example.RData")) 

library(tidyverse)
library(readxl)
data_excel <- read_excel(here::here("inst", "README_data",
                                    "PPTX_processed_with_imputation.xlsx"))


# Automatically extract data matrix from excel/csv table

# debug(extract_data)
data_extract <- extract_data(data_excel,
                             "Unique identifier")


# Explore data -----------------------------------------------------------------

report_info <- list(
  omics_data_type = "PPTX",
  data_description = "Old phosphoproteomics data with the missing two samples",
  data_collection_date = "February 2024",
  analyst_name = "Thomas Rauter",
  contact_info = "thomas.rauter@plus.ac.at",
  project_name = "DGTX")

# condition <- "Phase"
condition <- "Phase"
meta_batch_column <- "Reactor"
# meta_batch2_column <- "Reactor"
report_dir <- here::here("results", "explore_data")

data_with_errors <- data
meta_with_errors <- meta

data_with_errors$new_col <- rep(1, nrow(data))
meta_with_errors <- rbind(meta_with_errors, data.frame(Sample = 10, Reactor = "R3", Time.Point = "TP13",
                                                       Phase = "Stationary", Time = 240))

# debug(explore_data)
plots <- explore_data(data = data_with_errors,
                      meta = meta_with_errors,
                      condition = condition,
                      report_info = report_info,
                      meta_batch_column = meta_batch_column,
                      # meta_batch2_column = meta_batch2_column,
                      report_dir = report_dir)


# Prep input to hyperparams screen function ------------------------------------
data1 <- data
meta1 <- meta

data2 <- data[, -c(1, 2)]
meta2 <- meta[-c(1, 2),]

datas <- list(data1, data2)
datas_descr <- c("full_data", "outliers_removed")
metas <- list(meta1, meta2)
designs <- c("~ 1 + Phase*X + Reactor", "~ 1 + X + Reactor")
condition <- "Phase"
report_dir <- here::here("results", "hyperparams_screen_reports")
meta_batch_column = "Reactor"
pthresholds <- c(0.05, 0.1)

# Every row a combo to test.
spline_test_configs <- data.frame(spline_type = c("n", "n", "n", "n"),
                                  degree = c(NA, NA, NA, NA),
                                  dof = c(2L, 3L, 4L, 5L),
                                  knots = I(list(c(NA), c(NA), c(NA), c(NA))),
                                  bknots = I(list(c(NA), c(NA), c(NA), c(NA))))


# hyperparams screen limma -----------------------------------------------------
# debug(limma_hyperparams_screen)
# result <- limma_hyperparams_screen(datas,
#                                    datas_descr,
#                                    metas,
#                                    designs,
#                                    condition,
#                                    spline_test_configs,
#                                    report_info,
#                                    report_dir,
#                                    pthresholds,
#                                    meta_batch_column)


## Run limma splines -----------------------------------------------------------
design <- "~ 1 + Phase*X + Reactor"          # Chosen limma design
# design <- "~ 1 + X + Reactor"

spline_params = list(spline_type = c("n"),   # Chosen spline parameters
                     dof = c(2L))

# Run the limma spline analysis
result <- run_limma_splines(data, 
                            meta, 
                            design, 
                            spline_params = spline_params,
                            condition)

top_tables <- result$avrg_diff_conditions 

report_dir <- here::here("results", "limma_reports")

limma_report(run_limma_splines_result = result,
             report_info = report_info,
             report_dir = report_dir)

## Cluster hits ----------------------------------------------------------------
adj_pthresholds <- c(0.05, 0.05)   
clusters <- list(2L, 2L)   
report_dir <- here::here("results", "clustering_reports")

# debug(cluster_hits)
clustering_results <- cluster_hits(top_tables = top_tables, 
                                   data = data, 
                                   meta = meta, 
                                   design = design,
                                   condition = condition, 
                                   spline_params = spline_params,
                                   adj_pthresholds = adj_pthresholds,
                                   clusters = clusters,
                                   report_info = report_info,
                                   meta_batch_column = meta_batch_column,
                                   # meta_batch2_column = meta_batch2_column,
                                   report_dir = report_dir,
                                   report = TRUE)

clustering_results$all_levels_clustering[[2]]$clustered_hits