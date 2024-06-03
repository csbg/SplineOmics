# Load libraries ---------------------------------------------------------------
library(tidyverse)
library(readxl)
library(jsonlite)
library(SplineOmics)

# Load and process the data and meta -------------------------------------------
data_excel <- read_excel("/data/data.xlsx")
meta <- read_excel("/data/meta.xlsx")

data <- extract_data(data_excel, feature_names_columns_old)


# Read parameters from JSON file -----------------------------------------------

params <- fromJSON("/workspace/params.json")

# Extract parameters
report_info <- params$report_info
condition <- params$condition
meta_batch_column <- params$meta_batch_column
design <- params$design
spline_params <- params$spline_params
adj_pthresholds <- params$adj_pthresholds
clusters <- params$clusters


# Explore data -----------------------------------------------------------------

report_dir <- "/output/clustering_reports"

plots <- explore_data(data = data,
                      meta = meta,
                      condition = condition,
                      report_info = report_info,
                      meta_batch_column = meta_batch_column,
                      report_dir = report_dir


# Run limma splines -----------------------------------------------------------

# Run the limma spline analysis
result <- run_limma_splines(data, 
                            meta, 
                            design, 
                            spline_params = spline_params,
                            condition)

top_tables <- result$time_effect  

limma_report(run_limma_splines_result = result,
             report_info = report_info,
             report_dir = report_dir)


## Cluster hits ---------------------------------------------------------------- 

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
