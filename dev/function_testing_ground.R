rm(list = ls(all.names = TRUE))


# Setup ------------------------------------------------------------------------
library(devtools)
devtools::load_all()

library(conflicted)


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

condition <- "Phase"
meta_batch_column <- "Reactor"
report_dir <- here::here("results", "explore_data")

# debug(explore_data)
# plots <- explore_data(data = data,
#                       meta = meta,
#                       condition = condition,
#                       report_info = report_info,
#                       meta_batch_column = meta_batch_column,
#                       # meta_batch2_column = NULL,
#                       report_dir = report_dir)


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

top_tables1 <- result$time_effect
top_tables2 <- result$avrg_diff_conditions 
top_tables3 <- result$interaction_condition_time

report_dir <- here::here("results", "limma_reports")

# limma_report(run_limma_splines_result = result,
#              report_info = report_info,
#              report_dir = report_dir)

## Cluster hits ----------------------------------------------------------------
adj_pthresholds <- c(0.05, 0.05)   
clusters <- list(2L, 2L)   
report_dir <- here::here("results", "clustering_reports")

combo_list <- list(top_tables1, top_tables1)

#debug(cluster_hits)
clustering_results <- cluster_hits(top_tables = top_tables2, 
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

clustered_genes_exp <- 
  clustering_results$all_levels_clustering[[1]]$clustered_hits



# Perform gsea -----------------------------------------------------------------

gene_set_lib <- c("WikiPathways_2019_Human",
                  "NCI-Nature_2016",
                  "TRRUST_Transcription_Factors_2019",
                  "MSigDB_Hallmark_2020",
                  "GO_Cellular_Component_2018",
                  "CORUM",
                  "KEGG_2019_Human",
                  "TRANSFAC_and_JASPAR_PWMs",
                  "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
                  "GO_Biological_Process_2018",
                  "GO_Molecular_Function_2018",
                  "Human_Gene_Atlas"
)

download_enrichr_databases(gene_set_lib)

downloaded_dbs_filepath <- 
  here::here("data", "all_databases_08_04_2024-12_41_50.tsv")

clusterProfiler_params <- list(adj_p_value = 0.05,
                               pAdjustMethod = "BH",
                               minGSSize = 10,
                               maxGSSize = 500,
                               qvalueCutoff = 0.2)

result <- run_gsea(clustered_genes = clustered_genes_exp,
                   downloaded_dbs_filepath = downloaded_dbs_filepath,
                   params = clusterProfiler_params,)
















# Further section --------------------------------------------------------------

# plot exp stat diff ----
row_header <- 7011

# Extract the row and convert it to numeric vector
row_data <- as.numeric(data[row_header, ])

# Split the row into two series of 17 values each
series1 <- row_data[1:17]
series2 <- row_data[18:34]

# Get the first 17 Time values from 'meta'
time_values <- meta$Time[1:17]

# Combine both series into a single dataframe for plotting
plot_data <- data.frame(
  Time = c(time_values, time_values), # Duplicate Time for both series
  Value = c(series1, series2),
  Series = factor(c(rep("Series 1", 17), rep("Series 2", 17)))
)

# Plot without connecting lines
p <- ggplot(plot_data, aes(x = Time, y = Value, color = Series)) +
  geom_point() + # Only plot points, no connecting lines
  theme_minimal() +
  labs(title = paste("Comparison of Series for Row", row_header),
       x = "Time",
       y = "Value",
       color = "Series")

print(p)

# Plot Series 1 individually
p_series1 <- ggplot(subset(plot_data, Series == "Series 1"), aes(x = Time, y = Value)) +
  geom_point(color = "blue") +
  theme_minimal() +
  labs(title = "Series 1 Individual Plot",
       x = "Time",
       y = "Value")

print(p_series1) # Print Series 1 plot

# Plot Series 2 individually
p_series2 <- ggplot(subset(plot_data, Series == "Series 2"), aes(x = Time, y = Value)) +
  geom_point(color = "red") +
  theme_minimal() +
  labs(title = "Series 2 Individual Plot",
       x = "Time",
       y = "Value")

print(p_series2) # Print Series 2 plot



# Special results generation section (delete afterwards again) ----

## -------------------- Generate the final tables -----------------------------
library(WriteXLS)


ttlc <- ttslc[[1]]
sorted_indices <- ttlc$feature_index
sorted_annotation <- annotation
sorted_annotation <- sorted_annotation[sorted_indices, ]

desired_columns <- sorted_annotation[, c("Protein.Ids", "Genes", 
                                         "First.Protein.Description")]
ttlc_mod <- cbind(feature_index = sorted_indices, 
                                     desired_columns, ttlc)




## ---------------------- Generate the Excel file -------------------------------
timestamp <- format(Sys.time(), "%d_%m_%Y-%H_%M_%S")
filename <- paste0("limma_splines_PTX_exp_stat_diff_", timestamp, ".xlsx")
filepath <- here::here("results", filename)
dir_path <- dirname(filepath)
if (!dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)
}

WriteXLS(x = ttlc_mod, 
         ExcelFileName = filepath,  
         AdjWidth = TRUE, 
         BoldHeaderRow = TRUE, 
         FreezeRow = 1, 
         SheetNames = "exp_stat_diff",  
         na = "N/A")
