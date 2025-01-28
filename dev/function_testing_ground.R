rm(list = ls(all.names = TRUE))

# Setup ------------------------------------------------------------------------
library(devtools)
devtools::load_all()

library(readxl)
library(dplyr)


# Load the data ----------------------------------------------------------------

# data <- read_excel(here::here("dev" ,"data", "PPTX",
#                                     "2 PPTX new filtered (equal to no imputation).xlsx"))

# data <- read_excel(here::here("dev" ,"data", "PPTX",
#                                     "PPTX_processed_no_imputation.xlsx"))


# data_imputed <- read_excel(here::here("dev" ,"data", "PPTX",
#                                    "PPTX_processed_with_imputation.xlsx"))
# 
# meta <- read_excel(here::here("dev" ,"data",
#                               "Time_course_PPTX_old_metadata.xlsx"))

data_excel <- readRDS(system.file(
  "extdata",
  "proteomics_data.rds.xz",
  package = "SplineOmics"
))

meta <- read_excel(here::here(
  "inst",
  "extdata",
  "proteomics_meta.xlsx"
))

annotation <- data_excel %>%
  select(39:ncol(data_excel)) %>%
  dplyr::slice(-c(1:3))



# Get the gene list ------------------------------------------------------------

data_full <- as.data.frame(data_excel)
gene_column_name <- "Gene_symbol"
genes <- data_full[[gene_column_name]][4:nrow(data_full)]




# Automatically extract data matrix from excel/csv table -----------------------

# feature_name_columns <- c("Unique identifier",
#                           "Protein",
#                           "Positions within proteins")

# feature_name_columns <- c("Unique identifier",
#                           "Protein names",
#                           "Positions within proteins")

# feature_name_columns <- c("First.Protein.Description", "Protein.Ids")
feature_name_columns <- c("Gene_name")

# debug(extract_data)
data <- extract_data(
  data_full,
  feature_name_columns,
  user_prompt = FALSE
  )


# Plot the raw data ------------------------------------------------------------

# data <- data[4:3035, 2:35]  # 3036
# data <- as.matrix(data)
# rownames(data) <- rownames(data_imputed)
# # data <- data[1:10, 1:34]
# 
# meta <- meta[!meta$`Sample ID` %in% c(
#   "E10_TP09_Stationary",
#   "E12_TP09_Stationary"
# ), ]
# 
# # Extract first 18 columns from data
# data_exp <- data[, 1:18]
# 
# # Extract the remaining columns from data
# data_stat <- data[, 19:ncol(data)]
# 
# # Extract first 18 rows from meta
# meta_exp <- meta[1:18, ]
# 
# # Extract the remaining rows from meta
# meta_stat <- meta[19:nrow(meta), ]
# 
# make_scatter_plot_html(
#   data_exp,
#   meta_exp,
#   here::here("results", "PPTX_processed_no_imputation_scatter_plots.html"),
#   meta_replicate_column = "Reactor"
# )


# data <- data[4:3035, 2:37]  # 3036
# data <- as.matrix(data)
# rownames(data) <- rownames(data_imputed)
# # data <- data[1:10, 1:34]
# 
# # Extract first 18 columns from data
# data_exp <- data[, 1:18]
# 
# # Extract the remaining columns from data
# data_stat <- data[, 19:ncol(data)]
# 
# # Extract first 18 rows from meta
# meta_exp <- meta[1:18, ]
# 
# # Extract the remaining rows from meta
# meta_stat <- meta[19:nrow(meta), ]
# 
# make_scatter_plot_html(
#   data_exp,
#   meta_exp,
#   here::here("results", "2 PPTX new filtered (equal to no imputation)_scatter_plots.html"),
#   meta_replicate_column = "Reactor"
# )


# Remove outliers --------------------------------------------------------------

data <- data[, !(colnames(data) %in% c(  # Remove potential outliers
  "E12_TP05_Exponential",
  "E10_TP10_Stationary"
)
)]

meta <- meta[!meta$`Sample.ID` %in% c(
  "E12_TP05_Exponential",
  "E10_TP10_Stationary"
), ]


# Simulate RNA-seq data to test voom functionality -----------------------------

# generate_rnaseq_data <- function(n_genes = 1000, n_samples = 36) {
#   set.seed(123)  # For reproducibility
# 
#   # Define sample and gene names
#   gene_names <- paste0("Gene", 1:n_genes)
#   sample_names <- paste0("Sample", 1:n_samples)
# 
#   # Generate random raw RNA-seq counts (Poisson distributed)
#   # Base expression level with some variability
#   base_expression <- rpois(n_genes, lambda = 20)  # Baseline counts
#   counts_matrix <- sapply(1:n_samples, function(x) rpois(n_genes, lambda = base_expression))
# 
#   # Assign row and column names
#   rownames(counts_matrix) <- gene_names
#   colnames(counts_matrix) <- sample_names
# 
#   return(counts_matrix)
# }
# 
# # Example usage:
# n_genes <- 4162  # Adjust the number of genes as needed
# # data <- generate_rnaseq_data(n_genes = n_genes)
# 
# voom_obj <- preprocess_rna_seq_data(
#   raw_counts = data,
#   meta = meta,
#   spline_params = list(spline_type = c("n"),   # Chosen spline parameters
#                        dof = c(2L)),
#   design = "~ 1 + Phase*X + Reactor"
# )

# data <- voom_obj$E

# Explore data -----------------------------------------------------------------

report_info <- list(
  omics_data_type = "PPTX",
  data_description = "Phosphoproteomics data",
  data_collection_date = "Februar 2024",
  analyst_name = "Thomas Rauter",
  contact_info = "thomas.rauter@plus.ac.at",
  project_name = "DGTX"
)

splineomics <- create_splineomics(
  data = data,
  # rna_seq_data = voom_obj,
  meta = meta,
  mode = "integrated",
  annotation = annotation,
  feature_name_columns = feature_name_columns,
  report_info = report_info,
  condition = "Phase",
  meta_batch_column = "Reactor"
)

report_dir <- here::here("results", "explore_data")

# debug(explore_data)
# plots <- explore_data(
#   splineomics,
#   report_dir = report_dir,
#   report = TRUE
#   )


# Prep input to hyperparams screen function ------------------------------------
data1 <- data
meta1 <- meta

data2 <- data[, -c(1, 2)]
meta2 <- meta[-c(1, 2),]
# data2 <- data
# meta2 <- meta

datas <- list(data1, data2)
# rna_seq_datas <- list(voom_obj, voom_obj)  # Just to test it.
datas_descr <- c("full_data", "outliers_removed")
metas <- list(meta1, meta2)
designs <- c("~ 1 + Phase*X + Reactor", "~ 1 + Phase*X + Reactor")
modes <- c("integrated", "integrated")

report_dir <- here::here("results", "hyperparams_screen_reports")

pthresholds <- c(0.05, 0.01)

# Every row a combo to test.
spline_test_configs <- data.frame(
  spline_type = c("n", "n", "b", "b"),
  degree = c(NA, NA, 2L, 4L),
  dof = c(2L, 3L, 3L, 4L)
)


# hyperparams screen limma -----------------------------------------------------
# debug(screen_limma_hyperparams)
# screen_limma_hyperparams(
#   splineomics,
#   datas,
#   datas_descr,
#   metas,
#   designs,
#   modes,
#   spline_test_configs,
#   report_dir,
#   pthresholds,
#   # rna_seq_datas,
#   )


## Run limma splines -----------------------------------------------------------

splineomics <- update_splineomics(
  splineomics = splineomics,
  design = "~ 1 + Phase*Time + Reactor",
  # data = data2,
  # meta = meta2,
  spline_params = list(spline_type = c("n", "n"),   # Chosen spline parameters
                       dof = c(2L, 2L))
)

print(splineomics)


# debug(run_limma_splines)
splineomics <- run_limma_splines(
  splineomics
)


report_dir <- here::here("results", "limma_reports")

# plots <- create_limma_report(
#   splineomics,
#   adj_pthresh = 0.05,
#   report_dir = report_dir
# )


## Cluster hits ----------------------------------------------------------------
adj_pthresholds <- c(0.05, 0.05)   
clusters <- c(6L, 3L)   
report_dir <- here::here("results", "clustering_reports")

plot_info = list(
  y_axis_label = "log2 intensity",
  time_unit = "min",
  # treatment_labels = NA,
  # treatment_timepoints = NA
  treatment_labels = list("feeding"),
  treatment_timepoints = list(0)
)

plot_options = list(
  # meta_replicate_column = "Reactor",
  cluster_heatmap_columns = FALSE
)

# debug(cluster_hits)
clustering_results <- cluster_hits(
  splineomics = splineomics,
  adj_pthresholds = adj_pthresholds,
  clusters = clusters,
  genes = genes,
  plot_info = plot_info,
  plot_options = plot_options,
  report_dir = report_dir,
  report = TRUE,
  # adj_pthresh_avrg_diff_conditions = 0,
  adj_pthresh_interaction_condition_time = 0.25
)


# Perform gsea -----------------------------------------------------------------

gene_set_lib <- c(
  "WikiPathways_2019_Human",
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

report_dir <- here::here("results", "enrichr_databases")

download_enrichr_databases(gene_set_lib, output_dir = report_dir)
 
downloaded_dbs_filepath <- 
  here::here(
    "dev",
    "data",
    "all_databases_08_04_2024-12_41_50.tsv"
  )

databases <- readr::read_tsv(downloaded_dbs_filepath, col_types = readr::cols())

clusterProfiler_params <- list(
  pvalueCutoff = 0.05,
  # adj_p_value = 0.05,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.2
)

report_dir <- here::here("results", "gsea_reports")

# debug(run_gsea)
result <- run_gsea(
  levels_clustered_hits = clustering_results$clustered_hits_levels,
  databases = databases,
  clusterProfiler_params = clusterProfiler_params,
  report_info = report_info,
  report_dir = report_dir
)
