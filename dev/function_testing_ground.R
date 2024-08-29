
rm(list = ls(all.names = TRUE))

# Setup ------------------------------------------------------------------------
library(devtools)
devtools::load_all()

library(readxl)
library(dplyr)

# interactive_demo()


# Load the data ----------------------------------------------------------------

# data_excel <- read_excel(here::here("dev" ,"data", "PTX",
#                                     "PTX_processed_table.xlsx"))
# 
# meta <- read_excel(here::here("dev" ,"data",
#                               "Time_course_PTX_metadata.xlsx"))


data_excel <- readRDS(system.file(
  "extdata",
  "proteomics_data.rds",
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

# feature_name_columns <- c("First.Protein.Description", "Protein.Ids")
feature_name_columns <- c("Gene_name")

# debug(extract_data)
data <- extract_data(
  data_excel,
  feature_name_columns,
  user_prompt = FALSE
  )


# Remove outliers --------------------------------------------------------------

# data <- data[, !(colnames(data) %in% c(  # Remove potential outliers
#   "E12_TP05_Exponential",   
#   "E10_TP10_Stationary"
# )
# )]
# 
# meta <- meta[!meta$`Sample.ID` %in% c(
#   "E12_TP05_Exponential", 
#   "E10_TP10_Stationary"
# ), ]


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
# n_genes <- 7162  # Adjust the number of genes as needed
# data <- generate_rnaseq_data(n_genes = n_genes)

# Explore data -----------------------------------------------------------------

report_info <- list(
  omics_data_type = "PTX",
  data_description = "Proteomics data",
  data_collection_date = "February 2024",
  analyst_name = "Thomas Rauter",
  contact_info = "thomas.rauter@plus.ac.at",
  project_name = "DGTX"
)

splineomics <- create_splineomics(
  data = data,
  meta = meta,
  annotation = annotation,
  feature_name_columns = feature_name_columns,
  report_info = report_info,
  condition = "Phase",
  meta_batch_column = "Reactor",
  preprocess_rna_seq = FALSE
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

datas <- list(data1, data2)
datas_descr <- c("full_data", "outliers_removed")
metas <- list(meta1, meta2)
designs <- c("~ 1 + Phase*X + Reactor", "~ 1 + X + Reactor")

report_dir <- here::here("results", "hyperparams_screen_reports")

pthresholds <- c(0.05, 0.1)

# Every row a combo to test.
spline_test_configs <- data.frame(
  spline_type = c("n", "n", "n", "n"),
  degree = c(NA, NA, NA, NA),
  dof = c(2L, 3L, 4L, 5L),
  knots = I(list(c(NA), c(NA), c(NA), c(NA))),
  bknots = I(list(c(NA), c(NA), c(NA), c(NA)))
)


# hyperparams screen limma -----------------------------------------------------
# debug(screen_limma_hyperparams)
# screen_limma_hyperparams(
#   splineomics,
#   datas,
#   datas_descr,
#   metas,
#   designs,
#   spline_test_configs,
#   report_dir,
#   pthresholds,
#   )


## Run limma splines -----------------------------------------------------------

splineomics <- update_splineomics(
  splineomics = splineomics,
  design = "~ 1 + Phase*X + Reactor",
  data = data1,
  meta = meta1,
  spline_params = list(spline_type = c("n"),   # Chosen spline parameters
                       dof = c(2L))
)

# debug(run_limma_splines)
splineomics <- run_limma_splines(
  splineomics
)



# Extract the feature_nr column
feature_nr <- splineomics[["limma_splines_result"]][["time_effect"]][["Phase_Exponential"]][["feature_nr"]]

# Select the last 3000 entries
last_3000_feature_nr <- tail(feature_nr, 3000)

# Define the output file path for the RData file
output_rdata_file <- "inst/extdata/last_3000_feature_nr.RData"

# Save the list as RData
save(last_3000_feature_nr, file = output_rdata_file)




report_dir <- here::here("results", "limma_reports")

plots <- create_limma_report(
  splineomics,
  adj_pthresh = 0.1,
  report_dir = report_dir
)


## Cluster hits ----------------------------------------------------------------
adj_pthresholds <- c(0.05, 0.05)   
clusters <- c(6L, 3L)   
report_dir <- here::here("results", "clustering_reports")

plot_info = list(
  y_axis_label = "log2 intensity",
  time_unit = "min",
  treatment_labels = c("Feeding", "Test"),
  treatment_timepoints = c(10, 200)
)

# debug(cluster_hits)
clustering_results <- cluster_hits(
  splineomics = splineomics,
  analysis_type = "time_effect",
  adj_pthresholds = adj_pthresholds,
  clusters = clusters,
  genes = genes,
  plot_info = plot_info,
  report_dir = report_dir,
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

# download_enrichr_databases(gene_set_lib)
 
# # Get gene vector
# genes <- sub(" .*", "", genes)
# genes <- sub(";.*", "", genes)
# genes <- sub("_.*", "", genes)
# genes <- sub("-.*", "", genes)


downloaded_dbs_filepath <- 
  here::here(
    "dev",
    "data",
    "all_databases_08_04_2024-12_41_50.tsv"
  )

databases <- readr::read_tsv(downloaded_dbs_filepath, col_types = readr::cols())

clusterProfiler_params <- list(
  adj_p_value = 0.05,
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
