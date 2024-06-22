rm(list = ls(all.names = TRUE))

# Setup ------------------------------------------------------------------------
library(devtools)
devtools::load_all()

library(readxl)
library(conflicted)

# interactive_demo()



# Load the data ----------------------------------------------------------------

# data_excel <- read_excel(here::here("dev" ,"data", "PPTX",
#                                     "PPTX_processed_with_imputation.xlsx"))
# 
# meta <- read_excel(here::here("dev" ,"data",
#                               "Time_course_PPTX_old_metadata.xlsx"))


data_excel <- read_excel(here::here("inst" ,"extdata", "data.xlsx"))

meta <- read_excel(here::here("inst" ,"extdata", "meta.xlsx"))


# Get the gene list ------------------------------------------------------------

data <- as.data.frame(data_excel)
gene_column_name <- "Genes"
genes <- data[[gene_column_name]][4:nrow(data)]




# Automatically extract data matrix from excel/csv table -----------------------

# feature_name_columns <- c("Unique identifier",
#                           "Protein", 
#                           "Positions within proteins")

feature_name_columns <- c("First.Protein.Description")

# debug(extract_data)
data <- extract_data(data_excel,
                     feature_name_columns)


# Remove outliers --------------------------------------------------------------

# data <- data %>%
#   select(-E12_TP05_Exponential, -E10_TP10_Stationary)
# 
# meta <- meta %>%
#   filter(!Sample.ID %in% c("E12_TP05_Exponential", "E10_TP10_Stationary"))

# Explore data -----------------------------------------------------------------

report_info <- list(
  omics_data_type = "PTX",
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
#                       meta_batch2_column = NA,
#                       report_dir = report_dir,
#                       report = TRUE)


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
# debug(screen_limma_hyperparams)
# screen_limma_hyperparams(datas,
#                          datas_descr,
#                          metas,
#                          designs,
#                          condition,
#                          spline_test_configs,
#                          report_info,
#                          report_dir,
#                          pthresholds,
#                          meta_batch_column)


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

# create_limma_report(run_limma_splines_result = result,
#                     report_info = report_info,
#                     report_dir = report_dir)



## Cluster hits ----------------------------------------------------------------
adj_pthresholds <- c(0.05, 0.05)   
clusters <- list(2L, 2L)   
report_dir <- here::here("results", "clustering_reports")

combo_list <- list(top_tables1, top_tables1)

plot_info <- list(y_axis_label = "log2 intensity",
                  time_unit = "min",
                  treatment_labels = "feeding",
                  treatment_timepoints = 0)

# debug(cluster_hits)
clustering_results <- cluster_hits(top_tables = top_tables1, 
                                   data = data, 
                                   meta = meta, 
                                   design = design,
                                   condition = condition, 
                                   spline_params = spline_params,
                                   adj_pthresholds = adj_pthresholds,
                                   clusters = clusters,
                                   report_info = report_info,
                                   genes = genes,
                                   meta_batch_column = meta_batch_column,
                                   # meta_batch2_column = meta_batch2_column,
                                   plot_info = plot_info,
                                   report_dir = report_dir,
                                   report = TRUE)



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
                  "Human_Gene_Atlas")

# download_enrichr_databases(gene_set_lib)

# Get gene vector
data <- as.data.frame(data_excel)
gene_column_name <- "Gene Names"
genes <- data[[gene_column_name]][4:nrow(data)]
genes <- sub(" .*", "", genes)
genes <- sub(";.*", "", genes)
genes <- sub("_.*", "", genes)


downloaded_dbs_filepath <- 
  here::here("data", "all_databases_08_04_2024-12_41_50.tsv")

databases <- readr::read_tsv(downloaded_dbs_filepath, col_types = readr::cols())

clusterProfiler_params <- list(adj_p_value = 0.05,
                               pAdjustMethod = "BH",
                               minGSSize = 10,
                               maxGSSize = 500,
                               qvalueCutoff = 0.2)

report_dir <- here::here("results", "gsea_reports")

# debug(run_gsea)
result <- create_gsea_report(levels_clustered_hits = 
                                clustering_results$clustered_hits_levels,
                              genes = genes,
                              databases = databases,
                              params = clusterProfiler_params,
                              report_info = report_info,
                              report_dir = report_dir)

