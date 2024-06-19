## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(SplineOmics)
library(readxl)

## ----load the files-----------------------------------------------------------
data_excel <- 
  read_excel(system.file("extdata", "data.xlsx", package = "SplineOmics"))
meta <- read_excel(system.file("extdata", "meta.xlsx", package = "SplineOmics"))

print(data_excel)
print(meta)

## ----process inputs, eval = TRUE----------------------------------------------
data <- extract_data(data_excel, c("First.Protein.Description", "ID"))

## ----Load EDA arguments, eval = TRUE------------------------------------------
report_info <- list(
  omics_data_type = "PTX",
  data_description = "Proteomics data of CHO cells",
  data_collection_date = "February 2024",
  analyst_name = "Thomas Rauter",
  contact_info = "thomas.rauter@plus.ac.at",
  project_name = "DGTX")

condition <- "Phase"
meta_batch_column <- "Reactor"

report_dir <- here::here("results", "explore_data")

## ----Run EDA function, eval = FALSE-------------------------------------------
#  plots <- explore_data(data = data,
#                        meta = meta,
#                        condition = condition,
#                        report_info = report_info,
#                        meta_batch_column = meta_batch_column,
#                        report_dir = report_dir)
#  

## ----Load hyperparameter-screening args, eval = TRUE--------------------------
data1 <- data 
meta1 <- meta

data2 <- data[, !(colnames(data) %in% c("E12_TP05_Exponential", 
                                        "E10_TP10_Stationary"))]
meta2 <- meta[!meta$`Sample.ID` %in% c("E12_TP05_Exponential", 
                                       "E10_TP10_Stationary"), ]

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
                                  knots = I(list(c(NA), c(NA), c(NA), c(NA))),                                                        bknots = I(list(c(NA), c(NA), c(NA), c(NA))))

## ----Perform hyperparameter-screening, eval = FALSE---------------------------
#  screen_limma_hyperparams(datas,
#                           datas_descr,
#                           metas,
#                           designs,
#                           condition,
#                           spline_test_configs,
#                           report_info,
#                           report_dir,
#                           pthresholds,
#                           meta_batch_column)
#  

## ----limma spline analysis----------------------------------------------------

design <- "~ 1 + Phase*X + Reactor"     # Chosen limma design

spline_params = list(spline_type = c("n"),  # Chosen spline parameters
                     dof = c(2L))

# Run the limma spline analysis
result <- run_limma_splines(data2,   # Chosen version of the data
                            meta2,
                            design,
                            spline_params = spline_params,
                            condition)

top_tables1 <- result$time_effect 
print(top_tables1)

top_tables2 <- result$avrg_diff_conditions 
print(top_tables2)

top_tables3 <- result$interaction_condition_time 
print(top_tables3) 

## ----build limma report, eval = FALSE-----------------------------------------
#  report_dir <- here::here("results", "create_limma_reports")
#  
#  create_limma_report(run_limma_splines_result = result,
#                      report_info = report_info,
#                      report_dir = report_dir)

## ----cluster the hits---------------------------------------------------------
adj_pthresholds <- c(0.05, 0.05)   
clusters <- list(2L, 2L)   
report_dir <- here::here("results", "clustering_reports")

gene_column_name <- "Genes"
genes <- data_excel[[gene_column_name]][4:nrow(data_excel)]

clustering_results <- cluster_hits(top_tables = top_tables1, 
                                   data = data2, 
                                   meta = meta2, 
                                   design = design,
                                   condition = condition, 
                                   spline_params = spline_params,
                                   genes = genes,
                                   adj_pthresholds = adj_pthresholds,
                                   clusters = clusters,
                                   report_info = report_info,
                                   meta_batch_column = meta_batch_column,
                                   report_dir = report_dir,
                                   report = TRUE)

## ----download Enrichr databases, eval = FALSE---------------------------------
#  gene_set_lib <- c("WikiPathways_2019_Human",
#                    "NCI-Nature_2016",
#                    "TRRUST_Transcription_Factors_2019",
#                    "MSigDB_Hallmark_2020",
#                    "GO_Cellular_Component_2018",
#                    "CORUM",
#                    "KEGG_2019_Human",
#                    "TRANSFAC_and_JASPAR_PWMs",
#                    "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
#                    "GO_Biological_Process_2018",
#                    "GO_Molecular_Function_2018",
#                    "Human_Gene_Atlas")
#  
#  download_enrichr_databases(gene_set_lib)

## ----run GSEA, eval = FALSE---------------------------------------------------
#  # Get gene vector
#  gene_column_name <- "Genes"
#  genes <- annotation[[gene_column_name]][1:nrow(annotation)]
#  genes <- sub(" .*", "", genes)
#  genes <- sub(";.*", "", genes)
#  genes <- sub("_.*", "", genes)
#  genes <- sub("-.*", "", genes)
#  
#  downloaded_dbs_filepath <-
#    here::here("all_databases_08_04_2024-12_41_50.tsv")
#  databases <- readr::read_tsv(downloaded_dbs_filepath, col_types = readr::cols())
#  
#  clusterProfiler_params <- list(adj_p_value = 0.05,
#                                 pAdjustMethod = "BH",
#                                 minGSSize = 10,
#                                 maxGSSize = 500,
#                                 qvalueCutoff = 0.2)
#  
#  report_dir <- here::here("results", "gsea_reports")
#  
#  result <- create_gsea_report(levels_clustered_hits =
#                                    clustering_results$clustered_hits_levels,
#                               genes = genes,
#                               databases = databases,
#                               params = clusterProfiler_params,
#                               report_info = report_info,
#                               report_dir = report_dir)

