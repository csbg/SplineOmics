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
  read_excel(
    system.file(
      "extdata",
      "proteomics_data.xlsx",
      package = "SplineOmics"
      )
    )
meta <- read_excel(
  system.file(
    "extdata",
    "proteomics_meta.xlsx",
    package = "SplineOmics"
    )
  )

# Extract the annotation part from the dataframe.
first_na_col <- which(is.na(data_excel[1,]))[1]
annotation <- data_excel |>
  dplyr::select((first_na_col + 1):ncol(data_excel)) |>
  dplyr::slice(-c(1:3))

print(data_excel)
print(annotation)
print(meta)

## ----process inputs, eval = TRUE----------------------------------------------
data <- extract_data(
  data = data_excel,
  feature_name_columns = c("Gene_name"),
  user_prompt = FALSE
  )

## ----Load EDA arguments, eval = TRUE------------------------------------------
report_info <- list(
  omics_data_type = "PTX",
  data_description = "Proteomics data of CHO cells",
  data_collection_date = "February 2024",
  analyst_name = "Thomas Rauter",
  contact_info = "thomas.rauter@plus.ac.at",
  project_name = "DGTX"
  )

report_dir <- here::here(
  "results",
  "explore_data"
  )

## ----Create the SplineOmics object, eval = TRUE-------------------------------
splineomics <- create_splineomics(
  data = data,
  meta = meta,
  annotation = annotation,
  report_info = report_info,
  condition = "Phase",
  meta_batch_column = "Reactor"
)

## ----Run EDA function, eval = FALSE-------------------------------------------
#  plots <- explore_data(
#    splineomics = splineomics,
#    report_dir = report_dir
#    )
#  

## ----Load hyperparameter-screening args, eval = TRUE--------------------------
data1 <- data 
meta1 <- meta

data2 <- data[, !(colnames(data) %in% c(
  "E12_TP05_Exponential", 
  "E10_TP10_Stationary"
  )
  )]
meta2 <- meta[!meta$`Sample.ID` %in% c(
  "E12_TP05_Exponential", 
  "E10_TP10_Stationary"
  ), ]

datas <- list(data1, data2) 
datas_descr <- c(
  "full_data",
  "outliers_removed"
  ) 

metas <- list(meta1, meta2) 
designs <- c(
  "~ 1 + Phase*X + Reactor",
  "~ 1 + X + Reactor"
  ) 
condition <- "Phase" 
report_dir <- here::here(
  "results",
  "hyperparams_screen_reports"
  ) 
meta_batch_column = "Reactor" 
pthresholds <- c(
  0.05,
  0.1
  )

# Every row a combo to test.
spline_test_configs <- data.frame(
  spline_type = c("n", "n", "n", "n"),
  degree = c(NA, NA, NA, NA),
  dof = c(2L, 3L, 4L, 5L),
  knots = I(list(c(NA), c(NA), c(NA), c(NA))),                                                        bknots = I(list(c(NA), c(NA), c(NA), c(NA)))
  )

## ----Perform hyperparameter-screening, eval = FALSE---------------------------
#  screen_limma_hyperparams(
#    splineomics,
#    datas,
#    datas_descr,
#    metas,
#    designs,
#    spline_test_configs,
#    report_dir,
#    pthresholds,
#    )
#  

## ----Update the SplineOmics object, eval = TRUE-------------------------------
splineomics <- update_splineomics(
  splineomics = splineomics,
  design = "~ 1 + Phase*X + Reactor",
  data = data1,  # data1 and meta1 == data and meta that were loaded before. 
  meta = meta1,  # This is for illustration
  spline_params = list(
    spline_type = c("n"),   
    dof = c(2L)
    )
)


## ----limma spline analysis, eval = TRUE---------------------------------------

# Run the limma spline analysis
splineomics <- run_limma_splines(
  splineomics
  )

## ----build limma report, eval = FALSE-----------------------------------------
#  report_dir <- here::here(
#    "results",
#    "create_limma_reports"
#    )
#  
#  plots <- create_limma_report(
#    splineomics,
#    report_dir = report_dir
#    )

## ----cluster the hits, eval = FALSE-------------------------------------------
#  adj_pthresholds <- c(
#    0.05,
#    0.05
#    )
#  
#  clusters <- c(
#    6L,
#    3L
#    )
#  
#  report_dir <- here::here(
#    "results",
#    "clustering_reports"
#    )
#  
#  plot_info = list(
#    y_axis_label = "log2 intensity",
#    time_unit = "min",
#    treatment_labels = c("Feeding"),
#    treatment_timepoints = c(0)
#  )
#  
#  gene_column_name <- "Gene_symbol"
#  genes <- data_excel[[gene_column_name]][4:nrow(data_excel)]
#  
#  clustering_results <- cluster_hits(
#    splineomics = splineomics,
#    analysis_type = "time_effect",
#    adj_pthresholds = adj_pthresholds,
#    clusters = clusters,
#    genes = genes,
#    plot_info = plot_info,
#    report_dir = report_dir,
#    )

## ----download Enrichr databases, eval = FALSE---------------------------------
#  gene_set_lib <- c(
#    "WikiPathways_2019_Human",
#    "NCI-Nature_2016",
#    "TRRUST_Transcription_Factors_2019",
#    "MSigDB_Hallmark_2020",
#    "GO_Cellular_Component_2018",
#    "CORUM",
#    "KEGG_2019_Human",
#    "TRANSFAC_and_JASPAR_PWMs",
#    "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
#    "GO_Biological_Process_2018",
#    "GO_Molecular_Function_2018",
#    "Human_Gene_Atlas"
#    )
#  
#  download_enrichr_databases(gene_set_lib)

## ----run GSEA, eval = FALSE---------------------------------------------------
#  # Modify gene vector
#  genes <- sub(" .*", "", genes)
#  genes <- sub(";.*", "", genes)
#  genes <- sub("_.*", "", genes)
#  genes <- sub("-.*", "", genes)
#  
#  downloaded_dbs_filepath <-
#    here::here("all_databases_08_04_2024-12_41_50.tsv")
#  databases <- readr::read_tsv(
#    downloaded_dbs_filepath,
#    col_types = readr::cols()
#    )
#  
#  clusterProfiler_params <- list(
#    adj_p_value = 0.05,
#    pAdjustMethod = "BH",
#    minGSSize = 10,
#    maxGSSize = 500,
#    qvalueCutoff = 0.2
#    )
#  
#  report_dir <- here::here(
#    "results",
#    "gsea_reports"
#    )
#  
#  result <- create_gsea_report(
#    levels_clustered_hits = clustering_results$clustered_hits_levels,
#    genes = genes,
#    databases = databases,
#    params = clusterProfiler_params,
#    report_info = report_info,
#    report_dir = report_dir
#    )

