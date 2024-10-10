## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, eval = TRUE-------------------------------------------------------
library(SplineOmics)
library(readxl) # for loading Excel files
library(here) # For managing filepaths
library(dplyr) # For data manipulation

## ----load the files-----------------------------------------------------------
data <- readRDS(system.file(
  "extdata",
  "proteomics_data.rds",
  package = "SplineOmics"
))


meta <- read_excel(
  system.file(
    "extdata",
    "proteomics_meta.xlsx",
    package = "SplineOmics"
  )
)

# Extract the annotation part from the dataframe.
first_na_col <- which(is.na(data[1, ]))[1]
annotation <- data |>
  dplyr::select((first_na_col + 1):ncol(data)) |>
  dplyr::slice(-c(1:3))

print(head(data))
print(meta)
print(annotation)

## ----process inputs, eval = TRUE----------------------------------------------
data <- SplineOmics::extract_data(
  # The dataframe with the numbers on the left and info on the right.
  data = data,
  # Use this annotation column for the feature names.
  feature_name_columns = c("Gene_name"),
  # When TRUE, you must confirm that data is in the required format.
  user_prompt = FALSE
)

## ----Load EDA arguments, eval = TRUE------------------------------------------
# Those fields are mandatory, because we believe that when such a report is
# opened after half a year, those infos can be very helpful.
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
# splineomics now contains the SplineOmics object.
splineomics <- SplineOmics::create_splineomics(
  data = data,
  meta = meta,
  annotation = annotation,
  report_info = report_info,
  condition = "Phase", # Column of meta that contains the levels.
  meta_batch_column = "Reactor" # For batch effect removal
)

# Special print.SplineOmics function leads to selective printing
print(splineomics)

## ----Run EDA function, eval = FALSE-------------------------------------------
#  plots <- SplineOmics::explore_data(
#    splineomics = splineomics, # SplineOmics object
#    report_dir = report_dir
#  )

## ----Load hyperparameter-screening args, eval = TRUE--------------------------
data1 <- data
meta1 <- meta

# Remove the "outliers"
data2 <- data[, !(colnames(data) %in% c(
  "E12_TP05_Exponential",
  "E10_TP10_Stationary"
)
)]

# Adjust meta so that it matches data2
meta2 <- meta[!meta$Sample.ID %in% c(
  "E12_TP05_Exponential",
  "E10_TP10_Stationary"
), ]

# As mentioned above, all the values of one hyperparameter are stored
# and provided as a list.
datas <- list(data1, data2)

# This will be used to describe the versions of the data.
datas_descr <- c(
  "full_data",
  "outliers_removed"
)

metas <- list(meta1, meta2)

# Test two different limma designs
designs <- c(
  "~ 1 + Phase*X + Reactor",
  "~ 1 + X + Reactor"
)

# 'Integrated means' limma will use the full dataset to generate the results for
# each condition. 'Isolated' means limma will use only the respective part of
# the dataset for each condition. Designs that contain the condition column
# (here Phase) must have mode 'integrated', because the full data is needed to
# include the different conditions into the design formula.
modes <- c(
  "integrated",
  "isolated"
)

# Specify the meta "level" column
condition <- "Phase"

report_dir <- here::here(
  "results",
  "hyperparams_screen_reports"
)

# To remove the batch effect
meta_batch_column <- "Reactor"

# Test out two different p-value thresholds (inner hyperparameter)
adj_pthresholds <- c(
  0.05,
  0.1
)

# Create a dataframe with combinations of spline parameters to test
# (every row a combo to test)
spline_test_configs <- data.frame(
  # 'n' stands for natural cubic splines, b for B-splines.
  spline_type = c("n", "n", "b", "b"),
  # Degree is not applicable (NA) for natural splines.
  degree = c(NA, NA, 2L, 4L),
  # Degrees of freedom (DoF) to test.
  # Higher dof means spline can fit more complex patterns.
  dof = c(2L, 3L, 3L, 4L)
)

print(spline_test_configs)

## ----Perform hyperparameter-screening, eval = FALSE---------------------------
#  SplineOmics::screen_limma_hyperparams(
#    splineomics = splineomics,
#    datas = datas,
#    datas_descr = datas_descr,
#    metas = metas,
#    designs = designs,
#    modes = modes,
#    spline_test_configs = spline_test_configs,
#    report_dir = report_dir,
#    adj_pthresholds = adj_pthresholds,
#  )

## ----Update the SplineOmics object, eval = TRUE-------------------------------
splineomics <- SplineOmics::update_splineomics(
  splineomics = splineomics,
  design = "~ 1 + Phase*X + Reactor", # best design formula
  mode = "integrated", # means limma uses the full data for each condition.
  data = data2, # data without "outliers" was better
  meta = meta2,
  spline_params = list(
    spline_type = c("n"), # natural cubic splines (take these if unsure)
    dof = c(2L) # If you are unsure about which dof, start with 2 and increase
  )
)

## ----limma spline analysis, eval = TRUE---------------------------------------
splineomics <- SplineOmics::run_limma_splines(
  splineomics = splineomics
)

## ----build limma report, eval = FALSE-----------------------------------------
#  report_dir <- here::here(
#    "results",
#    "create_limma_reports"
#  )
#  
#  plots <- SplineOmics::create_limma_report(
#    splineomics = splineomics,
#    report_dir = report_dir
#  )

## ----cluster the hits, eval = FALSE-------------------------------------------
#  adj_pthresholds <- c( # 0.05 for both levels
#    0.05, # exponential
#    0.05 # stationary
#  )
#  
#  clusters <- c(
#    6L, # 6 clusters for the exponential phase level
#    3L # 3 clusters for the stationary phase level
#  )
#  
#  report_dir <- here::here(
#    "results",
#    "clustering_reports"
#  )
#  
#  plot_info <- list( # For the spline plots
#    y_axis_label = "log2 intensity",
#    time_unit = "min", # our measurements were in minutes
#    treatment_labels = list("feeding"), # add this for all conditions
#    treatment_timepoints = list(0) # Feeding was at 0 minutes.
#  )
#  
#  # Like this you can add individual treatment labels to your plots:
#  # treatment_labels = list(
#  #   exponential = "treatment 1", # One treatment in exp
#  #   stationary = c("treatment 2", "treatment 3")  # Two treatments in stat
#  #   additional_condition = NA  # No treatment in the hypothetical third condition
#  #   )
#  #
#  # treatment_timepoints = list(
#  #   exponential = 0,
#  #   stationary = C(100, 140),  # Two treatments also need two timepoints
#  #   additional_condition = NA
#  #   )
#  #
#  # or set a treatment for ALL conditions (still always make a list):
#  #
#  # treatment_labels = list("treatment")
#  # treatment_timepoints = list(120)
#  #
#  # or set multiple treatments for ALL conditions:
#  #
#  # treatment_labels = list(c("treatment1", "treatment2"))
#  # treatment_timepoints = list(c(120, 90))
#  
#  
#  # Get all the gene names. They are used for generating files
#  # which contents can be directly used as the input for the Enrichr webtool,
#  # if you prefer to manually perform the enrichment. Those files are
#  # embedded in the output HTML report and can be downloaded from there.
#  gene_column_name <- "Gene_symbol"
#  genes <- annotation[[gene_column_name]]
#  
#  plot_options <- list(
#    # When meta_replicate_column is not there, all datapoints are blue.
#    meta_replicate_column = "Reactor", # Colors the data points based on Reactor
#    cluster_heatmap_columns = FALSE # Per default FALSE, just for demonstration
#  )
#  
#  clustering_results <- SplineOmics::cluster_hits(
#    splineomics = splineomics,
#    adj_pthresholds = adj_pthresholds,
#    clusters = clusters,
#    genes = genes,
#    plot_info = plot_info,
#    plot_options = plot_options,
#    report_dir = report_dir,
#    adj_pthresh_avrg_diff_conditions = 0,
#    adj_pthresh_interaction_condition_time = 0.25
#  )

## ----download Enrichr databases, eval = FALSE---------------------------------
#  # Specify which databases you want to download from Enrichr
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
#  )
#  
#  SplineOmics::download_enrichr_databases(
#    gene_set_lib = gene_set_lib,
#    filename = "databases.tsv"
#  )

## ----prepare arguments for GSEA, eval = FALSE---------------------------------
#  # Specify the filepath of the TSV file with the database info
#  downloaded_dbs_filepath <- here::here("databases.tsv")
#  
#  # Load the file
#  databases <- read.delim(
#    downloaded_dbs_filepath,
#    sep = "\t",
#    stringsAsFactors = FALSE
#  )
#  
#  # Specify the clusterProfiler parameters
#  clusterProfiler_params <- list(
#    pvalueCutoff = 0.05,
#    pAdjustMethod = "BH",
#    minGSSize = 10,
#    maxGSSize = 500,
#    qvalueCutoff = 0.2
#  )
#  
#  report_dir <- here::here(
#    "results",
#    "gsea_reports"
#  )

## ----run GSEA, eval = FALSE---------------------------------------------------
#  result <- SplineOmics::run_gsea(
#    # A dataframe with three columns: feature, cluster, and gene. Feature contains
#    # the integer index of the feature, cluster the integer specifying the cluster
#    # number, and gene the string of the gene, such as "CLSTN2".
#    levels_clustered_hits = clustering_results$clustered_hits_levels,
#    databases = databases,
#    clusterProfiler_params = clusterProfiler_params,
#    report_info = report_info,
#    report_dir = report_dir
#  )

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

