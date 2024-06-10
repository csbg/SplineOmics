rm(list = ls(all.names = TRUE))    # Just to ensure that


# Arguments of script -----------------------------------------------------

omics_datatype <- "PTX"    # Choose PTX or PPTX

clusterProfiler_params <- list(adj_p_value = 0.05, pAdjustMethod = "BH",
                               minGSSize = 10, maxGSSize = 500,
                               qvalueCutoff = 0.2)

use_background = TRUE



# Function definitions ----------------------------------------------------


create_log_file <- function(base_dir, data, omics_datatype, phase,
                            log_file_df_list) {
  # Generate timestamp and log filename
  timestamp <- format(Sys.time(), "%d_%m_%Y__%H_%M_%S")
  log_filename <- paste0("/gsea_log_", "_", omics_datatype, "_",
                         phase, "_", timestamp, ".csv")

  # Ensure the base directory exists; create it if it doesn't
  if (!dir.exists(base_dir)) {
    dir.create(base_dir, recursive = TRUE)
  }

  # Complete filepath for the log file
  log_filepath <- paste0(base_dir, log_filename)

  # Open a connection to the file for writing
  con <- file(log_filepath, open = "w")

  for (datafile_name in data) {
    writeLines(sprintf("Data: %s", datafile_name), con)
  }

  # Iterate over each item in the list and write to the file
  for (name in names(log_file_df_list)) {
    writeLines(name, con)  # Write the name as the header
    write.table(log_file_df_list[[name]], con, sep = ",", row.names = FALSE)
    writeLines(c("", ""), con)  # Write two blank lines
  }

  # Close the file connection
  close(con)
}



# Source fun and load data -----------------------------------------------------

## Source functions ------------------------------------------------------------
fun_source_path <- here::here("src", "enrichment", "clusterProfiler_fun.R")
source(fun_source_path)

## Load data -------------------------------------------------------------------
# State the RData filenames with the input data.
if (omics_datatype == "PTX") {
  clustered_genes_filenames <-
    list(paste("ID4_28_03_2024-17_03_32_adj_p_value_below_0.05_clustered_hit_",
               "genes_exp_PTX.RData", sep = ""),
         paste("ID4_28_03_2024-17_03_32_adj_p_value_below_0.05_clustered_hit_",
               "genes_stat_PTX.RData", sep = ""))
} else if (omics_datatype == "PPTX") {
  clustered_genes_filenames <-
    list(paste0("ID4_20_03_2024-15_33_57_adj_p_value_below_0.05_clustered_hit_",
                "genes_exp_PPTX.RData"),
         paste0("ID4_20_03_2024-15_33_57_adj_p_value_below_0.05_clustered_hit_",
                "genes_stat_PPTX.RData"))
}

loaded_df_names <- c()

for (filename in clustered_genes_filenames) {
  filepath <- here::here("results", "clust_hit_lists", "storage", filename)
  loaded_names <- load(filepath)
  loaded_dataframe_name <- loaded_names[1]
  loaded_df_names <- c(loaded_df_names, loaded_dataframe_name)
}

# Background (for clusterProfiler)
if (use_background) {
  background <- toupper(all_genes_names)
} else {
  background <- NULL
}



## Download databases (only do it once) ----------------------------------------

# Specified databases
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

download_databases(gene_set_lib)



# Enrichment results generation ------------------------------------------------

phases <- list("exp", "stat")

downloaded_dbs_filepath <- here::here("databases",
                                      "all_databases_08_04_2024-12_41_50.tsv")
all_term2genes <- dbs_to_term2genes(downloaded_dbs_filepath)

params_str <- sapply(names(clusterProfiler_params), function(name) {
  paste(name, clusterProfiler_params[[name]], sep = ": ")
})
params_flat <- paste(params_str, collapse = ", ")
elements <- strsplit(params_flat, ", ")[[1]]

params_string <- ""

for (i in 1:length(elements)) {
  params_string <- paste0(params_string, elements[i],
                          ifelse(i %% 2 == 0, "\n", ", "))
}

params_string <- sub(", $", "", params_string)
params_string <- sub("\n$", "", params_string)


for (phase in phases) {
  plot_title <- paste(omics_datatype, phase, "\n",
                      "(< 2 genes -> filtered out,\n any cluster only top 5)",
                      "\nParameters: ", params_string, sep = " ")

  clustered_genes_df_name <- paste0("clustered_hit_genes_", phase)
  clustered_genes <- get(clustered_genes_df_name)

  result <- pcgsea(clustered_genes, all_term2genes, clusterProfiler_params,
                   plot_title, background)

  print(result$dotplot)

  base_dir <- here::here("results", "clust_hit_lists", "log_file")
  create_log_file(base_dir, clustered_genes_filenames, omics_datatype,
                  phase, result$raw_results)
}
