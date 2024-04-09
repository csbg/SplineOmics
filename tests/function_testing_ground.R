rm(list = ls(all.names = TRUE))
# options(error = recover)

# Setup ------------------------------------------------------------------------

## Source functions ---------------------------------
limma_hyperparams_screen_fun_path <- 
  here::here("R", "limma_hyperparams_screen.R")
source(limma_hyperparams_screen_fun_path)

run_limma_splines_fun_path <- here::here("R", "run_limma_splines.R")
source(run_limma_splines_fun_path)

cluster_hits_fun_path <- here::here("R", "cluster_hits.R")
source(cluster_hits_fun_path)



# Data loading and processing --------------------------------------------------

# Input to the whole package are the standardizes dataframes data, meta, and 
# annotation. data contains the raw data, meta the column descriptions of data
# (the sample descriptions), and annotation the row descriptions (the feature
# descriptions)
input_file_path <- here::here("data", "PTX_input_data.RData")
load(input_file_path) 



# Run limma splines and cluster hits --------------------------------------

## hyperparams screen limma ----
# limma_hyperparams_screen()

## Run limma splines ----
DoFs <- c(2L, 2L)

design <- "~ 1 + Phase*X + Reactor"
# design <- "~ 1 + X + Reactor"
group_factors <- c("Phase")
feature_names <- annotation$First.Protein.Description

# debug(run_limma_splines)
result <- run_limma_splines(data, meta, design, DoFs, group_factors,
                                feature_names, "integrated")

top_tables <- result$top_tables
ttslc <- result$ttslc

## Cluster hits ----
p_values <- c(0.05, 0.05)
clusters <- c(6L, 3L)
report_dir <- here::here("results")
data <- removeBatchEffect(x = data, batch = meta$Reactor)

clustering_results <- cluster_hits(top_tables, data, meta, group_factors, 
                                   p_values, clusters, report_dir)



clustering_results[[2]]$clustered_hits