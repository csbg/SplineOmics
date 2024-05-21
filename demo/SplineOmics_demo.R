library(SplineOmics)

# Proteomics example
input_file_path <- here::here("data", "timeseries_proteomics_example.RData")

# Phosphoproteomics example
# input_file_path <- 
#   here::here("data", "timeseries_phosphoproteomics_example.RData")

load(input_file_path) 



# hyperparams screen limma -----------------------------------------------------
# data must contain only the numeric values of the timeseries omics experiment.
# No other things, such as columns with feature descriptions are allowed.
data <- as.matrix(data)     
data1 <- data
meta1 <- meta

data2 <- data   # data2 is the same as data1, but just to show the functionality
meta2 <- meta

datas <- list(data1, data2)
datas_descr <- c("full_data", "outliers_removed")
metas <- list(meta1, meta2)

# limma design formulas
designs <- c("~ 1 + Phase*X + Reactor", "~ 1 + X + Reactor")

# Mode for the limma design formula 1 and 2, respectively. Integrated means
# limma identifies the hits of one level of the factor by taking into account
# the information from all the levels, and isolated means that a completely 
# separate limma run is carried out for each level (no interaction between the
# levels)
modes <- c("integrated", "isolated")

condition <- "Phase"                   # This is the factor of the experiment
feature_names <- annotation$First.Protein.Description

# The HTML reports will go in this dir (will be created automatically if it does
# not exist yet). Can be ommited, the default argument is the cwd.
report_dir <- here::here("results", "hyperparams_screen_reports")

# By specifying this, the pipeline will automatically remove the batch-effect.
# (Removing batch effect does not mean it is gone entirely, it means that the
# function from the R package limma, that was designed to tackle the batch
# effect, is run on the data and this batch effect removed data is then used
# for plotting, but not for the limma analysis.)
meta_batch_column = "Reactor"   
pthresholds <- c(0.05, 0.1)

# Every row a combo of hyperparameters to test.
spline_test_configs <- data.frame(spline_type = c("n", "n", "n", "n"),
                                  degree = c(NA, NA, NA, NA),
                                  dof = c(2L, 3L, 4L, 5L),
                                  knots = I(list(c(NA), c(NA), c(NA), c(NA))),
                                  bknots = I(list(c(NA), c(NA), c(NA), c(NA))))

# This info will be written on top of every report. It is mandatory, to enforce
# what in our opinion is a good organisation. Write in NA if you do not want
# to provide the info.
report_info <- list(
  omics_data_type = "PTX",
  data_description = "Timeseries Proteomics",
  data_collection_date = "2024-05-15",
  analyst_name = "Thomas Rauter",
  project_name = "DGTX",
  contact_info = "thomas.rauter@plus.ac.at"
)

# Test all the different hyperparameters and their combinations and generate
# the HTML reports with the results.
result <- limma_hyperparams_screen(datas,
                                   datas_descr,
                                   metas,
                                   designs,
                                   modes,
                                   condition,
                                   spline_test_configs,
                                   feature_names,
                                   report_info,
                                   report_dir,
                                   pthresholds,
                                   meta_batch_column)

# Now, in the workflow you would check out the generated HTMLs and investigate
# the impact of your hyperparameter choices. For example, you could have a look
# which degree of freedom (DoF) is best for your data (choose a DoF
# where the curve is not too wiggly (DoF too high), but is also not ignoring the 
# patterns in your data (DoF too low))

## Run limma splines -----------------------------------------------------------
design <- "~ 1 + Phase*X + Reactor"          # Chosen limma design
spline_params = list(spline_type = c("n"),   # Chosen spline parameters
                     dof = c(2L))

# Run the limma spline analysis
result <- run_limma_splines(data, 
                            meta, 
                            design, 
                            spline_params = spline_params,
                            condition,
                            feature_names, 
                            mode = "integrated")

# Contain the adj. p-value for all features (for the tested hypothesis if the
# features showed a change over the time or not)
top_tables <- result$top_tables  

# ttslc = top_tables_level_comparison (which features are significantly changed
# between the levels). These tables contain the adj. p-values for the tested
# hypothesis, if the features are the same for the levels or not.
ttslc_factor_only <- result$ttslc_factor_only    
ttslc_factor_time <- result$ttslc_factor_time    # Factor AND time


## Cluster hits ----------------------------------------------------------------
p_values <- c(0.05, 0.05)    # Chosen p-value threshold

# Try out cluster nr. until the clusters are "pure". Pure means the cluster 
# does not contain very distinctive shapes. Unfortunately, this is a little bit 
# of an art. Try to hit the sweet spot between cluster "purity" and amount of
# clusters (since more clusters means more complex subsequent analysis, you do
# not want to have very many clusters with only little hits in each cluster.)
clusters <- list(6L, 3L)   

# The resulting HTML reports will be placed here (path gets created if it does
# not exist yet)
report_dir <- here::here("results", "clustering_reports")

# This info will be written on top of the report. The fields until dataset_name 
# (included) are mandatory (the rest down there not), to enforce
# what in our opinion is a good organisation. Write in NA if you do not want
# to provide the info.
report_info <- list(
  omics_data_type = "PTX",
  data_description = "Timeseries Proteomics",
  data_collection_date = "2024-05-15",
  analyst_name = "Thomas Rauter",
  project_name = "DGTX",
  dataset_name = "CHO Cell Study",
  limma_design = design,
  method_description = "Spline analysis for time-series data.",
  results_summary = "Identified significant changes in the proteins.",
  conclusions = "Found the hits",
  contact_info = "thomas.rauter@plus.ac.at"
)

# Cluster the results from the limma spline analysis.
clustering_results <- cluster_hits(top_tables = top_tables, 
                                   data = data, 
                                   meta = meta, 
                                   condition = condition, 
                                   spline_params = spline_params,
                                   mode = "integrated",
                                   p_values = p_values,
                                   clusters = clusters,
                                   report_info = report_info,
                                   meta_batch_column = meta_batch_column,
                                   report_dir = report_dir
)

# This is the result in a df. This can be for example used as the input for
# automated gene set enrichment analysis with EnrichR or clusterProfiler.
clustering_results[[2]]$clustered_hits
