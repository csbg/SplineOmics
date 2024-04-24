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

general_fun_path <- here::here("R", "splinetime_general_fun.R")
source(general_fun_path)



# Data loading and processing --------------------------------------------------

# Input to the whole package are the standardizes dataframes data, meta, and 
# annotation. data contains the raw data, meta the column descriptions of data
# (the sample descriptions), and annotation the row descriptions (the feature
# descriptions)
input_file_path <- here::here("data", "PTX_input_data.RData")
load(input_file_path) 



# Run limma splines and cluster hits -------------------------------------------
data1 <- data
meta1 <- meta

data2 <- data
meta2 <- meta

datas <- list(data1, data2)
datas_descr <- c("full_data", "outliers_removed")
metas <- list(meta1, meta2)
designs <- c("~ 1 + Phase*X + Reactor", "~ 1 + X + Reactor")
modes <- c("integrated", "isolated")
condition <- "Phase"
DoFs <- c(2L, 3L, 4L, 5L)
feature_names <- annotation$First.Protein.Description
report_dir <- here::here("results", "jungfernflug2")
pthresholds <- c(0.05, 0.1)

# hyperparams screen limma -----------------------------------------------------
# debug(limma_hyperparams_screen)
result <- limma_hyperparams_screen(datas,
                                   datas_descr,
                                   metas,
                                   designs,
                                   modes,
                                   condition,
                                   DoFs,
                                   feature_names,
                                   report_dir,
                                   pthresholds)

## Run limma splines ----
DoFs <- c(2L, 2L)

design <- "~ 1 + Phase*X + Reactor"
# design <- "~ 1 + X + Reactor"

# debug(run_limma_splines)
result <- run_limma_splines(data, meta, design, DoFs, condition,
                                feature_names, "integrated")

top_tables <- result$top_tables
ttslc_factor_only <- result$ttslc_factor_only
ttslc_factor_time <- result$ttslc_factor_time

## Cluster hits ----
p_values <- c(0.05, 0.05)
clusters <- c(6L, 3L)
report_dir <- here::here("results")
data <- removeBatchEffect(x = data, batch = meta$Reactor)

clustering_results <- cluster_hits(top_tables, data, meta, condition, 
                                   p_values, clusters, report_dir)

clustering_results[[2]]$clustered_hits








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
