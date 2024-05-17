rm(list = ls(all.names = TRUE))
# options(error = recover)


# Setup ------------------------------------------------------------------------
# Import libraries

# limma_hyperparams_screen()
# # library(pheatmap)
# library(ggplot2)
# library(stats)
# library(purrr)
# library(furrr)
# library(tidyr)
# library(dplyr)
# library(patchwork)
# library(stringr)
# library(progress)
# library(here)
# library(knitr)
# library(kableExtra)
# 
# 
# # run_limma_splines()
# # library(limma)
# library(splines)
# library(purrr)
# 
# # cluster_hits()
# library(limma)
# library(splines)
# library(ggplot2)
# library(tidyr)
# library(dplyr)
# library(tibble)
# library(dendextend)
# library(RColorBrewer)
# library(patchwork)
# library(ComplexHeatmap)
# library(circlize)
# library(grid)
# library(cluster)
# 
# # utils
# library(ragg)

library(import)
import::from(ComplexHeatmap, Heatmap)
import::from(ComplexHeatmap, draw)
import::from(ComplexHeatmap, ht_opt)
import::from(RColorBrewer, brewer.pal)
import::from(RColorBrewer, brewer.pal.info)
import::from(base64enc, dataURI)
import::from(cluster, silhouette)
import::from(dendextend, as.ggdend)
import::from(dendextend, color_branches)
import::from(dplyr, all_of)
import::from(dplyr, arrange)
import::from(dplyr, bind_rows)
import::from(dplyr, filter)
import::from(dplyr, group_by)
import::from(dplyr, left_join)
import::from(dplyr, mutate)
import::from(dplyr, pull)
import::from(dplyr, relocate)
import::from(dplyr, select)
import::from(dplyr, summarize)
import::from(dplyr, vars)
import::from(grid, gpar)
import::from(ggplot2, aes)
import::from(ggplot2, annotate)
import::from(ggplot2, element_blank)
import::from(ggplot2, element_text)
import::from(ggplot2, expansion)
import::from(ggplot2, facet_wrap)
import::from(ggplot2, geom_col)
import::from(ggplot2, geom_line)
import::from(ggplot2, geom_point)
import::from(ggplot2, geom_text)
import::from(ggplot2, geom_vline)
import::from(ggplot2, ggplot)
import::from(ggplot2, ggsave)
import::from(ggplot2, ggtitle)
import::from(ggplot2, labs)
import::from(ggplot2, scale_color_brewer)
import::from(ggplot2, scale_colour_manual)
import::from(ggplot2, scale_x_continuous)
import::from(ggplot2, scale_y_continuous)
import::from(ggplot2, theme)
import::from(ggplot2, theme_bw)
import::from(ggplot2, theme_minimal)
import::from(ggplot2, xlab)
import::from(ggplot2, ylab)
import::from(ggplot2, unit)
import::from(grDevices, dev.off)
import::from(here, here)
import::from(htmltools, save_html)
import::from(kableExtra, kable_styling)
import::from(knitr, kable)
import::from(limma, eBayes)
import::from(limma, lmFit)
import::from(limma, topTable)
import::from(limma, removeBatchEffect)
import::from(magrittr, `%>%`)
import::from(patchwork, plot_annotation)
import::from(patchwork, wrap_plots)
import::from(pheatmap, pheatmap)
import::from(progress, progress_bar)
import::from(purrr, flatten_chr)
import::from(purrr, imap)
import::from(purrr, map)
import::from(purrr, map2)
import::from(purrr, map_chr)
import::from(purrr, map_int)
import::from(purrr, partial)
import::from(purrr, pmap)
import::from(purrr, set_names)
import::from(ragg, agg_png)
import::from(splines, bs)
import::from(splines, ns)
import::from(stats, as.dendrogram)
import::from(stats, as.formula)
import::from(stats, coef)
import::from(stats, cutree)
import::from(stats, dist)
import::from(stats, hclust)
import::from(stats, model.matrix)
import::from(stats, relevel)
import::from(stats, setNames)
import::from(stringr, str_extract)
import::from(tibble, column_to_rownames)
import::from(tibble, enframe)
import::from(tibble, rownames_to_column)
import::from(tidyr, as_tibble)
import::from(tidyr, expand_grid)
import::from(tidyr, pivot_longer)
import::from(tidyr, pivot_wider)
import::from(tidyr, replace_na)
import::from(tidyr, separate)
import::from(tidyr, unnest_longer)
import::from(tools, file_path_sans_ext)
import::from(utils, combn)
import::from(utils, tail)



## Source functions ---------------------------------
limma_hyperparams_screen_fun_path <-
  here::here("R", "limma_hyperparams_screen.R")
source(limma_hyperparams_screen_fun_path)

run_limma_splines_fun_path <- here::here("R", "run_limma_splines.R")
source(run_limma_splines_fun_path)

cluster_hits_fun_path <- here::here("R", "cluster_hits.R")
source(cluster_hits_fun_path)

general_fun_path <- here::here("R", "utils.R")
source(general_fun_path)



# Data loading and processing --------------------------------------------------

# Input to the whole package are the standardizes dataframes data, meta, and 
# annotation. data contains the raw data, meta the column descriptions of data
# (the sample descriptions), and annotation the row descriptions (the feature
# descriptions)
input_file_path <- here::here("data", "timeseries_proteomics_example.RData")
load(input_file_path) 

# Prep input to hyperparams screen function ------------------------------------
data <- as.matrix(data)
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
feature_names <- annotation$First.Protein.Description
report_dir <- here::here("results", "hyperparams_screen_reports")
meta_batch_column = "Reactor"
pthresholds <- c(0.05, 0.1)

# Every row a combo to test.
spline_test_configs <- data.frame(spline_type = c("n", "n", "n", "n"),
                                  degree = c(NA, NA, NA, NA),
                                  dof = c(2L, 3L, 4L, 5L),
                                  knots = I(list(c(NA), c(NA), c(NA), c(NA))),
                                  bknots = I(list(c(NA), c(NA), c(NA), c(NA))))

report_info <- list(
  omics_data_type = "PTX",
  data_description = "Gene expression levels over time.",
  data_collection_date = "2024-05-15",
  analyst_name = "Thomas Rauter",
  project_name = "DGTX",
  contact_info = "rauterthomas0@gmail.com"
)


# hyperparams screen limma -----------------------------------------------------
# debug(limma_hyperparams_screen)
# result <- limma_hyperparams_screen(datas,
#                                    datas_descr,
#                                    metas,
#                                    designs,
#                                    modes,
#                                    condition,
#                                    spline_test_configs,
#                                    feature_names,
#                                    report_info,
#                                    report_dir,
#                                    pthresholds,
#                                    meta_batch_column)


## Run limma splines -----------------------------------------------------------
DoFs <- c(2L, 2L)

design <- "~ 1 + Phase*X + Reactor"
# design <- "~ 1 + X + Reactor"


spline_params = list(spline_type = c("n"),
                     dof = c(2L))

# debug(run_limma_splines)
result <- run_limma_splines(data, 
                            meta, 
                            design, 
                            spline_params = spline_params,
                            condition,
                            feature_names, 
                            mode = "integrated")

top_tables <- result$top_tables
ttslc_factor_only <- result$ttslc_factor_only
ttslc_factor_time <- result$ttslc_factor_time


## Cluster hits ----------------------------------------------------------------
p_values <- c(0.05, 0.05)
clusters <- list(6L, 3L)
report_dir <- here::here("results", "clustering_reports")

report_info <- list(
  omics_data_type = "PTX",
  data_description = "Gene expression levels over time.",
  data_collection_date = "2024-05-15",
  analyst_name = "Thomas Rauter",
  project_name = "DGTX",
  dataset_name = "CHO Cell Study",
  limma_design = design,
  method_description = "Spline analysis for time-series data.",
  results_summary = "Identified significant changes in gene expression.",
  conclusions = "Potential biomarkers identified for further validation.",
  contact_info = "rauterthomas0@gmail.com"
)

# debug(cluster_hits)
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
