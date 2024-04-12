# libraries --------------
library(limma)
library(splines)
library(ggplot2)
library(tidyr)
library(dplyr)
library(tibble)
library(dendextend)
library(RColorBrewer)
library(patchwork)


# Define subfuntioncs --------






# Working section -----

rm(list = ls(all.names = TRUE))
# options(error = recover)

# Setup ------------------------------------------------------------------------

## Source functions ---------------------------------
limma_hyperparams_screen_fun_path <- 
  here::here("R", "limma_hyperparams_screen.R")
source(limma_hyperparams_screen_fun_path)

run_limma_splines_fun_path <- here::here("R", "run_limma_splines.R")
source(run_limma_splines_fun_path)

# cluster_hits_fun_path <- here::here("R", "cluster_hits.R")
# source(cluster_hits_fun_path)



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
metas <- list(meta1, meta2)
designs <- c("~ 1 + Phase*X + Reactor", "~ 1 + X + Reactor")
modes <- c("integrated", "isolated")
condition <- "Phase"
DoFs <- c(2L, 3L, 4L, 5L)
feature_names <- annotation$First.Protein.Description
pthresholds <- c(0.05, 0.1)

## hyperparams screen limma ----------------------------------------------------
top_tables_combos <- get_limma_combos_results(datas, 
                                              metas, 
                                              designs, 
                                              modes, 
                                              condition, 
                                              DoFs, 
                                              feature_names, 
                                              pthresholds, 
                                              padjust_method = "BH")

combo_pair_plots <- plot_limma_combos_results(top_tables_combos, 
                                              datas, 
                                              metas)

# gen_hyperparams_screen_reports(combo_pair_plots)

header_section <- '
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Page Title</title>
    <link rel="stylesheet" href="style.css">
    <meta name="description" content="A brief description of the page">
    <meta name="keywords" content="HTML, CSS, JavaScript">
    <meta name="author" content="Your Name">
    <!-- You can add more meta tags here -->
</head>
<body>

<!-- Page content goes here -->

</body>
</html>
'

convertPlotToBase64ImgTag <- function(plot, plot_nrows, width = 7, 
                                      base_height_per_row = 2.5, units = "in", 
                                      html_img_width="100%") {
  
  additional_height_per_row <- 2.1
  height <- base_height_per_row + (plot_nrows - 1) * 
    additional_height_per_row
  
  # Temporarily save the plot as an image
  img_file <- tempfile(fileext = ".png")
  ggsave(img_file, plot=plot, width = width, height = height, units = units, 
         limitsize = FALSE)
  
  # Convert the image to a Base64 string
  img_base64 <- base64enc::dataURI(file = img_file, mime = "image/png")
  
  # Delete the temporary image file
  unlink(img_file)
  
  # Return the HTML img tag with the Base64 string and a fixed width
  return(sprintf('<img src="%s" alt="Plot" style="width:%s;">', 
                 img_base64, html_img_width))
}


build_plot_report_html <- function(header_section, plots, rowcounts, path) {
  
  
  
  index <- 1
  for (plot in plots) {
    plot_nrows <- rowcounts[[index]]
    index <- index + 1
    img_tag <- convertPlotToBase64ImgTag(plot, plot_nrows)
    header_section <- paste(header_section, img_tag, sep="\n")
  }
  
  
  # Close the HTML document
  html_content <- paste(header_section, "</body></html>", sep="\n")
  
  
  dir_path <- dirname(path)
  
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  writeLines(html_content, path)
}


generate_report_html <- function(plot_list, plot_list_nrows, report_dir) {
  timestamp <- format(Sys.time(), "%d_%m_%Y-%H_%M_%S")
  
  omics_data_type <- "PTX"
  header_text <- paste("limma clustered hits", 
                       omics_data_type, timestamp, sep=" | ")
  
  if (omics_data_type == "PTX") {
    design_text <- "splines | X,data = exp AND stat | limma design: 1 + 
      Phase*X + Reactor | timepoints: E12_TP05_Exponential & 
      E10_TP10_Stationary removed<br>(Note: batch-corrected data used for 
      plotting
      individual blue datapoints (plots with the red spline, below))"
  } else if (omics_data_type == "PPTX") {
    design_text <- "splines | X,data = exp AND stat | limma design: 1 + 
      Phase*X + Reactor | timepoints: all<br>(Note: batch-corrected data used 
      for plotting
      individual blue datapoints (plots with the red spline, below))"
  }
  
  html_content <- paste(
    "<html><head><title>My Plots</title></head><body>",
    "<h1 style='color:red;'>", header_text, "</h1>",
    "<h2>Design</h2>",
    "<p>", design_text, "</p>",
    "</body></html>",
    sep=""
  )
  
  file_name <- sprintf("limma_hyperparams_screen_report_%s_%s.html",
                       omics_data_type, timestamp)
  
  output_file_path <- here::here(report_dir, file_name)
  
  build_plot_report_html(html_content, 
                         plot_list, 
                         plot_list_nrows, 
                         output_file_path)
}


process_combo_pair <- function(combo_pair) {
  
  plots <- list()
  plots_len <- integer(0)
  
  hitcomp <- combo_pair$hitcomp

  plots[[1]] <- hitcomp$vennheatmap
  plots[[2]] <- hitcomp$vennheatmap
  plots[[3]] <- hitcomp$barplot
  
  plots_len <- c(plots_len, 
                 2, 
                 hitcomp$vennheatmap_len, 
                 hitcomp$barplot_len)
  
  
  composites <- combo_pair$composites
  
  for (composite in composites) {
    # Append 'composite_plots' ggplot objects to 'plots'
    for (plot in composite$composite_plots) {
      plots[[length(plots) + 1]] <- plot
    }
    
    # Append 'composite_plots_len' values to 'plots_len'
    for (len in composite$composite_plots_len) {
      plots_len <- c(plots_len, len)
    }
  }
  
  generate_report_html(plots, plots_len, report_dir)
}

report_dir <- here::here("results", "hyperparams_screen_reports")

debugonce(process_combo_pair)
result <- map(combo_pair_plots, process_combo_pair)
