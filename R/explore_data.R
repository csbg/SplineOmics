# This function of the SplineOmics package automates the process of exploratory
# data analysis. Mandatory inputs are only data, meta, condition, and a list
# with report infos, that will be placed on top of the HTML reports with the 
# exploratory figures.



# Exported function: explore_data () -------------------------------------------


#' Generate Exploratory Plots
#'
#' @description
#' This function takes a data matrix, checks its validity, and generates a list 
#' of exploratory plots including density plots, boxplots, PCA plots, MDS plots, 
#' variance explained plots, and violin plots.
#'
#' @param data A dataframe containing the data.
#' @param meta A dataframe containing metadata corresponding to the `data`.
#' @param condition A character string specifying the condition (column of meta)
#' @param report_info A named list containing report information.
#' @param meta_batch_column A character string specifying the meta batch column.
#' @param meta_batch2_column A character string specifying the meta batch2
#'                           column.
#' @param report_dir A non-empty string specifying the report directory.
#'
#' @return A list of ggplot objects representing various exploratory plots.
#'
#' @export
#' 
explore_data <- function(data,
                         meta,
                         condition,
                         report_info,
                         meta_batch_column = NA,
                         meta_batch2_column = NA,
                         report_dir = here::here(),
                         report = TRUE) {
  
  # Control the function arguments
  args <- lapply(as.list(match.call()[-1]), eval, parent.frame())
  check_null_elements(args)
  input_control <- InputControl$new(args)
  input_control$auto_validate()

  data_report <- rownames_to_column(data, var = "Feature_name")
  
  matrix_and_feature_names <- process_data(data)
  data <- matrix_and_feature_names$data

  
  data_list <- list(data = data)
  
  if (!is.na(meta_batch_column)) {
    
    args <- list(
      x = data,
      batch = meta[[meta_batch_column]],
      group = meta[[condition]]
    )

    if (!is.na(meta_batch2_column)) {
      args$batch2 <- meta[[meta_batch2_column]]
    }

    batch_corrected_data <- do.call(removeBatchEffect, args)
    
    data_list$batch_corrected_data <- batch_corrected_data
  }
  
  all_plots <- list()
  report_info$meta_condition <- c(condition)
  report_info$meta_batch <- paste(meta_batch_column, 
                                  meta_batch2_column,
                                  sep = ", ")
  timestamp = format(Sys.time(), "%d_%m_%Y-%H_%M_%S")
  
  for (data_name in names(data_list)) {
    
    current_data <- data_list[[data_name]]
    plots_and_plots_sizes <- generate_explore_plots(current_data, 
                                                    meta, 
                                                    condition)
    
    if (report) {
      generate_report_html(plots = plots_and_plots_sizes$plots,
                           plots_sizes = plots_and_plots_sizes$plots_sizes,
                           report_info = report_info,
                           data = data_report,
                           meta = meta,
                           filename = paste0("explore_", data_name),
                           timestamp = timestamp,
                           report_dir = report_dir)
    }
    
    all_plots[[data_name]] <- plots_and_plots_sizes$plots
  }
  
  return(all_plots)
}


# Level 1 internal functions ---------------------------------------------------


#' #' Control Input Parameters for Data Exploration Functions
#' #'
#' #' @description
#' #' This function validates the input parameters for data exploration functions.
#' #' It checks the structure and data types of `data`, `meta`, `condition`,
#' #' `report_info`, `meta_batch_column`, `meta_batch2_column`, `report_dir`, 
#' #' and `report`. If any condition is not met, the function stops the script 
#' #' and produces an informative error message.
#' #'
#' #' @param data A dataframe containing the data to be explored.
#' #' @param meta A dataframe containing metadata associated with the data.
#' #' @param condition A character string specifying the condition to be explored.
#' #' @param report_info A list containing information for generating the report.
#' #' @param meta_batch_column A character string specifying the batch column in 
#' #'        the metadata.
#' #' @param meta_batch2_column A character string specifying the second batch 
#' #'        column in the metadata.
#' #' @param report_dir A character string specifying the directory where the 
#' #'        report will be saved.
#' #' @param report A Boolean indicating whether to generate a report (TRUE) or 
#' #'        not (FALSE).
#' #'
#' #' @return This function does not return a value. It stops with an error message 
#' #'         if the conditions are not met.
#' #'
#' control_inputs_explore_data <- function(data,
#'                                         meta,
#'                                         condition,
#'                                         report_info,
#'                                         meta_batch_column,
#'                                         meta_batch2_column,
#'                                         report_dir,
#'                                         report) {
#'   
#'   check_data_and_meta(data, 
#'                       meta, 
#'                       condition, 
#'                       meta_batch_column,
#'                       meta_batch2_column)
#'   
#'   check_report_info(report_info)
#'   
#'   check_and_create_report_dir(report_dir)
#'   
#'   
#'   if (report != TRUE && report != FALSE) {
#'     stop("report must be either Boolean TRUE or FALSE", call. = FALSE)
#'   }
#' }


#' Generate exploratory plots
#'
#' @description
#' This function generates various exploratory plots including density plots, 
#' box plots, violin plots, PCA plots, and correlation heatmaps based on the 
#' provided data and metadata.
#'
#' @param data A data frame or matrix containing the data to be plotted.
#' @param meta A data frame containing metadata associated with the data.
#' @param condition A string specifying the column in the metadata that contains 
#' the condition or grouping variable.
#'
#' @return A list containing two elements:
#' \describe{
#'   \item{plots}{A list of ggplot objects representing the generated plots.}
#'   \item{plots_sizes}{A vector of numeric values indicating the sizes of the 
#'   corresponding plots.}
#' }
#' 
generate_explore_plots <- function(data,
                                   meta,
                                   condition) {
  
  meta[[condition]] <- as.factor(meta[[condition]])
  
  plot_functions_and_sizes <- list(
    list(func = make_density_plots, size = 1),
    list(func = make_box_plots, size = 1.5),
    list(func = make_violin_plots, size = 1.5),
    list(func = plot_mean_correlation_with_time, size = 1.5),
    list(func = plot_lag1_differences, size = 1.5),
    list(func = plot_first_lag_autocorrelation, size = 1.5),
    list(func = plot_cv, size = 1.5),
    list(func = make_pca_plot, size = 1.5, flatten = FALSE),
    list(func = make_pca_variance_explained_plot, size = 1.5, flatten = FALSE),
    list(func = make_mds_plot, size = 1.5, flatten = FALSE),
    list(func = make_correlation_heatmaps, size = NULL)
  )
  
  apply_plot_function <- function(entry) {
    plot_result <- entry$func(data, meta, condition)
    list(plots = plot_result, size = entry$size, 
         flatten = if ("flatten" %in% names(entry)) entry$flatten else TRUE)
  }
  
  plot_results <- lapply(plot_functions_and_sizes, apply_plot_function)

  all_plots <- list()
  all_plots_sizes <- c()
  
  # Flatten the results and sizes conditionally
  for (result in plot_results) {
    if (is.null(result$size)) {
      # Special handling for make_correlation_heatmaps
      all_plots <- c(all_plots, result$plots$heatmaps)
      all_plots_sizes <- c(all_plots_sizes, result$plots$heatmaps_sizes)
    } else if (!is.null(result$flatten) && result$flatten == FALSE) {
      # Do not flatten the result, add it as a sublist
      all_plots <- c(all_plots, list(result$plots))
      all_plots_sizes <- c(all_plots_sizes, result$size)
    } else {
      # Flatten the result
      all_plots <- c(all_plots, result$plots)
      all_plots_sizes <- c(all_plots_sizes, rep(result$size, length(result$plots)))
    }
  }
  
  list(plots = all_plots, 
       plots_sizes = all_plots_sizes)
}


#' Build Explore Data Report
#'
#' @description
#' This function generates an HTML report containing a header section, table of 
#' contents, and a series of plots. Each plot is included in the report with 
#' specified sizes.
#'
#' @param header_section A string containing the HTML content for the header 
#' section of the report.
#' @param plots A list of ggplot objects representing the plots to be included 
#' in the report.
#' @param plots_sizes A list of sizes corresponding to each plot, defining the 
#' dimensions to be used when rendering the plots.
#' @param output_file_path A string specifying the file path where the HTML 
#' report will be saved.
#'
#' @return None. This function writes the HTML content to the specified file.
#' 
build_explore_data_report <- function(header_section, 
                                      plots, 
                                      plots_sizes, 
                                      output_file_path) {  
  
  level_count = (length(plots)-5)/8

  content_with_plots <- paste(header_section, "<!--TOC-->", sep="\n")
  
  toc <- "<div id='toc' style='text-align: center; display: block; margin: auto;
          width: 80%;'> 
        <h2 style='font-size: 40px;'>Table of Contents</h2>
        <ul style='display: inline-block; text-align: left;'>"
  
  section_header_style <- "font-size: 70px; color: #001F3F; text-align: center;"
  toc_style <- "font-size: 30px;"

  pb <- create_progress_bar(plots)
  plot_names <- c("Density Plots", 
                  "Boxplots", 
                  "Violin Plots",
                  "Mean Time Correlation",
                  "Lag1 Differences",
                  "First Lag Autocorrelation",
                  "Coefficient of Variation",
                  "PCA ", 
                  "PCA Variance explained",
                  "MDS",
                  "Correlation Heatmaps")
  
  plot_explanations <- c(
    "Density plots show the distribution of intensities or abundances across 
  samples. Peaks in the plot indicate common values, while the spread 
  indicates variability. Use this plot to identify the range and most 
  frequent values in your data.",
    
    "Boxplots display the spread and outliers of values for each sample. The 
  box represents the interquartile range (IQR), the line inside the box is 
  the median, and the whiskers extend to 1.5 * IQR from the quartiles. 
  Points outside this range are considered outliers. Use boxplots to 
  compare distributions and identify outliers.",
    
    "Violin plots combine boxplots and density plots to show the distribution 
  of values. They provide a summary of the data's range, central tendency, 
  and distribution shape. Use violin plots to understand the full 
  distribution and compare between groups.",
    
    "Mean Time Correlation plots summarize the correlation of each feature 
  with time, highlighting time-dependent trends. Positive correlations 
  indicate increasing trends, negative correlations indicate decreasing 
  trends. Use this plot to identify features that change over time.",
  
  "Purpose: Lag-1 Differences plots illustrate the changes in feature values 
  between successive time points. What It Shows: This plot helps identify 
  the variability and consistency in the changes of feature values over time.
  The mean lag-1 difference indicates the average change between time points.
  The standard deviation of the lag-1 differences indicates the variability 
  in these changes. How to Use It: Use this plot to understand the magnitude 
  and variability of changes in feature values over time. Identify features 
  with consistent increases or decreases and those with high variability 
  between time points.",  
  
    "First Lag Autocorrelation plots illustrate the temporal dependencies 
  within each feature by calculating the autocorrelation of the feature 
  with itself at a lag of one time point. This plot helps identify 
  features that have consistent patterns over time. A high positive 
  autocorrelation indicates that the feature values are similar to their 
  previous values, suggesting a trend or periodicity. A high negative 
  autocorrelation suggests an alternating pattern over time. Use this 
  plot to understand the persistence and cyclic behavior of the features 
  across time points, and to identify features that exhibit stable or 
  periodic patterns.",
  
  "Purpose: Coefficient of Variation (CV) plots illustrate the relative 
  variability of each feature by calculating the CV for each feature. 
  What It Shows: This plot helps identify features with high or low 
  relative variability. The mean CV indicates the average relative 
  variability across features. The standard deviation of CV indicates the 
  variability in the relative variability across features. How to Use It:
  Use this plot to understand the consistency and variability of feature 
  values. Identify features with consistently high or low variability 
  relative to their mean.",
    
  "PCA plots visualize the major trends and patterns in high-dimensional 
  data by reducing it to a few principal components. Points close to each 
  other are similar. Use PCA plots to identify clustering and variance 
  explained by the principal components.",
  
  "PCA Variance Explained plots show the amount of variance explained by 
  each principal component. The y-axis represents the percentage of total 
  variance explained, and the x-axis represents the principal components. 
  Use this plot to determine the number of components to consider.",
  
  "MDS plots display similarities or dissimilarities between samples in a 
  reduced dimension space. Points close to each other are more similar. 
  Use MDS plots to visualize the distance or similarity between samples.",
  
  "Correlation Heatmaps illustrate the pairwise correlation between all 
  samples. Colors represent the strength of correlation, with a color 
  gradient indicating positive or negative correlations. Use this plot to 
  identify highly correlated samples or groups."
  )

  toc_index <- 0
  toc_index_memory <- toc_index

  # Generate the sections and plots
  for (index in seq_along(plots)) {
    
    # Determine when to place headers based on the provided logic
    if (length(plots) == 7) {     # Just one level exists
      
      toc_index <- toc_index + 1
      
    } else if (index == 1 || 
               index == 2 + level_count || 
               index == 2 + 2 * level_count || 
               index == 2 + 3 * level_count ||
               index == 2 + 4 * level_count ||
               index == 2 + 5 * level_count ||
               index == 2 + 6 * level_count ||
               index == 2 + 7 * level_count || 
               index == 3 + 7 * level_count || 
               index == 4 + 7 * level_count ||
               index == 5 + 7 * level_count) {    # More than just one level
    
      toc_index <- toc_index + 1
    }
    
    
    if (toc_index != toc_index_memory) {
      
      section_id <- paste0("section_", toc_index)
      toc <- paste(toc, 
                   sprintf('<li style="%s"><a href="#%s">%s</a></li>', 
                           toc_style, 
                           section_id,
                           plot_names[toc_index]), 
                   sep = "\n")
      
      section_header <- sprintf('<h2 id="%s" style="%s">%s</h2>', 
                                section_id, 
                                section_header_style, 
                                plot_names[toc_index])
      
      plot_description <- sprintf('<p style="font-size: 2em;">%s</p>',
                                  plot_explanations[toc_index])
      
      
      content_with_plots <- paste(content_with_plots, 
                                  section_header, 
                                  plot_description, 
                                  sep = "\n")
      
      toc_index_memory <- toc_index
    }
    
    
    # Process each plot
    plot <- plots[[index]]
    plot_size <- plots_sizes[[index]]
    img_tag <- plot2base64(plot, plot_size)
    
    content_with_plots <- paste(content_with_plots, img_tag, sep = "\n")
    pb$tick()
  }
  
  # Close the Table of Contents
  toc <- paste(toc, "</ul></div>", sep="\n")
  
  # Insert the Table of Contents at the placeholder
  content_with_plots <- gsub("<!--TOC-->", toc, content_with_plots)
  
  # Append the final closing tags for the HTML body and document
  content_with_plots <- paste(content_with_plots, "</body></html>", sep="\n")
  
  # Ensure the directory exists
  dir_path <- dirname(output_file_path)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  # Write the HTML content to file
  writeLines(content_with_plots, output_file_path)
  cat("Report written to:", normalizePath(output_file_path), "\n")
}



# Level 2 internal functions ---------------------------------------------------


#' Generate Density Plot
#'
#' @description
#' This function generates a density plot for a given data matrix. The density
#' plot shows the distribution of the values in the data matrix.
#'
#' @param data A numeric matrix containing the data.
#' @param meta A dataframe containing the column meta data of data
#' @param condition The name of the factor column of meta for the experiment
#'
#' @return A ggplot object representing the density plot.
#'
#' @examples
#' \dontrun{
#' # Example usage
#' data <- matrix(rnorm(1000), ncol = 10)
#' density_plot <- make_density_plot(data)
#' print(density_plot)}
#'
#' @importFrom ggplot2 ggplot geom_density ggtitle aes
#' @importFrom reshape2 melt
#'
make_density_plots <- function(data, 
                               meta, 
                               condition) {
  
  density_plots <- list()
  
  # Melt the data to long format
  data_long <- reshape2::melt(as.data.frame(data), id.vars = NULL)
  
  if (length(unique(meta[[condition]])) > 1) {    # Only when > 2 levels
    # Create the overall density plot for all data
    overall_plot <- ggplot2::ggplot(data_long, 
                                    ggplot2::aes(x = !!rlang::sym("value"))) +
      ggplot2::geom_density(fill = "blue", alpha = 0.5) +
      ggplot2::ggtitle("All Levels")
    
    density_plots <- c(density_plots, list(overall_plot))
  }
  
  
  # Create density plots for each level of the condition
  levels <- unique(meta[[condition]])
  for (level in levels) {
    # Filter the data for the current level
    indices <- which(meta[[condition]] == level)
    data_level <- data[, indices, drop = FALSE]
    data_level_long <- reshape2::melt(as.data.frame(data_level), id.vars = NULL)
    
    # Create the density plot for the current level
    level_plot <- ggplot2::ggplot(data_level_long, 
                                  ggplot2::aes(x = !!rlang::sym("value"))) +
      ggplot2::geom_density(fill = "blue", alpha = 0.5) +
      ggplot2::ggtitle(paste("Level:", level))
    
    # Add the level plot to the list
    density_plots <- c(density_plots, list(level_plot))
  }
  
  return(density_plots)
}


#' Generate Boxplot
#'
#' @description
#' This function generates a boxplot for a given data matrix. The boxplot shows 
#' the distribution of the values in the data matrix across different variables, 
#' with log2 intensity on the y-axis.
#'
#' @param data A numeric matrix containing the data.
#' @param meta A dataframe containing the column meta data of data
#' @param condition The name of the factor column of meta for the experiment
#'
#' @return A ggplot object representing the boxplot.
#'
#' @examples
#' \dontrun{
#' # Example usage
#' data <- matrix(rnorm(1000), ncol = 10)
#' boxplot <- make_boxplot(data)
#' print(boxplot)}
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot labs theme element_text
#' @importFrom reshape2 melt
#'
make_box_plots <- function(data,
                           meta,
                           condition) {
  
  boxplots <- list()
  
  # Create boxplots for each level of the condition
  levels <- unique(meta[[condition]])
  for (level in levels) {
    # Filter the data for the current level
    indices <- which(meta[[condition]] == level)
    data_level <- data[, indices, drop = FALSE]
    data_level_long <- reshape2::melt(as.data.frame(data_level), id.vars = NULL)
    
    # Create the boxplot for the current level
    level_plot <- ggplot2::ggplot(data_level_long, 
                                  ggplot2::aes(x = !!rlang::sym("variable"),
                                               y = !!rlang::sym("value"))) +
      ggplot2::geom_boxplot(fill = "#AEC6CF", color = "black") +
      ggplot2::labs(x = "Variables", y = "log2 intensity", 
                    title = paste("Level:", level)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1,
                                                         size = 6))
    
    # Add the level plot to the list
    boxplots <- c(boxplots, list(level_plot))
  }
  
  return(boxplots)
}


#' Generate Violin Plot
#'
#' @description
#' This function generates a violin plot for a given data matrix. The violin 
#' plot shows the distribution of the values in the data matrix across different 
#' variables, with each variable's distribution displayed as a separate violin.
#'
#' @param data A numeric matrix containing the data.
#' @param meta A dataframe containing the column meta data of data
#' @param condition The name of the factor column of meta for the experiment
#'
#' @return A ggplot object representing the violin plot.
#'
#' @examples
#' \dontrun{
#' # Example usage
#' data <- matrix(rnorm(1000), ncol = 10)
#' violin_plot <- make_violin_plot(data)
#' print(violin_plot)}
#' 
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes geom_violin theme labs
#' @importFrom grid unit
#' 
make_violin_plots <- function(data,
                              meta,
                              condition) {
  
  violin_plots <- list()
  
  # Create violin plots for each level of the condition
  levels <- unique(meta[[condition]])
  for (level in levels) {
    # Filter the data for the current level
    indices <- which(meta[[condition]] == level)
    data_level <- data[, indices, drop = FALSE]
    data_level_long <- reshape2::melt(as.data.frame(data_level), id.vars = NULL)
    
    # Create the violin plot for the current level
    level_plot <- ggplot2::ggplot(data_level_long, 
                                  ggplot2::aes(x = !!rlang::sym("variable"),
                                               y = !!rlang::sym("value"))) +
      ggplot2::geom_violin(trim = FALSE, fill = "#77DD77", 
                           color = "black") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1,
                                                         size = 6),
                     plot.margin = grid::unit(c(0, 0, 0, 0), "cm")) +
      ggplot2::labs(x = "Timepoint", y = "Value", 
                    title = paste("Level:", level))
    
    violin_plots <- c(violin_plots, list(level_plot))
  }
  
  return(violin_plots)
}


#' Mean Correlation with Time Plot
#'
#' @description
#' This function takes a data frame with time series data 
#' (rows as features and columns as samples) 
#' and a meta table with sample information including time points, computes 
#' the correlation of each 
#' feature with time, and plots the distribution of these correlations.
#'
#' @param data A data frame where rows are features and columns are samples.
#' @param meta A data frame with sample metadata. Must contain a column "Time".
#'
#' @return A ggplot2 object showing the distribution of mean correlations
#'  with time.
#'  
plot_mean_correlation_with_time <- function(data,
                                            meta,
                                            condition) {

  plot_list <- list()
  
  # Loop through each level of the condition
  for (cond in unique(meta[[condition]])) {
    # Subset the data and meta for the current condition
    condition_indices <- which(meta[[condition]] == cond)
    data_subset <- data[, condition_indices]
    time_subset <- meta$Time[condition_indices]
    
    # Compute correlation of each feature with time
    correlations <- apply(data_subset, 1, function(feature) {
      cor(feature, time_subset, use = "complete.obs")
    })
    
    # Create a data frame for plotting, ensuring row names are set
    if (is.null(rownames(data))) {
      rownames(data) <- paste0("Feature", 1:nrow(data))
    }
    
    cor_data <- data.frame(Feature = rownames(data), Correlation = correlations)
    
    # Generate the plot
    p <- ggplot2::ggplot(cor_data, aes(x = Correlation)) +
      ggplot2::geom_histogram(binwidth = 0.05, fill = "#bcbd22", 
                              color = "black") +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = paste("Level:", cond),
           x = "Correlation with Time",
           y = "Count of Features")
    
    # Add the plot to the list
    plot_list[[cond]] <- p
  }
  
  return(plot_list)
}


#' First Lag Autocorrelation Coefficients Plot
#'
#' @description
#' This function takes a data frame with time series data 
#' (rows as features and columns as samples),
#' a meta table with sample information including time points and conditions, 
#' computes the first lag
#' autocorrelation for each feature for each condition level, and plots the 
#' distribution of these
#' autocorrelation coefficients.
#'
#' @param data A data frame where rows are features and columns are samples.
#' @param meta A data frame with sample metadata. Must contain a column "Time" 
#' and the condition column.
#' @param condition The name of the column in the meta table that contains the 
#' condition information.
#'
#' @return A list of ggplot2 objects, each showing the distribution of first 
#' lag autocorrelation coefficients for one condition.
#'
plot_first_lag_autocorrelation <- function(data,
                                           meta, 
                                           condition) {
  
  # Initialize a list to store the plots
  plot_list <- list()
  
  # Loop through each level of the condition
  for (cond in unique(meta[[condition]])) {
    # Subset the data and meta for the current condition
    condition_indices <- which(meta[[condition]] == cond)
    data_subset <- data[, condition_indices]
    time_subset <- meta$Time[condition_indices]
    
    # Compute first lag autocorrelation of each feature
    autocorrelations <- apply(data_subset, 1, function(feature) {
      # Compute first lag difference
      lag_diff <- diff(feature)
      # Compute autocorrelation
      acf(lag_diff, plot = FALSE)$acf[2]
    })
    
    # Calculate mean and standard deviation of autocorrelations
    mean_autocorrelation <- mean(autocorrelations, na.rm = TRUE)
    std_autocorrelation <- sd(autocorrelations, na.rm = TRUE)
    
    # Create a data frame for plotting
    cor_data <- data.frame(Feature = 1:nrow(data),
                           Autocorrelation = autocorrelations)
    
    # Generate the plot
    p <- ggplot2::ggplot(cor_data, aes(x = Autocorrelation)) +
      ggplot2::geom_histogram(binwidth = 0.05, fill = "#9467bd",
                              color = "black") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = element_text(size = 13),         # Title text size
        axis.title.x = element_text(size = 10),       # X-axis title text size
        axis.title.y = element_text(size = 10),       # Y-axis title text size
        axis.text.x = element_text(size = 7),        # X-axis text size
        axis.text.y = element_text(size = 7)         # Y-axis text size
      ) +
      ggplot2::labs(title = paste("Level:", cond),
           x = "Autocorrelation Coefficient",
           y = "Count of Features",
           subtitle = paste("Mean:", round(mean_autocorrelation, 3), 
                            "SD:", round(std_autocorrelation, 3)))
    
    # Add the plot to the list
    plot_list[[cond]] <- p
  }
  
  return(plot_list)
}


#' Lag-1 Differences Plot
#'
#' @description
#' This function takes a data frame with time series data 
#' (rows as features and columns as samples),
#' a meta table with sample information including time points and conditions, 
#' computes the lag-1
#' differences for each feature for each condition level, and plots the 
#' distribution of these
#' differences.
#'
#' @param data A data frame where rows are features and columns are samples.
#' @param meta A data frame with sample metadata. Must contain a column "Time" 
#' and the condition column.
#' @param condition The name of the column in the meta table that contains the 
#' condition information.
#'
#' @return A list of ggplot2 objects, each showing the distribution of lag-1 
#' differences for one condition.
#'
plot_lag1_differences <- function(data, 
                                  meta, 
                                  condition) {

  plot_list <- list()
  
  # Loop through each level of the condition
  for (cond in unique(meta[[condition]])) {
    # Subset the data and meta for the current condition
    condition_indices <- which(meta[[condition]] == cond)
    data_subset <- data[, condition_indices]
    
    # Compute lag-1 differences of each feature
    lag1_differences <- t(apply(data_subset, 1, 
                                function(feature) {diff(feature)}))
    
    # Calculate mean and stdev of lag-1 differences for each feature
    mean_lag1_diff <- apply(lag1_differences, 1, mean, na.rm = TRUE)
    std_lag1_diff <- apply(lag1_differences, 1, sd, na.rm = TRUE)

    # Create a data frame for plotting
    diff_data <- data.frame(
      Feature = 1:nrow(data),
      Mean_Lag1_Difference = mean_lag1_diff,
      Std_Lag1_Difference = std_lag1_diff
    )
    
    # Generate the plot
    p <- ggplot2::ggplot(diff_data, aes(x = Mean_Lag1_Difference)) +
      ggplot2::geom_histogram(binwidth = 0.05, fill = "#ff7f0e",
                              color = "black") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = element_text(size = 13),         # Title text size
        axis.title.x = element_text(size = 10),       # X-axis title text size
        axis.title.y = element_text(size = 10),       # Y-axis title text size
        axis.text.x = element_text(size = 7),        # X-axis text size
        axis.text.y = element_text(size = 7)         # Y-axis text size
      ) +
      ggplot2::labs(title = paste("Level:", cond),
           x = "Mean Lag-1 Difference",
           y = "Count of Features",
           subtitle = paste("Mean:", 
                            round(mean(mean_lag1_diff, na.rm = TRUE), 3), 
                            "SD:", round(sd(mean_lag1_diff, na.rm = TRUE), 3)))
    
    plot_list[[cond]] <- p
  }
  
  return(plot_list)
}


#' Coefficient of Variation (CV) Plot
#'
#' @description
#' This function takes a data frame with time series data 
#' (rows as features and columns as samples),
#' a meta table with sample information including time points and conditions, 
#' computes the coefficient
#' of variation (CV) for each feature for each condition level, and plots the 
#' distribution of these
#' CVs.
#'
#' @param data A data frame where rows are features and columns are samples.
#' @param meta A data frame with sample metadata. Must contain a column "Time"
#'  and the condition column.
#' @param condition The name of the column in the meta table that contains the 
#' condition information.
#'
#' @return A list of ggplot2 objects, each showing the distribution of CVs for
#'  one condition.
#'
plot_cv <- function(data, 
                    meta, 
                    condition) {
  
  plot_list <- list()
  
  for (cond in unique(meta[[condition]])) {
    condition_indices <- which(meta[[condition]] == cond)
    data_subset <- data[, condition_indices]
    
    # Compute CV of each feature
    cvs <- apply(data_subset, 1, function(feature) {
      sd(feature) / mean(feature)
    })
    
    # Calculate mean and standard deviation of CVs
    mean_cv <- mean(cvs, na.rm = TRUE)
    std_cv <- sd(cvs, na.rm = TRUE)
    
    # Create a data frame for plotting
    cv_data <- data.frame(
      Feature = seq_len(nrow(data)),
      CV = cvs
    )
    
    p <- ggplot2::ggplot(cv_data, aes(x = .data$CV)) +
      ggplot2::geom_histogram(binwidth = 0.05, fill = "#e377c2",
                              color = "black") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = element_text(size = 13),         
        axis.title.x = element_text(size = 10),       
        axis.title.y = element_text(size = 10),       
        axis.text.x = element_text(size = 7),        
        axis.text.y = element_text(size = 7)        
      ) +
      ggplot2::labs(title = paste("Level:", cond),
           x = "Coefficient of Variation (CV)",
           y = "Count of Features",
           subtitle = paste("Mean CV:", round(mean_cv, 3), 
                            "SD CV:", round(std_cv, 3)))
    
    plot_list[[cond]] <- p
  }
  
  return(plot_list)
}


#' Generate PCA Plot with Dynamic Coloring
#'
#' @description
#' This function generates a PCA plot from a data matrix, dynamically coloring
#' the points based on the levels of a specified factor in the metadata.
#'
#' @param data A numeric matrix containing the data.
#' @param meta A dataframe containing the metadata.
#' @param condition The column name in the metadata dataframe that contains the
#' factor levels for coloring the PCA plot.
#'
#' @return A ggplot object representing the PCA plot.
#'
#' @examples
#' \dontrun{
#' # Example usage
#' data <- matrix(rnorm(1000), ncol = 10)
#' meta <- data.frame(Sample = colnames(data), Condition = rep(c("A", "B"), 
#' each = 5))
#' plots <- make_pca_plot(data, meta, "Condition")
#' print(plots)}
#' 
#' @importFrom stats prcomp
#' @importFrom ggplot2 ggplot geom_point xlim xlab ylab ggtitle theme_minimal
#' theme
#' @importFrom ggrepel geom_text_repel
#' 
make_pca_plot <- function(data, 
                          meta, 
                          condition) {
  
  pc <- stats::prcomp(t(data))
  pca_df <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2])
  
  # Get the labels and colors from the metadata
  pca_df$Labels <- colnames(data)
  
  pca_df$Levels <- meta[[condition]]
  
  # Calculate the variance explained
  variance_explained <- pc$sdev^2 / sum(pc$sdev^2)
  percent_variance_explained <- round(variance_explained * 100, digits = 1)
  
  # Calculate the x-axis range and extend it
  x_range <- range(pc$x[, 1])
  extended_x_max <- x_range[2] + (x_range[2] - x_range[1]) * 0.2
  
  # Create the PCA plot
  pca_plot <- ggplot2::ggplot(pca_df, aes(x = !!rlang::sym("PC1"), 
                                          y = !!rlang::sym("PC2"), 
                                          color = !!rlang::sym("Levels"))) +
    ggplot2::geom_point() +
    ggrepel::geom_text_repel(aes(label = !!rlang::sym("Labels")), 
                             box.padding = 0.35, 
                             point.padding = 0.5, 
                             max.overlaps = Inf,
                             size = 2) +
    ggplot2::xlim(x_range[1], extended_x_max) +
    ggplot2::xlab(paste("PC1 -", percent_variance_explained[1], "% variance")) +
    ggplot2::ylab(paste("PC2 -", percent_variance_explained[2], "% variance")) +
    ggplot2::labs(color = "Level") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = element_text(hjust = 0.5))
}


#' Generate PCA Variance Explained Plot
#'
#' @description
#' This function generates a bar plot showing the variance explained by the 
#' first seven principal components of a given data matrix.
#'
#' @param data A numeric matrix containing the data.
#'
#' @return A ggplot object representing the variance explained by the principal 
#' components.
#'
#' @examples
#' \dontrun{
#' # Example usage
#' data <- matrix(rnorm(1000), ncol = 10)
#' variance_plot <- make_pca_variance_explained_plot(data)
#' print(variance_plot)}
#'
#' @importFrom stats prcomp
#' @importFrom limma plotMDS
#' @importFrom ggplot2 ggplot aes geom_col geom_text xlab ylab ggtitle 
#' theme_minimal
#' 
make_pca_variance_explained_plot <- function(data,
                                             meta,
                                             condition) {
  
  # meta and condition are not used, just taken because the functions are 
  # called from another function that passes these arguments.
  
  # Perform PCA
  pc <- stats::prcomp(t(data))
  
  # Calculate the variance explained
  variance_explained <- pc$sdev^2 / sum(pc$sdev^2)
  percent_variance_explained <- round(variance_explained * 100, digits = 1)
  
  # Create a dataframe for the first 7 principal components
  var_explained <- data.frame(
    PC = paste0("PC", 1:7),
    Variance = percent_variance_explained[1:7]
  )
  
  # Create the variance explained plot
  variance_plot <- ggplot2::ggplot(var_explained, 
                                   ggplot2::aes(x = !!rlang::sym("PC"), 
                                                y = !!rlang::sym("Variance"))) +
    ggplot2::geom_col(fill = "#AEC6CF", color = "black") +
    ggplot2::geom_text(aes(label = round(!!rlang::sym("Variance"), digits = 2)), 
                       vjust = -0.2) +
    ggplot2::xlab("Principal Component") +
    ggplot2::ylab("Variance Explained (%)") +
    ggplot2::theme_minimal() +
    ggplot2::ylim(0, max(var_explained$Variance) + 5)  # Space above highest bar
  
  return(variance_plot)
}


#' Generate MDS Plot
#'
#' @description
#' This function generates a multidimensional scaling (MDS) plot for a given 
#' data matrix. The MDS plot visualizes the similarities or dissimilarities 
#' between samples in the data matrix.
#'
#' @param data A numeric matrix containing the data.
#'
#' @return A ggplot object representing the MDS plot.
#'
#' @examples
#' # Example usage
#' data <- matrix(rnorm(1000), ncol = 10)
#' mds_plot <- make_mds_plot(data)
#' print(mds_plot)
#' 
#' @importFrom limma plotMDS
#' @importFrom ggplot2 ggplot geom_point ggtitle theme_minimal
#' @importFrom ggrepel geom_text_repel
#' 
make_mds_plot <- function(data,
                          meta,
                          condition) {

  mds <- limma::plotMDS(x = data, plot = FALSE)
  
  # Extract MDS coordinates
  mds_df <- data.frame(Dim1 = mds$x, 
                       Dim2 = mds$y, 
                       Labels = colnames(data))
  
  mds_df$Levels <- meta[[condition]]
  
  # Generate the MDS plot using ggplot2 and ggrepel
  mds_plot <- ggplot2::ggplot(mds_df, 
                              ggplot2::aes(x = .data$Dim1, 
                                           y = .data$Dim2, 
                                           label = .data$Labels,
                                           color = .data$Levels)) +
    ggplot2::geom_point() +
    ggrepel::geom_text_repel(box.padding = 0.35, 
                             point.padding = 0.5, 
                             max.overlaps = Inf,
                             size = 2) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Dimension 1",
                  y = "Dimension 2",
                  color = "Level") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
}


#' Generate Correlation Heatmaps
#'
#' @description
#' This function generates correlation heatmaps using Spearman correlation for 
#' a given data matrix. It creates a combined heatmap for all levels and 
#' individual heatmaps for each level specified in the condition column of the 
#' metadata.
#'
#' @param data A numeric matrix containing the data.
#' @param meta A dataframe containing the metadata.
#' @param condition The column name in the metadata dataframe that contains the 
#' factor levels for generating individual heatmaps.
#'
#' @return A list of `ComplexHeatmap` heatmap objects representing the 
#' correlation 
#' heatmaps.
#'
#' @examples
#' \dontrun{
#' # Example usage
#' data <- matrix(rnorm(1000), ncol = 10)
#' colnames(data) <- paste0("Sample", 1:10)
#' meta <- data.frame(Sample = colnames(data), Condition = rep(c("A", "B"), 
#' each = 5))
#' heatmaps <- make_correlation_heatmaps(data, meta, "Condition")
#' for (heatmap in heatmaps) {
#'   draw(heatmap, heatmap_legend_side = "right")
#' }}
#'
#' @importFrom stats cor
#' @importFrom circlize colorRamp2
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom grDevices colorRampPalette
#' 
make_correlation_heatmaps <- function(data, 
                                      meta, 
                                      condition) {
  
  heatmaps <- list()
  heatmaps_sizes <- c()
  
  if (length(unique(meta[[condition]])) > 1) {    # Only when > 2 levels
    # Create combined correlation heatmap
    corr_all <- stats::cor(data, method = "spearman", 
                           use = "pairwise.complete.obs")
    # Remove perfect correlations for better visualization
    diag(corr_all) <- NA
    
    # Define the color function for the heatmap based on the range of corr_all
    breaks <- seq(min(corr_all, na.rm = TRUE), 
                  max(corr_all, na.rm = TRUE), length.out = 100)
    col_fun <- 
      circlize::colorRamp2(breaks, 
                           colorRampPalette(c("blue", "white", "red")
                           )(length(breaks)))
    
    heatmap_all <- ComplexHeatmap::Heatmap(
      corr_all,
      col = col_fun,
      name = "Correlation",
      column_title = "All Levels",
      heatmap_legend_param = list(
        title = "Spearman Correlation",
        title_gp = gpar(fontsize = 6),  
        labels_gp = gpar(fontsize = 6),  
        title_position = "leftcenter-rot" 
      ),
      na_col = "grey",
      row_names_gp = gpar(fontsize = 6),
      column_names_gp = gpar(fontsize = 6),
      column_names_rot = 60
    )
    heatmaps <- c(heatmaps, list(heatmap_all))
  }
  
  # Custom scaling logic for the HTML report
  heatmap_all_size <- max(1.5 * length(meta[[condition]]) / 25,
                          1)
  heatmaps_sizes <- c(heatmaps_sizes, heatmap_all_size)
  
  # Create correlation heatmaps for each level of the condition
  levels <- unique(meta[[condition]])
  for (level in levels) {
    # Filter the data for the current level
    indices <- which(meta[[condition]] == level)
    data_level <- data[, indices, drop = FALSE]
    
    # Compute Spearman correlation
    corr_level <- stats::cor(data_level, method = "spearman", 
                             use = "pairwise.complete.obs")
    # Remove perfect correlations for better visualization
    diag(corr_level) <- NA
    
    # Define the color function for the heatmap based on the range of corr_level
    breaks_level <- seq(min(corr_level, na.rm = TRUE), 
                        max(corr_level, na.rm = TRUE), length.out = 100)
    
    col_fun_level <- 
      circlize::colorRamp2(breaks_level, 
                           colorRampPalette(c("blue", "white", "red")
                           )(length(breaks_level)))
    
    # Create the correlation heatmap for the current level
    heatmap_level <- ComplexHeatmap::Heatmap(
      corr_level,
      col = col_fun_level,
      name = "Correlation",
      column_title = paste("Level:", level),
      heatmap_legend_param = list(
        title = "Spearman Correlation",
        title_gp = gpar(fontsize = 6),  
        labels_gp = gpar(fontsize = 6), 
        title_position = "leftcenter-rot"  
      ),
      na_col = "grey",
      row_names_gp = gpar(fontsize = 6),
      column_names_gp = gpar(fontsize = 6),
      column_names_rot = 60
    )
    
    heatmaps <- c(heatmaps, list(heatmap_level))
    
    # Custom scaling logic for the HTML report
    nr_level_timepoints <- sum(meta[[condition]] == level)
    heatmap_level_size <- max(1.5 * nr_level_timepoints / 17,
                              1)
    heatmaps_sizes <- c(heatmaps_sizes, heatmap_level_size)
  }
  
  list(heatmaps = heatmaps,
       heatmaps_sizes = heatmaps_sizes)
}
