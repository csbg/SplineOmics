#' explore_data()
#'
#' @description
#' This function automatically generates exploratory data analysis (EDA) plots 
#' from the provided data. These include density plots, boxplots, PCA plots, 
#' MDS plots, variance explained plots, violin plots, mean correlation with 
#' time, first lag autocorrelation, lag1 differences, and coefficient of 
#' variation. The function returns all EDA plots as a list and, by default, 
#' creates an HTML report containing the plots, saving it to the specified 
#' report directory.
#'
#' @param splineomics A SplineOmics object, containing the data, meta,
#'                    condition, report_info, meta_batch_column, and
#'                    meta_batch2_column.
#' @param report_dir A non-empty string specifying the report directory. The
#'                   output HTML reports will be placed there. Default is the
#'                   current working dir, determined with the here library:
#'                   here::here().
#' @param report A Boolean TRUE or FALSE value, specifying if a report should
#'               be generated or not. A report is generated per default, but
#'               when only the plots as plot objects inside R are desired, this
#'               argument can be set to FALSE.
#'
#' @return A list of ggplot objects representing various exploratory plots.
#'
#' @export
#'
explore_data <- function(
    splineomics,
    report_dir = here::here(),
    report = TRUE
    ) {
  
  report_dir <- normalizePath(
    report_dir,
    mustWork = FALSE
  )

  check_splineomics_elements(
    splineomics = splineomics,
    func_type = "explore_data"
  )

  # Control the function arguments
  args <- lapply(as.list(match.call()[-1]), eval, parent.frame())
  check_null_elements(args)
  input_control <- InputControl$new(args)
  input_control$auto_validate()


  data <- splineomics[["data"]]
  meta <- splineomics[["meta"]]
  annotation <- splineomics[["annotation"]]
  report_info <- splineomics[["report_info"]]
  condition <- splineomics[["condition"]]
  report_info <- splineomics[["report_info"]]
  meta_batch_column <- splineomics[["meta_batch_column"]]
  meta_batch2_column <- splineomics[["meta_batch2_column"]]

  data_list <- list(data = data)

  if (!is.null(meta_batch_column)) {
    args <- list(
      x = data,
      batch = meta[[meta_batch_column]],
      group = meta[[condition]]
    )

    if (!is.null(meta_batch2_column)) {
      args$batch2 <- meta[[meta_batch2_column]]
    }

    batch_corrected_data <- do.call(
      removeBatchEffect,
      args
      )

    data_list$batch_corrected_data <- batch_corrected_data
  }

  all_plots <- list()
  report_info$meta_condition <- c(condition)
  report_info[["plot_data_batch_correction"]] <- paste(
    meta_batch_column,
    meta_batch2_column,
    sep = ", "
  )
  timestamp <- format(Sys.time(), "%d_%m_%Y-%H_%M_%S")

  for (data_name in names(data_list)) {
    current_data <- data_list[[data_name]]
    plots_and_plots_sizes <- generate_explore_plots(
      current_data,
      meta,
      condition
    )

    if (report) {
      generate_report_html(
        plots = plots_and_plots_sizes$plots,
        plots_sizes = plots_and_plots_sizes$plots_sizes,
        report_info = report_info,
        data = data <- bind_data_with_annotation(data, annotation),
        meta = meta,
        filename = paste0("explore_", data_name),
        timestamp = timestamp,
        report_dir = report_dir
      )
    }

    all_plots[[data_name]] <- plots_and_plots_sizes$plots
  }

  print_info_message(
    message_prefix = "Exploratory data analysis",
    report_dir = report_dir
  )

  return(all_plots)
}


# Level 1 internal functions ---------------------------------------------------


#' Generate exploratory plots
#' 
#' @noRd
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
generate_explore_plots <- function(
    data,
    meta,
    condition) {
  meta[[condition]] <- as.factor(meta[[condition]])

  plot_functions_and_sizes <- list(
    list(func = make_density_plots, size = 1),
    list(func = make_violin_box_plots, size = 1.5),
    list(func = make_pca_plot, size = 1.5, flatten = FALSE),
    list(func = make_mds_plot, size = 1.5, flatten = FALSE),
    list(func = make_correlation_heatmaps, size = NULL),
    list(func = plot_mean_correlation_with_time, size = 1.5),
    list(func = plot_lag1_differences, size = 1.5),
    list(func = plot_first_lag_autocorrelation, size = 1.5),
    list(func = plot_cv, size = 1.5)
  )

  apply_plot_function <- function(entry) {
    plot_result <- entry$func(data, meta, condition)
    list(
      plots = plot_result, size = entry$size,
      flatten = if ("flatten" %in% names(entry)) entry$flatten else TRUE
    )
  }

  plot_results <- lapply(
    plot_functions_and_sizes,
    apply_plot_function
    )

  all_plots <- list()
  all_plots_sizes <- c()

  # Flatten the results and sizes conditionally
  for (result in plot_results) {
    all_plots <- c(all_plots, "section_break")
    all_plots_sizes <- c(all_plots_sizes, NA)

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
      all_plots_sizes <- c(
        all_plots_sizes,
        rep(
          result$size,
          length(result$plots)
        )
      )
    }
  }

  list(
    plots = all_plots,
    plots_sizes = all_plots_sizes
  )
}


#' Build Explore Data Report
#' 
#' @noRd
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
#' @param report_info A named list containg the report info fields. Here used
#'                    for the email hotkey functionality.
#' @param output_file_path A string specifying the file path where the HTML
#' report will be saved.
#'
#' @importFrom purrr discard
#'
#' @return None. This function writes the HTML content to the specified file.
#'
build_explore_data_report <- function(
    header_section,
    plots,
    plots_sizes,
    report_info,
    output_file_path) {
  html_content <- paste(header_section, "<!--TOC-->", sep = "\n")

  toc <- create_toc()

  styles <- define_html_styles()
  section_header_style <- styles$section_header_style
  toc_style <- styles$toc_style

  just_plots <- plots |> purrr::discard(~ is.character(.))
  pb <- create_progress_bar(just_plots)

  plot_names <- c(
    "Density Plots",
    "Violin Box Plots",
    "PCA ",
    "MDS",
    "Correlation Heatmaps",
    "Mean Time Correlation",
    "Lag1 Differences",
    "First Lag Autocorrelation",
    "Coefficient of Variation"
  )

  plot_explanations <- get_explore_plots_explanations()

  major_headers <- c(
    "Distribution and Variability Analysis",
    "Dimensionality Reduction and Clustering",
    "Time Series Analysis"
  )

  major_header_style <-
    "font-size: 6em; font-family: Arial, sans-serif; display: inline-block;"

  toc_index <- 0
  toc_index_memory <- toc_index
  major_header_index <- 0

  # Generate the sections and plots
  for (index in seq_along(plots)) {
    if (is.character(plots[[index]])) { # Section break

      toc_index <- toc_index + 1

      if (
        toc_index == 1 || # positions of major headers.
          toc_index == 3 ||
          toc_index == 6
      ) {
        major_header_index <- major_header_index + 1

        section_id <- paste0(
          "section_major_",
          major_header_index
        )
        toc <- paste(
          toc,
          sprintf(
            '<li style="%s"><a href="#%s">%s</a></li>',
            toc_style,
            section_id,
            major_headers[major_header_index]
          ),
          sep = "\n"
        )

        group_header <- sprintf(
          '<div style="text-align: center;"><h1 id="%s"
          style="%s">%s</h1></div>',
          section_id,
          major_header_style,
          major_headers[major_header_index]
        )

        if (toc_index == 3) {
          group_header <- paste(
            group_header,
            '<p style="font-size: 40px;">If you are unsure
            which dimensionality reduction plot to consult,
            choose PCA.</p>',
            sep = "\n"
          )
        }

        html_content <- paste(
          html_content,
          group_header,
          sep = "\n"
        )
      }
      next
    }


    if (toc_index != toc_index_memory) {
      section_id <- paste0("section_", toc_index)
      toc <-
        paste(
          toc,
          sprintf(
            paste0(
              '<li style="margin-left: 20px; font-size: 30px;">',
              '<a href="#%s">%s</a></li>'
            ),
            section_id,
            plot_names[toc_index]
          ),
          sep = "\n"
        )

      section_header <- sprintf(
        '<h2 id="%s" style="%s">%s</h2>',
        section_id,
        section_header_style,
        plot_names[toc_index]
      )

      plot_description <- sprintf(
        '<p style="font-size: 2em;">%s</p>',
        plot_explanations[toc_index]
      )

      html_content <- paste(
        html_content,
        section_header,
        plot_description,
        sep = "\n"
      )

      toc_index_memory <- toc_index
    }

    # Process each plot
    plot <- plots[[index]]
    plot_size <- plots_sizes[[index]]
    img_tag <- plot2base64(plot, plot_size)

    html_content <- paste(
      html_content,
      img_tag,
      "<hr>", # Add horizontal line after each plot
      sep = "\n"
    )
    pb$tick()
  }

  generate_and_write_html(
    toc = toc,
    html_content = html_content,
    report_info = report_info,
    output_file_path = output_file_path
  )
}



# Level 2 internal functions ---------------------------------------------------


#' Generate Density Plot
#' 
#' @noRd
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
#' @importFrom ggplot2 ggplot geom_density ggtitle aes theme
#' @importFrom dplyr everything
#'
make_density_plots <- function(
    data,
    meta,
    condition) {
  custom_theme <- ggplot2::theme(
    panel.background = ggplot2::element_rect(fill = "white", color = "white"),
    plot.background = ggplot2::element_rect(fill = "white", color = "white"),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_line(color = "grey"),
    panel.grid.minor = ggplot2::element_blank()
  )

  density_plots <- list()

  # Melt the data to long format using tidyr
  data_long <- tidyr::pivot_longer(
    as.data.frame(data),
    cols = dplyr::everything(),
    names_to = "variable",
    values_to = "value"
  )


  if (length(unique(meta[[condition]])) > 1) { # Only when > 2 levels
    # Create the overall density plot for all data
    overall_plot <- ggplot2::ggplot(
      data_long,
      ggplot2::aes(x = !!rlang::sym("value"))
    ) +
      ggplot2::geom_density(
        fill = "blue",
        alpha = 0.5
      ) +
      ggplot2::ggtitle("All Levels") +
      custom_theme

    density_plots <- c(
      density_plots,
      list(overall_plot)
    )
  }


  # Create density plots for each level of the condition
  levels <- unique(meta[[condition]])

  for (level in levels) {
    # Filter the data for the current level
    indices <- which(meta[[condition]] == level)

    data_level <- data[, indices, drop = FALSE]

    data_level_long <- tidyr::pivot_longer(
      as.data.frame(data_level),
      cols = everything(),
      names_to = "variable",
      values_to = "value"
    )

    # Create the density plot for the current level
    level_plot <- ggplot2::ggplot(
      data_level_long,
      ggplot2::aes(x = !!rlang::sym("value"))
    ) +
      ggplot2::geom_density(fill = "blue", alpha = 0.5) +
      ggplot2::ggtitle(paste("Level:", level)) +
      custom_theme

    # Add the level plot to the list
    density_plots <- c(
      density_plots,
      list(level_plot)
    )
  }

  return(density_plots)
}



#' Generate Violin Box Plot
#' 
#' @noRd
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
#' @importFrom ggplot2 ggplot aes geom_violin theme labs
#' @importFrom grid unit
#'
make_violin_box_plots <- function(
    data,
    meta,
    condition) {
  custom_theme <- ggplot2::theme(
    panel.background = ggplot2::element_rect(fill = "white", color = "white"),
    plot.background = ggplot2::element_rect(fill = "white", color = "white"),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_line(color = "grey"),
    panel.grid.minor = ggplot2::element_blank()
  )

  plots <- list()

  # Create plots for each level of the condition
  levels <- unique(meta[[condition]])

  for (level in levels) {
    # Filter the data for the current level
    indices <- which(meta[[condition]] == level)

    data_level <- data[, indices, drop = FALSE]

    data_level_long <- tidyr::pivot_longer(
      as.data.frame(data_level),
      cols = everything(),
      names_to = "variable",
      values_to = "value"
    )

    # Create the violin plot with boxplot overlay for the current level
    level_plot <- ggplot2::ggplot(
      data_level_long,
      ggplot2::aes(
        x = !!rlang::sym("variable"),
        y = !!rlang::sym("value")
      )
    ) +
      ggplot2::geom_violin(trim = FALSE, fill = "#77DD77", color = "black") +
      ggplot2::geom_boxplot(
        width = 0.1, fill = "white",
        color = "black", outlier.shape = NA
      ) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(
          angle = 60, hjust = 1,
          size = 6
        ),
        plot.margin = grid::unit(c(0, 0, 0, 0), "cm")
      ) +
      ggplot2::labs(
        x = "Timepoint", y = "Value",
        title = paste("Level:", level)
      ) +
      custom_theme

    # Add the level plot to the list
    plots <- c(plots, list(level_plot))
  }

  return(plots)
}


#' Mean Correlation with Time Plot
#' 
#' @noRd
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
#' @param condition The column of the meta dataframe containign the levels that
#'                  separate the experiment.
#'
#' @return A ggplot2 object showing the distribution of mean correlations
#'  with time.
#'
#'  @importFrom rlang .data
#'
plot_mean_correlation_with_time <- function(
    data,
    meta,
    condition
    ) {
  
  plot_list <- list()

  # Loop through each level of the condition
  for (cond in unique(meta[[condition]])) {
    # Subset the data and meta for the current condition
    condition_indices <- which(meta[[condition]] == cond)
    data_subset <- data[, condition_indices]
    time_subset <- meta$Time[condition_indices]

    # Compute correlation of each feature with time
    correlations <- apply(data_subset, 1, function(feature) {
      # Check if either feature or time_subset has zero variance
      if (sd(feature, na.rm = TRUE) == 0 
          || sd(time_subset, na.rm = TRUE) == 0) {
        return(NA) # Assign NA if variance is zero
      } else {
        return(cor(
          feature,
          time_subset,
          use = "complete.obs"
          ))
      }
    })

    # Create a data frame for plotting, ensuring row names are set
    if (is.null(rownames(data))) {
      rownames(data) <- paste0("Feature", 1:nrow(data))
    }

    cor_data <- data.frame(
      Feature = rownames(data),
      Correlation = correlations
    )
    
    # Remove non-finite values from the correlation data
    cor_data <- cor_data[is.finite(cor_data$Correlation), ]

    # Generate the plot
    p <- ggplot2::ggplot(cor_data, aes(x = .data$Correlation)) +
      ggplot2::geom_histogram(
        binwidth = 0.05, fill = "#bcbd22",
        color = "black"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = paste("Level:", cond),
        x = "Correlation with Time",
        y = "Count of Features"
      )

    # Add the plot to the list
    plot_list[[cond]] <- p
  }

  return(plot_list)
}


#' First Lag Autocorrelation Coefficients Plot
#' 
#' @noRd
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
#' @importFrom stats acf
#' @importFrom rlang .data
#'
plot_first_lag_autocorrelation <- function(
    data,
    meta,
    condition
    ) {
  
  # Initialize a list to store the plots
  plot_list <- list()

  # Loop through each level of the condition
  for (cond in unique(meta[[condition]])) {
    # Subset the data and meta for the current condition
    condition_indices <- which(meta[[condition]] == cond)
    data_subset <- data[, condition_indices]
    # time_subset <- meta$Time[condition_indices]

    # Compute first lag autocorrelation of each feature
    autocorrelations <- apply(data_subset, 1, function(feature) {
      # Compute first lag difference
      lag_diff <- diff(feature)
      # Compute autocorrelation
      stats::acf(lag_diff, plot = FALSE)$acf[2]
    })
    
    # Remove non-finite values from autocorrelations
    valid_indices <- which(is.finite(autocorrelations)) # Find valid rows
    autocorrelations <- autocorrelations[valid_indices]
    data_subset <- data_subset[valid_indices, ] # Subset data to valid rows

    # Calculate mean and standard deviation of autocorrelations
    mean_autocorrelation <- mean(autocorrelations, na.rm = TRUE)
    std_autocorrelation <- sd(autocorrelations, na.rm = TRUE)

    # Create a data frame for plotting
    cor_data <- data.frame(
      Feature = valid_indices,
      Autocorrelation = autocorrelations
    )

    # Generate the plot
    p <- ggplot2::ggplot(cor_data, aes(x = .data$Autocorrelation)) +
      ggplot2::geom_histogram(
        binwidth = 0.05, fill = "#9467bd",
        color = "black"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = element_text(size = 13), # Title text size
        axis.title.x = element_text(size = 10), # X-axis title text size
        axis.title.y = element_text(size = 10), # Y-axis title text size
        axis.text.x = element_text(size = 7), # X-axis text size
        axis.text.y = element_text(size = 7) # Y-axis text size
      ) +
      ggplot2::labs(
        title = paste("Level:", cond),
        x = "Autocorrelation Coefficient",
        y = "Count of Features",
        subtitle = paste(
          "Mean:", round(mean_autocorrelation, 3),
          "SD:", round(std_autocorrelation, 3)
        )
      )

    # Add the plot to the list
    plot_list[[cond]] <- p
  }

  return(plot_list)
}


#' Lag-1 Differences Plot
#' 
#' @noRd
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
#' @importFrom stats sd
#' @importFrom rlang .data
#'
plot_lag1_differences <- function(
    data,
    meta,
    condition) {
  plot_list <- list()

  # Loop through each level of the condition
  for (cond in unique(meta[[condition]])) {
    # Subset the data and meta for the current condition
    condition_indices <- which(meta[[condition]] == cond)
    data_subset <- data[, condition_indices]

    # Compute absolute lag-1 differences of each feature
    lag1_differences <- t(apply(
      data_subset, 1,
      function(feature) {
        abs(diff(feature))
      }
    ))

    # Normalize lag-1 differences by the mean of the feature values
    feature_means <- apply(
      data_subset,
      1,
      mean,
      na.rm = TRUE
    )

    normalized_lag1_differences <- lag1_differences / feature_means

    # Calculate mean and stdev of normalized lag-1 differences for each feature
    mean_lag1_diff <- apply(
      normalized_lag1_differences,
      1,
      mean,
      na.rm = TRUE
    )

    std_lag1_diff <- apply(
      normalized_lag1_differences,
      1,
      stats::sd,
      na.rm = TRUE
    )
    
    # Identify rows with finite mean_lag1_diff
    valid_indices <- which(is.finite(mean_lag1_diff))
    
    # Subset all related data to ensure consistency
    mean_lag1_diff <- mean_lag1_diff[valid_indices]
    std_lag1_diff <- std_lag1_diff[valid_indices]
    
    binwidth <- max(0.05, diff(range(mean_lag1_diff, na.rm = TRUE)) / 30)

    # Create a data frame for plotting
    diff_data <- data.frame(
      Feature = valid_indices,
      Mean_Lag1_Difference = mean_lag1_diff,
      Std_Lag1_Difference = std_lag1_diff
    )

    # Generate the plot
    p <- ggplot2::ggplot(
      diff_data,
      aes(x = .data$Mean_Lag1_Difference)
    ) +
      ggplot2::geom_histogram(
        binwidth = binwidth,
        fill = "#ff7f0e",
        color = "black"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 13),
        axis.title.x = ggplot2::element_text(size = 10),
        axis.title.y = ggplot2::element_text(size = 10),
        axis.text.x = ggplot2::element_text(size = 7),
        axis.text.y = ggplot2::element_text(size = 7)
      ) +
      ggplot2::labs(
        title = paste("Level:", cond),
        x = "Mean Normalized Absolute Lag-1 Difference",
        y = "Count of Features",
        subtitle = paste(
          "Mean:",
          round(mean(mean_lag1_diff, na.rm = TRUE), 3),
          "SD:", round(stats::sd(mean_lag1_diff, na.rm = TRUE), 3)
        )
      )

    plot_list[[cond]] <- p
  }

  return(plot_list)
}



#' Coefficient of Variation (CV) Plot
#' 
#' @noRd
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
plot_cv <- function(
    data,
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
    std_cv <- stats::sd(cvs, na.rm = TRUE)
    
    # Remove non-finite values from CVs and adjust related data
    valid_indices <- which(is.finite(cvs)) # Find valid rows
    cvs <- cvs[valid_indices]             # Subset CVs to valid rows
    data_subset <- data_subset[valid_indices, ] # Subset data to valid rows
    
    binwidth <- max(0.01, diff(range(cvs, na.rm = TRUE)) / 30)
    
    # Create a data frame for plotting
    cv_data <- data.frame(
      Feature = valid_indices,
      CV = cvs
    )

    p <- ggplot2::ggplot(cv_data, aes(x = .data$CV)) +
      ggplot2::geom_histogram(
        binwidth = binwidth,
        fill = "#e377c2",
        color = "black"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = element_text(size = 13),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7)
      ) +
      ggplot2::labs(
        title = paste("Level:", cond),
        x = "Coefficient of Variation (CV)",
        y = "Count of Features",
        subtitle = paste(
          "Mean CV:",
          round(mean_cv, 3),
          "SD CV:",
          round(std_cv, 3)
        )
      )

    plot_list[[cond]] <- p
  }

  return(plot_list)
}


#' Generate PCA Plot with Dynamic Coloring
#' 
#' @noRd
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
#' @importFrom stats prcomp
#' @importFrom ggplot2 ggplot geom_point xlim xlab ylab ggtitle theme_minimal
#' theme
#' @importFrom ggrepel geom_text_repel
#'
make_pca_plot <- function(
    data,
    meta,
    condition) {
  # Perform PCA
  pc <- stats::prcomp(t(data))
  pca_df <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2])

  # Add labels and levels from the metadata
  pca_df$Labels <- colnames(data)
  pca_df$Levels <- meta[[condition]]

  # Add time column for alpha transparency
  pca_df$Time <- meta$Time

  # Normalize the Time column to a 0-1 range for alpha values
  pca_df$Alpha <- scales::rescale(pca_df$Time, to = c(0.3, 1))

  # Calculate the variance explained
  variance_explained <- pc$sdev^2 / sum(pc$sdev^2)
  percent_variance_explained <- round(variance_explained * 100, digits = 1)

  # Extend the x-axis range
  x_range <- range(pc$x[, 1])
  extended_x_max <- x_range[2] + (x_range[2] - x_range[1]) * 0.2

  # Create the PCA plot
  pca_plot <- ggplot2::ggplot(pca_df, aes(
    x = !!rlang::sym("PC1"),
    y = !!rlang::sym("PC2"),
    color = !!rlang::sym("Levels"),
    alpha = !!rlang::sym("Alpha")
  )) +
    ggplot2::geom_point() +
    ggrepel::geom_text_repel(aes(label = !!rlang::sym("Labels")),
      box.padding = 0.35,
      point.padding = 0.5,
      max.overlaps = Inf,
      size = 2
    ) +
    ggplot2::xlim(x_range[1], extended_x_max) +
    ggplot2::xlab(paste("PC1 -", percent_variance_explained[1], "% variance")) +
    ggplot2::ylab(paste("PC2 -", percent_variance_explained[2], "% variance")) +
    ggplot2::labs(color = condition) +
    ggplot2::scale_alpha(range = c(0.3, 1), guide = "none") + # Hide alpha
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = element_text(hjust = 0.5))

  return(pca_plot)
}



#' Generate MDS Plot
#' 
#' @noRd
#'
#' @description
#' This function generates a multidimensional scaling (MDS) plot for a given
#' data matrix. The MDS plot visualizes the similarities or dissimilarities
#' between samples in the data matrix.
#'
#' @param data A numeric matrix containing the data.
#' @param meta A dataframe, containign the meta information of data.
#' @param condition The column of the meta dataframe containign the levels that
#'                  separate the experiment.
#'
#' @return A ggplot object representing the MDS plot.
#'
#' @importFrom limma plotMDS
#' @importFrom ggplot2 ggplot geom_point ggtitle theme_minimal
#' @importFrom ggrepel geom_text_repel
#'
make_mds_plot <- function(
    data,
    meta,
    condition) {
  # Perform MDS using limma's plotMDS function
  mds <- limma::plotMDS(
    x = data,
    plot = FALSE
  )

  # Extract MDS coordinates
  mds_df <- data.frame(
    Dim1 = mds$x,
    Dim2 = mds$y,
    Labels = colnames(data)
  )

  # Add condition levels and time information from metadata
  mds_df$Levels <- meta[[condition]]
  mds_df$Time <- meta$Time

  # Normalize the Time column to a 0-1 range for alpha values
  mds_df$Alpha <- scales::rescale(mds_df$Time, to = c(0.3, 1))

  # Generate the MDS plot using ggplot2 and ggrepel
  mds_plot <- ggplot2::ggplot(
    mds_df,
    ggplot2::aes(
      x = .data$Dim1,
      y = .data$Dim2,
      label = .data$Labels,
      color = .data$Levels,
      alpha = .data$Alpha # Use alpha for transparency
    )
  ) +
    ggplot2::geom_point() +
    ggrepel::geom_text_repel(
      box.padding = 0.35,
      point.padding = 0.5,
      max.overlaps = Inf,
      size = 2
    ) +
    ggplot2::scale_alpha(range = c(0.3, 1), guide = "none") + # Hide alpha
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "Dimension 1",
      y = "Dimension 2",
      color = condition
    ) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  return(mds_plot)
}


#' Generate Correlation Heatmaps
#' 
#' @noRd
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
#' @importFrom stats cor
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom grDevices colorRampPalette
#'
make_correlation_heatmaps <- function(
    data,
    meta,
    condition) {
  heatmaps <- list()
  heatmaps_sizes <- c()

  if (length(unique(meta[[condition]])) > 1) { # Only when > 2 levels
    # Create combined correlation heatmap
    corr_all <- stats::cor(data,
      method = "spearman",
      use = "pairwise.complete.obs"
    )
    # Remove perfect correlations for better visualization
    diag(corr_all) <- NA

    # Define the color function for the heatmap based on the range of corr_all
    breaks <- seq(min(corr_all, na.rm = TRUE),
      max(corr_all, na.rm = TRUE),
      length.out = 100
    )

    col_fun <- colorRampPalette(c("blue", "white", "red"))

    heatmap_all <- ComplexHeatmap::Heatmap(
      corr_all,
      col = col_fun(100),
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
  heatmap_all_size <- max(
    1.5 * length(meta[[condition]]) / 25,
    1
  )
  heatmaps_sizes <- c(heatmaps_sizes, heatmap_all_size)

  # Create correlation heatmaps for each level of the condition
  levels <- unique(meta[[condition]])
  for (level in levels) {
    # Filter the data for the current level
    indices <- which(meta[[condition]] == level)
    data_level <- data[, indices, drop = FALSE]

    # Compute Spearman correlation
    corr_level <- stats::cor(data_level,
      method = "spearman",
      use = "pairwise.complete.obs"
    )
    # Remove perfect correlations for better visualization
    diag(corr_level) <- NA

    # Define the color function for the heatmap based on the range of corr_level
    breaks_level <- seq(min(corr_level, na.rm = TRUE),
      max(corr_level, na.rm = TRUE),
      length.out = 100
    )

    col_fun_level <- colorRampPalette(c("blue", "white", "red"))

    # Create the correlation heatmap for the current level
    heatmap_level <- ComplexHeatmap::Heatmap(
      corr_level,
      col = col_fun_level(100),
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
    heatmap_level_size <- max(
      1.5 * nr_level_timepoints / 17,
      1
    )
    heatmaps_sizes <- c(heatmaps_sizes, heatmap_level_size)
  }

  list(
    heatmaps = heatmaps,
    heatmaps_sizes = heatmaps_sizes
  )
}


#' Get Plot Explanations
#' 
#' @noRd
#'
#' @description
#' This function returns a vector of text explanations for various types of
#' plots. These explanations are used in HTML reports to describe the plots.
#'
#' @return A character vector containing explanations for different plot types.
#'
#' @details
#' The explanations cover a variety of plots, including density plots, boxplots,
#' violin plots, mean time correlation plots, lag-1 differences plots, first lag
#' autocorrelation plots, coefficient of variation (CV) plots, PCA plots, PCA
#' variance explained plots, MDS plots, and correlation heatmaps. Each
#' explanation provides insights on what the plot shows and how to interpret it.
#'
get_explore_plots_explanations <- function() {
  plot_explanations_file <- system.file(
    "descriptions",
    "explore_plot_explanations.txt",
    package = "SplineOmics"
  )
  plot_explanations <- readLines(plot_explanations_file)
}
