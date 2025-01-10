#' Generate scatter plot report for rach feature with embedded base64 images
#'
#' @param data A matrix with features as rows and samples as columns. 
#' Row names should be feature names.
#' @param meta A data frame with the meta information. Must contain a numeric
#'  column "Time".
#' @param output_file The name of the HTML output file.
#'
#' @import ggplot2
#' @import rmarkdown
#' @importFrom progress progress_bar
#'
make_scatter_plot_html <- function(
    data,
    meta,
    output_file = "scatter_plots_base64.html",
    meta_replicate_column = NULL
) {
  # Check if meta contains the Time column
  if (!"Time" %in% colnames(meta)) {
    stop("The meta data must contain a 'Time' column.")
  }
  
  # Ensure the number of columns in 'data' matches the number of rows in 'meta'
  if (ncol(data) != nrow(meta)) {
    stop(
      "The number of columns in 'data' must match the number of rows in 'meta'."
    )
  }
  
  # Check if meta_replicate_column exists in the meta data
  if (!is.null(meta_replicate_column) && !(meta_replicate_column %in% colnames(meta))) {
    stop(paste("The meta data must contain the column:", meta_replicate_column))
  }
  
  # Sort data by feature names
  data <- data[order(rownames(data)), ]
  
  # Create a temporary RMarkdown file for the report
  rmd_file <- base::tempfile(fileext = ".Rmd")
  
  # Initialize content list
  content <- c(
    "---",
    "title: 'Feature Scatter Plots'",
    "output: html_document",
    "---",
    "",
    "## Scatter Plots by Feature",
    ""
  )
  
  # Total number of iterations
  total_iterations <- nrow(data)
  
  # Initialize progress bar with iteration count
  pb <- progress::progress_bar$new(
    format = "  Generating plots [:bar] :current/:total (:percent) in :elapsed",
    total = total_iterations, clear = FALSE, width = 60
  )
  
  # Loop through each feature (row) and generate plots
  for (i in 1:total_iterations) {
    feature_name <- rownames(data)[i]
    feature_values <- data[i, ]
    
    # Check if all values are NA for the feature
    if (all(is.na(feature_values))) {
      content <- c(
        content,
        paste0(
          "### ",
          feature_name,
          " (no plot shown because all values are NA)\n"
        )
      )
    } else {
      # Filter out NA values from feature_values and corresponding Time values
      valid_idx <- !is.na(feature_values)
      valid_values <- feature_values[valid_idx]
      valid_time <- meta$Time[valid_idx]
      
      # Create a valid data frame for ggplot
      valid_data <- data.frame(Time = valid_time, Value = valid_values)
      valid_data$Value <- as.numeric(valid_data$Value)
      
      # Check if replicates are provided and set up coloring
      if (!is.null(meta_replicate_column)) {
        valid_replicates <- meta[[meta_replicate_column]][valid_idx]
        valid_data$Replicate <- factor(valid_replicates)
        aes_params <- ggplot2::aes(x = Time, y = Value, color = Replicate)
        color_label <- meta_replicate_column
      } else {
        aes_params <- ggplot2::aes(x = Time, y = Value)  # No coloring
        color_label <- NULL
      }
      
      p <- ggplot2::ggplot(valid_data, aes_params) +
        ggplot2::geom_point(size = 0.5, alpha = 0.7, shape = 16, na.rm = TRUE) + 
        ggplot2::scale_y_continuous(
          labels = scales::label_number(signif = TRUE, digits = 4),
          breaks = scales::pretty_breaks(n = 5)
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          axis.text = ggplot2::element_text(size = 5),
          axis.title = ggplot2::element_text(size = 6),
          plot.margin = ggplot2::margin(1, 1, 1, 1),
          panel.grid.major.y = ggplot2::element_line(
            size = 0.2,
            linetype = "solid",
            color = "gray"
          ),
          panel.grid.minor.y = ggplot2::element_blank(),
          legend.text = ggplot2::element_text(size = 3),
          legend.title = ggplot2::element_text(size = 3),
          legend.key.size = ggplot2::unit(0.2, "cm")
        )
      
      
      # Add a legend if replicates are provided
      if (!is.null(meta_replicate_column)) {
        p <- p + ggplot2::labs(color = meta_replicate_column)
      }
      
      
      # Convert the plot to base64 using the provided plot2base64 function
      plot_base64 <- plot2base64(
        plot = p,
        height = 1,
        width = 2,
        base_height_per_row = 1
      )
      
      # Add the feature name and plot to the content
      content <- c(
        content,
        paste0(
          "### ",
          feature_name
        ),
        plot_base64,
        "\n"
      )
    }
    
    # Update progress bar
    pb$tick()
  }
  
  # Write the content to the RMarkdown file
  base::writeLines(
    content,
    con = rmd_file
  )
  
  # Render the RMarkdown file to HTML
  rmarkdown::render(
    rmd_file,
    output_file = output_file
  )
  
  # Notify the user
  base::message(
    "Report generated: ",
    output_file
  )
}
