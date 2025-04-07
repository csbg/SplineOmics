#' make_scatter_plot_html()
#' 
#' @description
#' This function is used to make scatter plots for the raw data of all the 
#' features. It generates an HTML report in the fashion of the other functions
#' of the SplineOmics package which contains all the scatter plots.
#'
#' @param data A matrix with features as rows and samples as columns. 
#' Row names should be feature names.
#' @param meta A data frame with the meta information. Must contain a numeric
#'  column "Time".
#' @param output_file The name of the HTML output file.
#' @param meta_replicate_column Column name of the column in meta that contains
#' the info about the replicates, such as reactor.
#'
#' @import ggplot2
#' @seealso \code{\link[rmarkdown]{render}}
#' @importFrom progress progress_bar
#' @examples
#' \dontrun{
#' # Example Data
#' data <- matrix(rnorm(50), nrow = 5)
#' meta <- data.frame(Time = seq(1, 10, length.out = 10))
#'
#' # Generate HTML report (only if you want to test it)
#' make_scatter_plot_html(data, meta, "scatter_report.html")
#' }
#' @export
#'
make_scatter_plots_html <- function(
    data,
    meta,
    output_file = "scatter_report",  # no .html suffix anymore
    meta_replicate_column = NULL,
    features_per_file = 500
) {
  
  if (!"Time" %in% colnames(meta)) {
    stop_call_false("The meta data must contain a 'Time' column.")
  }
  
  if (ncol(data) != nrow(meta)) {
    stop_call_false(
      "The number of columns in 'data' must match the number of rows in 'meta'."
    )
  }
  
  if (!is.null(meta_replicate_column) &&
      !(meta_replicate_column %in% colnames(meta))) {
    stop_call_false(paste(
      "The meta data must contain the column:",
      meta_replicate_column
    ))
  }
  
  data <- data[order(rownames(data)), ]
  total_features <- nrow(data)
  total_chunks <- ceiling(total_features / features_per_file)
  
  message("Total features: ", total_features)
  message(
    "Generating ",
    total_chunks,
    " HTML reports in chunks of ",
    features_per_file
    )
  
  for (chunk_index in seq_len(total_chunks)) {
    feature_start <- (chunk_index - 1) * features_per_file + 1
    feature_end <- min(chunk_index * features_per_file, total_features)
    
    chunk_data <- data[feature_start:feature_end, , drop = FALSE]
    
    rmd_file <- tempfile(fileext = ".Rmd")
    con <- file(rmd_file, open = "wt")
    
    # Write RMarkdown header
    writeLines(c(
      "---",
      paste0("title: 'Feature Scatter Plots (Chunk ", chunk_index, ")'"),
      "output: html_document",
      "---",
      "",
      paste0("## Features ", feature_start, " to ", feature_end),
      ""
    ), con)
    
    pb <- progress::progress_bar$new(
      format = paste0(
        "  Chunk ",
        chunk_index,
        " [:bar] :current/:total (:percent) in :elapsed"
        ),
      total = nrow(chunk_data), clear = FALSE, width = 60
    )
    
    for (i in seq_len(nrow(chunk_data))) {
      feature_name <- rownames(chunk_data)[i]
      feature_values <- chunk_data[i, ]
      
      if (all(is.na(feature_values))) {
        writeLines(
          paste0(
            "### ",
            feature_name,
            " (no plot shown because all values are NA)\n"
            ),
          con
          )
      } else {
        valid_idx <- !is.na(feature_values)
        valid_values <- feature_values[valid_idx]
        valid_time <- meta$Time[valid_idx]
        valid_data <- data.frame(Time = valid_time, Value = valid_values)
        valid_data$Value <- as.numeric(valid_data$Value)
        
        if (!is.null(meta_replicate_column)) {
          valid_replicates <- meta[[meta_replicate_column]][valid_idx]
          valid_data$Replicate <- factor(valid_replicates)
          aes_params <- ggplot2::aes(x = Time, y = Value, color = Replicate)
        } else {
          aes_params <- ggplot2::aes(x = Time, y = Value)
        }
        
        p <- ggplot2::ggplot(valid_data, aes_params) +
          ggplot2::geom_point(
            size = 0.5,
            alpha = 0.7,
            shape = 16,
            na.rm = TRUE
            ) +
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
        
        if (!is.null(meta_replicate_column)) {
          p <- p + ggplot2::labs(color = meta_replicate_column)
        }
        
        plot_base64 <- plot2base64(
          plot = p,
          height = 1,
          width = 2,
          base_height_per_row = 1
        )

        writeLines(paste0("### ", feature_name), con)
        writeLines("```{=html}", con)
        writeLines(plot_base64, con)
        writeLines("```", con)
        writeLines("", con)
        
      }
      pb$tick()
    }
    
    close(con)
    
    chunk_output <- paste0(output_file, "_chunk_", chunk_index, ".html")
    rmarkdown::render(rmd_file, output_file = chunk_output, envir = new.env())
    message("✔ Rendered: ", chunk_output)
  }
  
  message("✅ All reports generated successfully.")
}

