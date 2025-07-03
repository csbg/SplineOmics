#' create_limma_report()
#'
#' @description
#' Generates an HTML report based on the results of a limma analysis with
#' splines.
#' The report includes various plots and sections summarizing the analysis
#' results for time effects, average differences between conditions,
#' and interaction effects between condition and time.
#'
#' @param splineomics An S3 object of class `SplineOmics` that contains all the
#' necessary data and parameters for the analysis, including:
#' \itemize{
#'   \item \code{limma_splines_result}: A list containing top tables from
#'    differential expression analysis for the three different limma results.
#'   \item \code{meta}: A data frame with sample metadata. Must contain a column
#'                      "Time".
#'   \item \code{condition}: A character string specifying the column name in
#'                          the metadata (\code{meta}) that defines groups
#'                          for analysis. This column contains levels such as
#'                           "exponential" and "stationary" for phases, or
#'                           "drug" and "no_drug" for treatments.
#'   \item \code{annotation}: A data frame containing feature information,
#'                            such as gene and protein names, associated with
#'                            the expression data.
#'   \item \code{report_info}: A list containing metadata about the analysis
#'                             for reporting purposes.
#' }
#' @param adj_pthresh A numeric value specifying the adjusted p-value threshold
#' for significance. Default is 0.05. Must be > 0 and < 1.
#' @param report_dir A string specifying the directory where the report should
#' be saved. Default is the current working directory.
#'
#' @return A list of plots included in the generated HTML report.
#'
#' @importFrom here here
#'
#' @export
#'
create_limma_report <- function(
    splineomics,
    adj_pthresh = 0.05,
    report_dir = here::here()
    ) {
  
  report_dir <- normalizePath(
    report_dir,
    mustWork = FALSE
  )

  check_splineomics_elements(
    splineomics = splineomics,
    func_type = "create_limma_report"
  )

  # Control the function arguments
  args <- lapply(as.list(match.call()[-1]), eval, parent.frame())
  check_null_elements(args)
  input_control <- InputControl$new(args)
  input_control$auto_validate()

  limma_splines_result <- splineomics[["limma_splines_result"]]
  meta <- splineomics[["meta"]]
  condition <- splineomics[["condition"]]
  annotation <- splineomics[["annotation"]]
  report_info <- splineomics[["report_info"]]
  
  # Put them in there under those names, so that the report generation fun
  # can access them directly like this.
  design <- splineomics[["design"]] 
  effects <- extract_effects(design)
  report_info[["Fixed effects"]] <- effects[["fixed_effects"]]
  report_info[["Random effects"]] <- effects[["random_effects"]] 

  # Get the top_tables of the three limma analysis categories
  time_effect <- limma_splines_result$time_effect
  avrg_diff_conditions <- limma_splines_result$avrg_diff_conditions
  interaction_condition_time <- limma_splines_result$interaction_condition_time

  plots <- list()
  plots_sizes <- list()
  section_headers_info <- list()


  result <- generate_time_effect_plots(
    time_effect,
    adj_pthresh
  )

  plots <- c(
    plots,
    result$plots
  )
  plots_sizes <- c(
    plots_sizes,
    result$plots_sizes
  )
  section_headers_info <- c(
    section_headers_info,
    result$section_headers_info
  )


  # length == 0 when there was just one level or no interaction effect
  if (length(avrg_diff_conditions) > 0) {
    result <- generate_avrg_diff_plots(
      avrg_diff_conditions,
      adj_pthresh
    )

    plots <- c(
      plots,
      result$plots
    )
    plots_sizes <- c(
      plots_sizes,
      result$plots_sizes
    )
    section_headers_info <- c(
      section_headers_info,
      result$section_headers_info
    )
  }

  # length == 0 when there was just one level or no interaction effect
  if (length(interaction_condition_time) > 0) {
    result <- generate_interaction_plots(
      interaction_condition_time,
      adj_pthresh
    )

    plots <- c(
      plots,
      result$plots
    )
    plots_sizes <- c(
      plots_sizes,
      result$plots_sizes
    )
    section_headers_info <- c(
      section_headers_info,
      result$section_headers_info
    )
  }

  all_top_tables <- c(
    time_effect,
    avrg_diff_conditions,
    interaction_condition_time
  )

  unique_values <- unique(meta[[condition]])
  new_names <- vapply(
    names(all_top_tables),
    shorten_names,
    unique_values = unique_values,
    FUN.VALUE = character(1)
  )
  names(all_top_tables) <- new_names

  if (!is.null(annotation)) {
    # Add annotation info into the top_tables
    for (index in seq_along(all_top_tables)) {
      all_top_tables[[index]] <- merge_top_table_with_annotation(
        top_table = all_top_tables[[index]],
        annotation = annotation
      )
    }
  }

  generate_report_html(
    plots,
    plots_sizes,
    report_info,
    topTables = all_top_tables,
    level_headers_info = section_headers_info,
    report_type = "create_limma_report",
    filename = "limma_report",
    report_dir = report_dir
  )

  print_info_message(
    message_prefix = "Limma report generation",
    report_dir = report_dir
  )

  return(plots)
}



# Level 1 internal function definitions ----------------------------------------


#' Generate Plots for Time Effect
#' 
#' @noRd
#'
#' @description
#' Creates p-value histograms for each time effect in the LIMMA analysis. This
#' function is used internally in the `create_limma_report` function.
#'
#' @param time_effect A list of top tables from the LIMMA analysis representing
#' the time effects.
#' @param adj_pthresh A numeric value specifying the adjusted p-value threshold
#' for significance.
#'
#' @return A list containing the plots and their sizes, as well as the
#' section header information.
#'
generate_time_effect_plots <- function(
    time_effect,
    adj_pthresh) {
  plots <- list("Time Effect")
  plots_sizes <- c(999)

  header_info <- list(header_name = "Time Effect")
  section_headers_info <- list(header_info)

  for (i in seq_along(time_effect)) {
    element_name <- names(time_effect)[i]
    top_table <- time_effect[[i]]

    title <- paste("P-Value Histogram:", element_name)

    p_value_hist <- create_p_value_histogram(
      top_table = top_table,
      title = title
    )

    plots <- c(
      plots,
      list(p_value_hist)
    )
    plots_sizes <- c(plots_sizes, 1)
  }

  list(
    plots = plots,
    plots_sizes = plots_sizes,
    section_headers_info = section_headers_info
  )
}


#' Generate Plots for Average Difference Conditions
#' 
#' @noRd
#'
#' @description
#' Creates p-value histograms and volcano plots for each condition in the
#' average difference conditions. This function is used internally in the
#' `create_limma_report` function.
#'
#' @param avrg_diff_conditions A list of top tables from the LIMMA analysis
#' representing the average difference between conditions.
#' @param adj_pthresh A numeric value specifying the adjusted p-value threshold
#' for significance.
#'
#' @return A list containing the plots and their sizes, as well as the
#' section header information.
#'
generate_avrg_diff_plots <- function(
    avrg_diff_conditions,
    adj_pthresh) {
  plots <- list("Average Difference Conditions")
  plots_sizes <- c(999)

  header_info <- list(header_name = "Average Difference Conditions")
  section_headers_info <- list(header_info)

  for (i in seq_along(avrg_diff_conditions)) {
    element_name <- names(avrg_diff_conditions)[i]
    top_table <- avrg_diff_conditions[[i]]

    comparison <- remove_prefix(element_name, "avrg_diff_")
    title <- paste("P-Value Histogram:", comparison)

    p_value_hist <- create_p_value_histogram(
      top_table = top_table,
      title = title
    )

    compared_levels <- strsplit(comparison, "_vs_")[[1]]

    volcano_plot <- create_volcano_plot(
      top_table = top_table,
      adj_pthresh = adj_pthresh,
      compared_levels
    )

    plots <- c(plots, list(p_value_hist), list(volcano_plot))
    plots_sizes <- c(plots_sizes, 1, 1.5)
  }

  list(
    plots = plots,
    plots_sizes = plots_sizes,
    section_headers_info = section_headers_info
  )
}


#' Generate Plots for Interaction of Condition and Time
#' 
#' @noRd
#'
#' @description
#' Creates p-value histograms for each interaction condition in the
#' interaction of condition and time. This function is used internally in the
#' `create_limma_report` function.
#'
#' @param interaction_condition_time A list of top tables from the LIMMA
#' analysis
#' representing the interaction effects between condition and time.
#' @param adj_pthresh A numeric value specifying the adjusted p-value threshold
#' for significance.
#'
#' @return A list containing the plots and their sizes, as well as the
#' section header information.
#'
generate_interaction_plots <- function(
    interaction_condition_time,
    adj_pthresh) {
  plots <- list("Interaction of Condition and Time")
  plots_sizes <- c(999)

  header_info <- list(header_name = "Interaction of Condition and Time")
  section_headers_info <- list(header_info)

  for (i in seq_along(interaction_condition_time)) {
    element_name <- names(interaction_condition_time)[i]
    title <- paste(element_name, ": Adj. p-value vs. F-statistic")
    top_table <- interaction_condition_time[[i]]

    comparison <- remove_prefix(element_name, "time_interaction_")
    title <- paste("P-Value Histogram:", comparison)

    p_value_hist <- create_p_value_histogram(
      top_table = top_table,
      title = title
    )

    plots <- c(plots, list(p_value_hist))
    plots_sizes <- c(plots_sizes, 1)
  }

  list(
    plots = plots,
    plots_sizes = plots_sizes,
    section_headers_info = section_headers_info
  )
}


#' Shorten Names
#' 
#' @noRd
#'
#' @description
#' Replaces occurrences of unique values within a name with their first three
#' characters. This function is useful for abbreviating long condition names
#' in a dataset.
#'
#' @param name A string representing the name to be shortened.
#' @param unique_values A vector of unique values whose abbreviations will
#' replace their occurrences in the name.
#'
#' @return A string with the unique values replaced by their abbreviations.
#'
shorten_names <- function(
    name,
    unique_values) {
  for (val in unique_values) {
    short_val <- substr(val, 1, 3)
    name <- gsub(val, short_val, name, fixed = TRUE)
  }
  return(name)
}


#' Build Cluster Hits Report
#' 
#' @noRd
#'
#' @description
#' Generates an HTML report for clustered hits, including plots and
#' spline parameter details, with a table of contents.
#'
#' @param header_section A character string containing the HTML header section.
#' @param plots A list of ggplot2 plot objects.
#' @param plots_sizes A list of integers specifying the size of each plot.
#' @param level_headers_info A list of header information for each level.
#' @param report_info A named list containg the report info fields. Here used
#'                    for the email hotkey functionality.
#' @param output_file_path A character string specifying the path to save the
#'                         HTML report.
#'
#' @return No return value, called for side effects.
#'
#' @seealso
#' \code{\link{plot2base64}}, \code{\link{create_progress_bar}}
#'
build_create_limma_report <- function(
    header_section,
    plots,
    plots_sizes,
    level_headers_info,
    report_info,
    output_file_path = here::here()
    ) {
  
  # Read the text file and split it into blocks
  descriptions_path <- system.file(
    "descriptions",
    "create_limma_report_html_plot_descriptions.txt",
    package = "SplineOmics"
  )
  text_blocks <- readLines(descriptions_path)
  text_blocks <- split(
    text_blocks,
    cumsum(text_blocks == "")
  ) # Split by empty lines

  # Remove empty elements created by split
  text_blocks <- Filter(function(x) length(x) > 0, text_blocks)


  html_content <- paste(
    header_section,
    "<!--TOC-->",
    sep = "\n"
  )

  toc <- create_toc()

  styles <- define_html_styles()
  section_header_style <- styles$section_header_style
  toc_style <- styles$toc_style

  current_header_index <- 1
  level_headers_info <- Filter(Negate(is.null), level_headers_info)

  pb <- create_progress_bar(plots)
  # Generate the sections and plots
  for (index in seq_along(plots)) {
    if (current_header_index <= length(level_headers_info)) {
      header_info <- level_headers_info[[current_header_index]]

      # means jump to next section
      if (any(is(plots[[index]], "character"))) {
        section_header <- sprintf(
          "<h2 style='%s' id='section%d'>%s</h2>",
          section_header_style,
          index,
          header_info$header_name
        )

        html_content <- paste(
          html_content,
          section_header,
          sep = "\n"
        )

        # Add text block after the main headers.
        if (current_header_index <= length(text_blocks)) {
          # Convert text block to a single string with HTML paragraph format
          text_block_html <- paste(
            "<p style='font-size: 200%;'>",
            paste(
              text_blocks[[current_header_index]],
              collapse = " "
            ),
            "</p>",
            sep = ""
          )

          html_content <- paste(
            html_content,
            text_block_html,
            sep = "\n"
          )
        }

        hits_info <- sprintf(
          "<p style='text-align: center; font-size: 30px;'>"
        )

        html_content <- paste(
          html_content,
          hits_info,
          sep = "\n"
        )

        toc_entry <- sprintf(
          "<li style='%s'><a href='#section%d'>%s</a></li>",
          toc_style,
          index,
          header_info[[1]]
        )

        toc <- paste(
          toc,
          toc_entry,
          sep = "\n"
        )

        current_header_index <- current_header_index + 1

        pb$tick()
        next
      }
    }

    result <- process_plots(
      plots_element = plots[[index]],
      plots_size = plots_sizes[[index]],
      html_content = html_content,
      toc = toc,
      header_index = current_header_index
    )
    html_content <- result$html_content
    toc <- result$toc

    pb$tick()
  }

  generate_and_write_html(
    toc = toc,
    html_content = html_content,
    report_info = report_info,
    output_file_path = output_file_path
  )
}



# Level 2 internal function definitions ----------------------------------------


#' Create a p-value histogram from a limma top_table
#' 
#' @noRd
#'
#' @description
#' This function generates a histogram of the unadjusted p-values from a
#' limma top_table.
#'
#' @param top_table A data frame containing the limma top_table with
#'                  a column named `P.Value` for unadjusted p-values.
#' @param title A character string for the title of the histogram.
#'
#' @return A ggplot2 object representing the histogram of unadjusted p-values.
#'
#' @importFrom ggplot2 ggplot geom_histogram labs theme_minimal
#' @importFrom rlang .data
#'
create_p_value_histogram <- function(
    top_table,
    title = "P-Value Histogram") {
  # Check if the top_table has a P.Value column
  if (!"P.Value" %in% colnames(top_table)) {
    stop("The top_table must contain a column named 'P.Value'.")
  }

  # Create the histogram
  p <- ggplot2::ggplot(
    top_table,
    # aes(x = P.Value)
    aes(x = .data$P.Value)
  ) +
    ggplot2::geom_histogram(
      binwidth = 0.01,
      fill = "orange",
      color = "black", alpha = 0.7
    ) +
    ggplot2::labs(
      title = title,
      x = "Unadjusted P-Value",
      y = "Frequency"
    ) +
    ggplot2::theme_minimal()

  return(p)
}


#' Create a Volcano Plot
#' 
#' @noRd
#'
#' @description
#' This function creates a volcano plot from a limma top table, plotting
#' log fold changes against the negative log10 of adjusted p-values.
#'
#' @param top_table A data frame from limma containing 'logFC' and 'adj.P.Val'
#' columns.
#' @param adj_pthresh A numeric value for the adjusted p-value threshold.
#' @param compared_levels A character vector of length 2 specifying the
#' compared levels.
#'
#' @return A ggplot object representing the volcano plot.
#'
#' @importFrom ggplot2 ggplot aes geom_point theme_minimal labs geom_hline
#'                     annotate scale_color_manual
#' @importFrom rlang .data
#'
create_volcano_plot <- function(
    top_table,
    adj_pthresh,
    compared_levels
    ) {

  logFC <- NULL
  adj.P.Val <- NULL

  top_table <- top_table |>
    dplyr::mutate(
      Regulation = ifelse(
        logFC > 0,
        compared_levels[2],
        compared_levels[1]
      ),
      Alpha = ifelse(
        adj.P.Val < adj_pthresh,
        1,
        0.5
      )
    )

  # Create a named vector for the colors
  colors <- c("blue", "darkgrey") |>
    setNames(c(
      compared_levels[2],
      compared_levels[1]
    ))

  # Calculate the number of hits
  num_hits <- sum(top_table$adj.P.Val < adj_pthresh)
  total_tests <- nrow(top_table)
  percent_hits <- (num_hits / total_tests) * 100

  volcano_plot <- ggplot2::ggplot(
    data = top_table,
    mapping = ggplot2::aes(
      x = logFC,
      y = -log10(adj.P.Val),
      color = .data$Regulation,
      alpha = .data$Alpha
    )
  ) +
    # Add a shaded region below the adjusted p-value threshold
    ggplot2::annotate(
      "rect",
      xmin = -Inf,
      xmax = Inf,
      ymin = 0,
      ymax = -log10(adj_pthresh),
      fill = "gray90" # Light gray color for visual distinction
    ) +
    ggplot2::geom_point() +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = paste(
        "Volcano Plot:",
        compared_levels[1],
        "vs",
        compared_levels[2]
      ),
      x = paste(
        "Log2 Fold Change (",
        compared_levels[2],
        " / ",
        compared_levels[1],
        ")",
        sep = ""
      ),
      y = "-Log10 Adjusted P-value"
    ) +
    # Add a dashed red line for the adjusted p-value threshold
    ggplot2::geom_hline(
      yintercept = -log10(adj_pthresh),
      linetype = "dashed",
      color = "red"
    ) +
    # Add a label for the adjusted p-value threshold
    ggplot2::annotate(
      geom = "text",
      x = Inf,
      y = -log10(adj_pthresh),
      label = paste(
        "adj.p-thresh:",
        adj_pthresh
      ),
      hjust = 1.1,
      vjust = -1.5,
      color = "red"
    ) +
    ggplot2::scale_color_manual(
      values = colors
    ) +
    ggplot2::scale_alpha(
      range = c(0.5, 1)
    ) +
    ggplot2::annotate(
      geom = "text",
      x = max(top_table$logFC) * 0.8,
      y = Inf,
      label = paste(
        compared_levels[2],
        "higher"
      ),
      hjust = 1.1,
      vjust = 1.1,
      color = "blue",
      size = 3
    ) +
    ggplot2::annotate(
      geom = "text",
      x = min(top_table$logFC) * 0.8,
      y = Inf,
      label = paste(
        compared_levels[1],
        "higher"
      ),
      hjust = -0.1,
      vjust = 1.1,
      color = "darkgrey",
      size = 3
    ) +
    ggplot2::annotate(
      geom = "text",
      x = 0,
      y = Inf,
      label = paste(
        "Hits:",
        num_hits,
        "(", round(percent_hits, 2), "%)"
      ),
      hjust = 0.5,
      vjust = 3,
      color = "black",
      size = 3
    ) +
    ggplot2::theme(
      legend.position = "none"
    )

  return(volcano_plot)
}


#' Remove Prefix from String
#' 
#' @noRd
#'
#' @description
#' Removes a specified prefix from the beginning of a string. This function
#' is useful for cleaning or standardizing strings by removing known prefixes.
#'
#' @param string A string from which the prefix should be removed.
#' @param prefix A string representing the prefix to be removed.
#'
#' @return A string with the prefix removed.
#'
remove_prefix <- function(
    string,
    prefix) {
  pattern <- paste0("^", prefix)
  result <- sub(pattern, "", string)
}
