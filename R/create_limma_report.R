# The function create_limma_report() takes the top_tables of the three different 
# categories (within level time diff, between level average diff, and 
# between level average and time diff) and makes histogram and vulcano plots 
# and places them into a nice HTML report.



# Exported function: create_limma_report() --------------------------------------------


create_limma_report <- function(run_limma_splines_result,
                                report_info,
                                adj_pthresh = 0.05,
                                report_dir = here::here()) {
  
  # Control the function arguments
  args <- lapply(as.list(match.call()[-1]), eval, parent.frame())
  check_null_elements(args)
  input_control <- InputControl$new(args)
  input_control$auto_validate()
  
  # Get the top_tables of the three limma analysis categories
  time_effect <- run_limma_splines_result$time_effect
  avrg_diff_conditions <- run_limma_splines_result$avrg_diff_conditions
  interaction_condition_time <- 
    run_limma_splines_result$interaction_condition_time
  
  plots <- list()
  plots_sizes <- list()
  section_headers_info <- list()
  ###
  
  
  plots <- c(plots, list("Time Effect"))
  plots_sizes <- c(plots_sizes, 999)
  
  header_info <- list(header_name = "Time Effect")
  section_headers_info <- c(section_headers_info, list(header_info))
  
  for (i in seq_along(time_effect)) {
    element_name <- names(time_effect)[i]
    top_table <- time_effect[[i]]
    
    title <- paste("P-Value Histogram:", element_name)
    
    p_value_hist <- create_p_value_histogram(top_table = top_table,
                                             adj_pthresh = adj_pthresh,
                                             title = title)
    
    plots <- c(plots, list(p_value_hist))
    plots_sizes <- c(plots_sizes, 1)
  }
  ###
  
  # length == 0 when there was just one level or no interaction effect
  if (length(avrg_diff_conditions) > 0) {
    
    plots <- c(plots, list("Average Difference Conditions"))
    plots_sizes <- c(plots_sizes, 999)
    
    header_info <- list(header_name = "Average Difference Conditions")
    section_headers_info <- c(section_headers_info, list(header_info))
    
    for (i in seq_along(avrg_diff_conditions)) {
      element_name <- names(avrg_diff_conditions)[i]
      top_table <- avrg_diff_conditions[[i]]
      
      comparison <- remove_prefix(element_name, "avrg_diff_")
      title <- paste("P-Value Histogram:", comparison)
      
      p_value_hist <- create_p_value_histogram(top_table = top_table,
                                               adj_pthresh = adj_pthresh,
                                               title = title)
      
      compared_levels <- stringr::str_split(comparison, "_vs_")[[1]]
      
      volcano_plot <- create_volcano_plot(top_table = top_table,
                                          adj_pthresh = adj_pthresh,
                                          compared_levels)
      
      plots <- c(plots, list(p_value_hist), list(volcano_plot))
      plots_sizes <- c(plots_sizes, 1, 1.5)
    }
  }

  # length == 0 when there was just one level or no interaction effect
  if (length(interaction_condition_time) > 0) {
  
    plots <- c(plots, list("Interaction of Condition and Time"))
    plots_sizes <- c(plots_sizes, 999)
    
    header_info <- list(header_name = "Interaction of Condition and Time")
    section_headers_info <- c(section_headers_info, list(header_info))
    
    for (i in seq_along(interaction_condition_time)) {
      element_name <- names(interaction_condition_time)[i]
      title <- paste(element_name, ": Adj. p-value vs. F-statistic")
      top_table <- interaction_condition_time[[i]]
      
      comparison <- remove_prefix(element_name, "time_interaction_")
      title <- paste("P-Value Histogram:", comparison)
      
      p_value_hist <- create_p_value_histogram(top_table = top_table,
                                               adj_pthresh = adj_pthresh,
                                               title = title)
      
      plots <- c(plots, list(p_value_hist))
      plots_sizes <- c(plots_sizes, 1)
    }
  }
  
  all_top_tables <- c(time_effect,
                      avrg_diff_conditions,
                      interaction_condition_time)
  
  unique_values <- unique(meta[[condition]])
  new_names <- sapply(names(all_top_tables),
                      shorten_names,
                      unique_values = unique_values)
  names(all_top_tables) <- new_names
  
  generate_report_html(plots, 
                       plots_sizes, 
                       report_info,
                       data = all_top_tables,
                       level_headers_info = section_headers_info,
                       report_type = "create_limma_report",
                       filename = "create_limma_report",
                       report_dir = report_dir)
  
  return(plots)
}



# Level 1 internal function definitions ----------------------------------------


#' Create a p-value histogram from a limma top_table
#'
#' @description
#' This function generates a histogram of the unadjusted p-values from a
#' limma top_table.
#'
#' @param top_table A data frame containing the limma top_table with
#'                  a column named `P.Value` for unadjusted p-values.
#' @param adj_pthresh A numeric value for the adjusted p-value threshold
#'                    (not used in this function, included for consistency).
#' @param title A character string for the title of the histogram.
#'
#' @return A ggplot2 object representing the histogram of unadjusted p-values.
#'
#' @importFrom ggplot2 ggplot geom_histogram labs theme_minimal
#'
create_p_value_histogram <- function(top_table,
                                     adj_pthresh = 0.05,
                                     title = "P-Value Histogram") {
  
  # Check if the top_table has a P.Value column
  if (!"P.Value" %in% colnames(top_table)) {
    stop("The top_table must contain a column named 'P.Value'.")
  }
  
  # Create the histogram
  p <- ggplot2::ggplot(top_table, aes(x = P.Value)) +
    ggplot2::geom_histogram(binwidth = 0.05,
                   fill = "orange",
                   color = "black", alpha = 0.7) +
    ggplot2::labs(title = title, x = "Unadjusted P-Value", y = "Frequency") +
    ggplot2::theme_minimal()
  
  return(p)
}


#' Create a Volcano Plot
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
#' 
create_volcano_plot <- function(top_table, 
                                adj_pthresh, 
                                compared_levels) {
  
  # Add a column for coloring points based on logFC
  top_table$Regulation <- ifelse(top_table$logFC > 0, 
                                 compared_levels[2], compared_levels[1])
  
  # Create a named vector for the colors
  colors <- c("blue", "darkgrey")
  names(colors) <- c(compared_levels[2], compared_levels[1])

  # Calculate the number of hits
  num_hits <- sum(top_table$adj.P.Val < adj_pthresh)
  
  volcano_plot <- ggplot2::ggplot(top_table, aes(x = logFC, 
                                                 y = -log10(adj.P.Val), 
                                                 color = Regulation)) +
    ggplot2::geom_point() +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = paste("Volcano Plot:",
                                compared_levels[1], "vs", compared_levels[2]),
                  x = paste("Log Fold Change (",
                            compared_levels[2], " / ",
                            compared_levels[1], ")", sep = ""),
                  y = "-Log10 Adjusted P-value") +
    ggplot2::geom_hline(yintercept = -log10(adj_pthresh), 
                        linetype = "dashed", color = "red") +
    ggplot2::annotate("text", x = Inf, y = -log10(adj_pthresh), 
                      label = paste("Adj.P.Val:", adj_pthresh), 
                      hjust = 1.1, vjust = -1.5, color = "red") +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::annotate("text", x = max(top_table$logFC) * 0.8, y = Inf, 
                      label = paste(compared_levels[2], "higher"), 
                      hjust = 1.1, vjust = 1.1, color = "blue", size = 3) +
    ggplot2::annotate("text", x = min(top_table$logFC) * 0.8, y = Inf, 
                      label = paste(compared_levels[1], "higher"), 
                      hjust = -0.1, vjust = 1.1, color = "darkgrey", size = 3) +
    ggplot2::annotate("text", x = Inf, y = Inf, 
                      label = paste("Hits:", num_hits), 
                      hjust = 1.1, vjust = 2, color = "black", size = 3) +
    ggplot2::theme(legend.position = "none")
}


remove_prefix <- function(string, 
                          prefix) {
  
  pattern <- paste0("^", prefix)
  result <- sub(pattern, "", string)
}


shorten_names <- function(name,
                          unique_values) {
  
  for (val in unique_values) {
    short_val <- substr(val, 1, 3)
    name <- gsub(val, short_val, name, fixed = TRUE)
  }
  return(name)
}


#' Build Cluster Hits Report
#'
#' @description
#' Generates an HTML report for clustered hits, including plots and 
#' spline parameter details, with a table of contents.
#'
#' @param header_section A character string containing the HTML header section.
#' @param plots A list of ggplot2 plot objects.
#' @param plots_sizes A list of integers specifying the size of each plot.
#' @param level_headers_info A list of header information for each level.
#' @param output_file_path A character string specifying the path to save the 
#'                         HTML report.
#'
#' @return No return value, called for side effects.
#'
#' @seealso
#' \code{\link{plot2base64}}, \code{\link{create_progress_bar}}
#' 
build_create_limma_report <- function(header_section, 
                                      plots, 
                                      plots_sizes, 
                                      level_headers_info,
                                      output_file_path = here::here()) {  
  
  html_content <- paste(header_section, "<!--TOC-->", sep = "\n")
  
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
      if (any(class(plots[[index]]) == "character")) {  
        
        section_header <- sprintf("<h2 style='%s' id='section%d'>%s</h2>", 
                                  section_header_style, 
                                  index, 
                                  header_info$header_name)
        
        html_content <- paste(html_content, section_header, 
                                    sep="\n")

        hits_info <- sprintf(
            "<p style='text-align: center; font-size: 30px;'>"
        )
        
        html_content <- paste(html_content, hits_info, sep="\n")
        
        toc_entry <- sprintf("<li style='%s'><a href='#section%d'>%s</a></li>", 
                             toc_style, index, header_info[[1]])
        toc <- paste(toc, toc_entry, sep="\n")
        
        current_header_index <- current_header_index + 1
        
        pb$tick()
        next
      } 
    }

    html_content <- process_plots(plots = plots,
                                  plots_sizes = plots_sizes,
                                  index = index,
                                  html_content = html_content)
    
    pb$tick()
  }
  
  # Close the Table of Contents
  toc <- paste(toc, "</ul></div>", sep="\n")
  
  # Insert the Table of Contents at the placeholder
  html_content <- gsub("<!--TOC-->", toc, html_content)
  
  # Append the final closing tags for the HTML body and document
  html_content <- paste(html_content, "</body></html>", sep="\n")
  
  # Ensure the directory exists
  dir_path <- dirname(output_file_path)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }

  writeLines(html_content, output_file_path)
}
