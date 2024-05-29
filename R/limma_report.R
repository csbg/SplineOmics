# Exported function: limma_report() ---------------------------------------


limma_report <- function(run_limma_splines_result,
                         report_info,
                         adj_pthresh = 0.05,
                         report_dir = here::here()) {
  
  # Get the top_tables of the three limma analysis categories
  within_level <- run_limma_splines_result$within_level_temporal_pattern
  between_level_mean_diff <- run_limma_splines_result$between_level_mean_diff
  between_level_mean_and_temporal_diff <- 
    run_limma_splines_result$between_level_mean_and_temporal_diff
  
  plots <- list()
  plots_sizes <- list()
  section_headers_info <- list()
  ###
  
  
  plots <- c(plots, list("Within Level"))
  plots_sizes <- c(plots_sizes, 999)
  
  header_info <- list(header_name = "Within Level")
  section_headers_info <- c(section_headers_info, list(header_info))
  
  for (i in seq_along(within_level)) {
    element_name <- names(within_level)[i]
    title <- paste(element_name, ": Adj. p-value vs. F-statistic")
    top_table <- within_level[[i]]
    
    f_stat_plot <- create_f_stat_plot(top_table,
                                      adj_pthresh = 0.05,
                                      title = title)
    
    plots <- c(plots, list(f_stat_plot))
    plots_sizes <- c(plots_sizes, 1)
  }
  ###
  
  
  plots <- c(plots, list("Between Level Condition Only"))
  plots_sizes <- c(plots_sizes, 999)
  
  header_info <- list(header_name = "Between Level Condition Only")
  section_headers_info <- c(section_headers_info, list(header_info))
  
  for (i in seq_along(between_level_mean_diff)) {
    element_name <- names(between_level_mean_diff)[i]
    top_table <- between_level_mean_diff[[i]]
    compared_levels <- stringr::str_split(element_name, "_vs_")[[1]]
    
    volcano_plot <- create_volcano_plot(top_table,
                                        adj_pthresh = 0.05,
                                        compared_levels)
    
    plots <- c(plots, list(volcano_plot))
    plots_sizes <- c(plots_sizes, 1.5)
  }
  ###
  
  
  plots <- c(plots, list("Between Level Condition & Time"))
  plots_sizes <- c(plots_sizes, 999)
  
  header_info <- list(header_name = "Between Level Condition & Time")
  section_headers_info <- c(section_headers_info, list(header_info))
  
  for (i in seq_along(between_level_mean_and_temporal_diff)) {
    element_name <- names(between_level_mean_and_temporal_diff)[i]
    title <- paste(element_name, ": Adj. p-value vs. F-statistic")
    top_table <- between_level_mean_and_temporal_diff[[i]]
    
    f_stat_plot <- create_f_stat_plot(top_table,
                                      adj_pthresh = 0.05,
                                      title = title)
    
    plots <- c(plots, list(f_stat_plot))
    plots_sizes <- c(plots_sizes, 1)
  }
  ###
  
  all_top_tables <- c(within_level,
                      between_level_mean_diff,
                      between_level_mean_and_temporal_diff)
  
  generate_report_html(plots, 
                       plots_sizes, 
                       report_info,
                       data = all_top_tables,
                       level_headers_info = section_headers_info,
                       report_type = "limma_report",
                       filename = "limma_report",
                       report_dir = report_dir)
  
  return(plots)
}



# Level 1 internal function definitions ----------------------------------------


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


#' Create an F-statistic plot
#'
#' @description
#' This function creates a scatter plot of F-statistics against the 
#' -log10 adjusted p-values and includes a horizontal dashed line indicating 
#' the adjusted p-value threshold.
#'
#' @param top_table A data frame containing F-statistics and adjusted p-values.
#' @param adj_pthresh A numeric value for the adjusted p-value threshold.
#' @param title A character string for the plot title.
#'
#' @return A ggplot object.
#'
#' @importFrom ggplot2 ggplot aes geom_point theme_minimal geom_hline labs 
#'             annotate
#'             
create_f_stat_plot <- function(top_table, 
                               adj_pthresh, 
                               title) {
  
  num_hits <- sum(top_table$adj.P.Val < adj_pthresh)
  
  y_max <- max(-log10(top_table$adj.P.Val))
  y_hits_position <- y_max * 0.9  
  
  ggplot2::ggplot(top_table, aes(x = F, y = -log10(adj.P.Val))) +
    ggplot2::geom_point(color = "darkgrey", alpha = 0.5) +
    ggplot2::theme_minimal() +
    ggplot2::geom_hline(yintercept = -log10(adj_pthresh), 
                        linetype = "dashed", color = "red") +
    ggplot2::labs(title = title,
                  x = "F-statistic",
                  y = "-Log10 Adjusted P-value") +
    ggplot2::annotate("text", x = Inf, y = -log10(adj_pthresh), 
                      label = paste("Adj.P.Val:", adj_pthresh),
                      hjust = 1.1, vjust = 1.5, color = "red") +
    ggplot2::annotate("text", x = min(top_table$F), y = y_hits_position,
                      label = paste("Hits:", num_hits), 
                      hjust = -0.05, vjust = 1.1, color = "black", size = 3)
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
#' @param spline_params A list of spline parameters.
#' @param mode A character string specifying the mode 
#'            ('isolated' or 'integrated').
#' @param output_file_path A character string specifying the path to save the 
#'                         HTML report.
#'
#' @return No return value, called for side effects.
#'
#' @seealso
#' \code{\link{plot2base64}}, \code{\link{create_progress_bar}}
#' 
build_limma_report <- function(header_section, 
                               plots, 
                               plots_sizes, 
                               level_headers_info,
                               output_file_path) {  
  
  content_with_plots <- paste(header_section, "<!--TOC-->", sep="\n")
  
  toc <- "<div id='toc' style='text-align: center; display: block; margin: auto;
          width: 80%;'> 
        <h2 style='font-size: 40px;'>Table of Contents</h2>
        <ul style='display: inline-block; text-align: left;'>"
  
  section_header_style <- "font-size: 70px; color: #001F3F; text-align: center;"
  toc_style <- "font-size: 30px;"
  
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
        
        content_with_plots <- paste(content_with_plots, section_header, 
                                    sep="\n")

        hits_info <- sprintf(
            "<p style='text-align: center; font-size: 30px;'>"
        )
        
        content_with_plots <- paste(content_with_plots, hits_info, sep="\n")
        
        toc_entry <- sprintf("<li style='%s'><a href='#section%d'>%s</a></li>", 
                             toc_style, index, header_info[[1]])
        toc <- paste(toc, toc_entry, sep="\n")
        
        current_header_index <- current_header_index + 1
        
        pb$tick()
        next
      } 
    }
    
    # Process each plot
    plot <- plots[[index]]
    plot_size <- plots_sizes[[index]]
    img_tag <- plot2base64(plot, plot_size)
    content_with_plots <- paste(content_with_plots, img_tag, sep="\n")
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
}



