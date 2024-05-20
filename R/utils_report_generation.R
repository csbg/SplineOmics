#' utils scripts contains shared functions that are used by at least two package 
#' functions of the SplineOmics package.

# Level 1 internal functions ---------------------------------------------------


#' Generate Report HTML
#'
#' @description
#' Generates an HTML report with the provided plots, spline parameters, and 
#' report information.
#'
#' @param plots A list of ggplot2 plot objects.
#' @param plots_sizes A list of integers specifying the size of each plot.
#' @param level_headers_info A list of header information for each level.
#' @param spline_params A list of spline parameters.
#' @param report_info A named list containing report information.
#' @param report_type A character string specifying the report type 
#'                    ('limma_hyperparams_screen' or 'cluster_hits').
#' @param mode A character string specifying the mode 
#'            ('isolated' or 'integrated').
#' @param filename A character string specifying the filename for the report.
#' @param timestamp A timestamp to include in the report filename.
#' @param report_dir A character string specifying the report directory.
#'
#' @return No return value, called for side effects.
#'
#' @examples
#' plots <- list(ggplot2::ggplot(mtcars, ggplot2::aes(mpg, cyl)) + 
#'                               ggplot2::geom_point())
#' plots_sizes <- list(2)
#' level_headers_info <- list(list("Level 1", 0))
#' spline_params <- list(spline_type = "n", dof = 3)
#' report_info <- list(
#'   omics_data_type = "genomics",
#'   data_description = "Sample description",
#'   data_collection_date = "2023-01-01",
#'   analyst_name = "John Doe",
#'   project_name = "Project XYZ"
#' )
#' generate_report_html(plots, plots_sizes, level_headers_info, spline_params,
#'                      report_info, 
#'                      report_type = "limma_hyperparams_screen", mode = NA, 
#'                      filename = "report", 
#'                      timestamp = format(Sys.time(), "%d_%m_%Y-%H_%M_%S"), 
#'                      report_dir = here::here())
#'
#' @seealso
#' \code{\link{build_hyperparams_screen_report}}, 
#' \code{\link{build_cluster_hits_report}}
#' 
#' @importFrom here here
#' @importFrom tools file_path_sans_ext
#' @importFrom grDevices dev.off
#' @importFrom htmltools save_html
#' 
generate_report_html <- function(plots, 
                                 plots_sizes, 
                                 level_headers_info,
                                 report_info,
                                 spline_params = NA,
                                 report_type = "explore_data",
                                 mode = NA,
                                 filename = "report",
                                 timestamp = format(Sys.time(), 
                                                    "%d_%m_%Y-%H_%M_%S"),
                                 report_dir = here::here()) {
  
  if (report_type == "explore_data") {    
    title <- "explore data"
  } else if (report_type == "limma_hyperparams_screen") {
    title <- "hyperparams screen"
  } else if (report_type == "cluster_hits") {                         
    title <- "clustered hits"
  } else {
    stop(paste("report_type must be explore_hits, limma_hyperparams_screen,", 
               "or cluster_hits"),
         call. = FALSE)
  }
  
  header_text <- paste(title, 
                       report_info$omics_data_type, 
                       timestamp, sep=" | ")
  
  header_section <- paste(
    "<html><head><title>", title, "</title>",
    "<style>",
    "table {",
    "  font-size: 30px;", 
    "}",
    "td {",
    "  padding: 8px;",  # Adds padding to table cells for better readability
    "  text-align: left;",  # Ensures text in cells is left-aligned
    "}",
    "td:first-child {",
    "  text-align: right;",  # Right aligns the first column
    "  color: blue;",        # Sets the color of the first column to blue
    "}",
    "h1 {",
    "  color: #333333;",  # Dark gray color for the title
    "}",
    "hr {",
    "  margin-top: 20px;",  # Space above the line
    "  margin-bottom: 20px;",  # Space below the line
    "}",
    "</style>",
    "</head><body>",
    "<h1>", header_text, "</h1>",
    "<table>",
    sep=""
  )
  
  all_fields <- c("omics_data_type",
                  "data_description", 
                  "data_collection_date",
                  "analyst_name", 
                  "project_name",
                  "dataset_name", 
                  "limma_design",
                  "method_description",
                  "results_summary", 
                  "conclusions",
                  "contact_info")
  
  max_field_length <- max(nchar(gsub("_", " ", all_fields)))
  
  # Add information from the report_info list to the header section
  for (field in all_fields) {
    value <- ifelse(is.null(report_info[[field]]), "NA", report_info[[field]])
    field_display <- sprintf("%-*s", max_field_length, gsub("_", " ", field))
    header_section <- paste(header_section,
                            sprintf('<tr><td style="text-align: right; 
                                    color:blue; padding-right: 5px;">%s 
                                    :</td><td>%s</td></tr>',
                                    field_display, value),
                            sep = "\n")
  }
  
  # Close the table
  header_section <- paste(header_section, "</table>", sep = "\n")
  
  
  file_name <- sprintf("%s_%s_%s.html",
                       filename,
                       report_info$omics_data_type, 
                       timestamp)
  
  output_file_path <- here::here(report_dir, file_name)
  
  if (report_type == "explore_data") {
    
    build_explore_data_report(header_section = header_section, 
                              plots = plots, 
                              plots_sizes = plots_sizes, 
                              output_file_path = output_file_path)
    
  } else if (report_type == "limma_hyperparams_screen") {
    
    build_hyperparams_screen_report(header_section = header_section, 
                                    plots = plots, 
                                    plots_sizes = plots_sizes, 
                                    output_file_path = output_file_path)
    
  }
  else {           # report_type == "cluster_hits"
    build_cluster_hits_report(header_section = header_section, 
                              plots = plots, 
                              plots_sizes = plots_sizes,
                              level_headers_info = level_headers_info,
                              spline_params = spline_params,
                              mode = mode,
                              output_file_path = output_file_path)
  }
}


# Level 2 internal functions ---------------------------------------------------


#' Convert Plot to Base64
#'
#' @description
#' Converts a ggplot2 plot to a Base64-encoded PNG image and returns an HTML 
#' img tag for embedding in a report.
#'
#' @param plot A ggplot2 plot object.
#' @param plot_nrows An integer specifying the number of rows in the plot.
#' @param width A numeric value specifying the width of the plot in inches.
#' @param base_height_per_row A numeric value specifying the base height per
#'  row in inches.
#' @param units A character string specifying the units for the width and 
#' height.
#' @param html_img_width A character string specifying the width of the image
#'  in HTML.
#'
#' @return A character string containing an HTML img tag with the Base64-encoded
#'  plot.
#'
#' @examples
#' library(ggplot2)
#' plot <- ggplot(mtcars, aes(mpg, cyl)) + geom_point()
#' plot2base64(plot, plot_nrows = 1)
#'
#' @seealso
#' \code{\link{ragg::agg_png}}, \code{\link{base64enc::dataURI}}
#' 
#' @importFrom ggplot2 ggsave
#' @importFrom base64enc dataURI
#' @importFrom ragg agg_png
#' 
plot2base64 <- function(plot, 
                        plot_nrows, 
                        width = 7, 
                        base_height_per_row = 2.5, 
                        units = "in", 
                        html_img_width = "100%") {
  
  additional_height_per_row <- 2.1
  height <- base_height_per_row + (plot_nrows - 1) * 
    additional_height_per_row
  
  # Create a graphical device in memory using ragg
  img_file <- tempfile(fileext = ".png")
  ragg::agg_png(filename = img_file, 
                width = width, 
                height = height, 
                units = units, 
                res = 300)
  
  # Draw the plot
  print(plot)
  
  # Turn off the device
  dev.off()
  
  # Convert the image to a Base64 string
  img_base64 <- base64enc::dataURI(file = img_file, mime = "image/png")
  
  # Delete the temporary image file
  unlink(img_file)
  
  # Return the HTML img tag with the Base64 string and a fixed width
  return(sprintf('<img src="%s" alt="Plot" style="width:%s;">', 
                 img_base64, html_img_width))
}


