#' utils scripts contains shared functions that are used by at least two package 
#' functions of the SplineOmics package. The level separation is only valid
#' internally in this script, and has no connection to the script level of the
#' respective exported functions scripts.

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
#' \dontrun{
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
#'                      report_dir = here::here())}
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
                                 report_info,
                                 data,
                                 meta = NA,
                                 level_headers_info = NA,
                                 spline_params = NA,
                                 report_type = "explore_data",
                                 mode = NA,
                                 filename = "report",
                                 timestamp = format(Sys.time(), 
                                                    "%d_%m_%Y-%H_%M_%S"),
                                 report_dir = here::here()) {

  if (report_type == "explore_data") {
    if (filename == "explore_data") {
      title <- "explore data"
    } else {
      title <- "explore batch-corrected data"
    }
  } else if (report_type == "limma_hyperparams_screen") {
    title <- paste("hyperparams screen |", filename)
  } else if (report_type == "limma_report") {
    title <- "limma report"
  } else if (report_type == "cluster_hits") {                         
    title <- "clustered hits | within level"
  } else {
    stop(paste("report_type must be explore_hits, limma_hyperparams_screen,", 
               "limma_report, or cluster_hits"),
         call. = FALSE)
  }
  
  header_text <- paste(title, 
                       report_info$omics_data_type, 
                       timestamp, sep=" | ")
  
  if (Sys.getenv("DEVTOOLS_LOAD") == "true") {
    logo_path <- file.path("inst", "extdata", "SplineOmics_logo.png")
  } else {
    logo_path <- system.file("extdata", "SplineOmics_logo.png",
                             package = "SplineOmics")
  }
  
  logo_base64 <- base64enc::dataURI(file = logo_path, mime = "image/png")

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
    "  display: flex;",
    "  align-items: center;",  # Aligns items vertically in the center
    "  justify-content: space-between;", # Ensures the logo is on the far right
    "}",
    ".logo {",
    "  position: absolute;",  # Position the logo absolutely
    "  top: 0;",              # Align with the top of the h1
    "  right: 0;",            # Align with the right of the h1
    "  width: 400px;",  # Adjust the width to make the logo smaller
    "  height: auto;",  # Maintain aspect ratio
    "}",
    "hr {",
    "  margin-top: 20px;",  # Space above the line
    "  margin-bottom: 20px;",  # Space below the line
    "}",
    "</style>",
    "</head><body>",
    "<h1>", header_text, "<img src='",
    logo_base64, "' alt='Logo' class='logo'></h1>",
    "<table>",
    sep=""
  )

  all_fields <- c("omics_data_type",
                  "data_description", 
                  "data_collection_date",
                  "data",
                  "meta",
                  "meta_condition",
                  "meta_batch",
                  "limma_design",
                  "analyst_name", 
                  "contact_info",
                  "project_name",
                  "method_description",
                  "results_summary", 
                  "conclusions")
  
  max_field_length <- max(nchar(gsub("_", " ", all_fields)))

  for (field in all_fields) {
    if (field == "data") {
      value <- sprintf('<a href="%s" download="data.xlsx">
                        <button>Download data.xlsx</button></a>', 
                       encode_df_to_base64(data))
      
    } else if (field == "meta" &&
               !is.null(meta) &&
               is.data.frame(meta) &&
               !any(is.na(meta))) {
      
      value <- sprintf('<a href="%s" download="meta.xlsx">
                        <button>Download meta.xlsx</button></a>', 
                       encode_df_to_base64(meta))
    } else {
      value <- ifelse(is.null(report_info[[field]]), "NA", report_info[[field]])
    }
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
    
  } else if (report_type == "limma_report") {
    
    build_limma_report(header_section = header_section,
                       plots = plots,
                       plots_sizes = plots_sizes,
                       level_headers_info = level_headers_info,
                       output_file_path = output_file_path)
    
  } else {           # report_type == "cluster_hits"
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


#' Encode DataFrame to Base64 for HTML Embedding
#'
#' @description
#' This function takes a dataframe as input and returns a base64 encoded 
#' CSV object. The encoded object can be embedded into an HTML document 
#' directly, with a button to download the file without pointing to a 
#' local file.
#'
#' @param df A dataframe to be encoded.
#' 
#' @return A character string containing the base64 encoded CSV data.
#' 
# encode_df_to_base64 <- function(df) {
#   # Convert dataframe to CSV
#   temp_file <- tempfile(fileext = ".csv")
#   write.csv(df, temp_file, row.names = FALSE)
#   
#   # Read the CSV and encode to base64
#   csv_content <- readBin(temp_file, "raw", file.info(temp_file)$size)
#   base64_csv <- base64enc::base64encode(csv_content)
#   
#   # Create the data URI scheme
#   data_uri <- paste0("data:text/csv;base64,", base64_csv)
#   
#   # Remove the temporary file
#   unlink(temp_file)
#   
#   return(data_uri)
# }
encode_df_to_base64 <- function(df) {
  
  temp_file <- tempfile(fileext = ".xlsx")
  wb <- openxlsx::createWorkbook()
  
  if (is.data.frame(df)) {
    # Convert single dataframe to Excel with one sheet
    sheet_name <- "Sheet1"
    openxlsx::addWorksheet(wb, sheet_name)
    openxlsx::writeData(wb, sheet = sheet_name, df)
  } else if (is.list(df) && all(sapply(df, is.data.frame))) {
    # Convert list of dataframes to Excel with multiple sheets
    sheet_names <- make.unique(names(df))
    
    for (i in seq_along(sheet_names)) {
      openxlsx::addWorksheet(wb, sheet_names[i])
      openxlsx::writeData(wb, sheet = sheet_names[i], df[[i]])
    }
  } else {
    stop("Input must be a dataframe or a list of dataframes.")
  }
  
  openxlsx::saveWorkbook(wb, temp_file, overwrite = TRUE)
  
  # Read the file and encode to base64
  file_content <- readBin(temp_file, "raw", file.info(temp_file)$size)
  base64_file <- base64enc::base64encode(file_content)
  
  # Determine MIME type
  mime_type <- "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
  
  # Create the data URI scheme
  data_uri <- paste0("data:", mime_type, ";base64,", base64_file)
  
  # Remove the temporary file
  unlink(temp_file)
  
  return(data_uri)
}


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
#' \dontrun{
#' library(ggplot2)
#' plot <- ggplot(mtcars, aes(mpg, cyl)) + geom_point()
#' plot2base64(plot, plot_nrows = 1)}
#'
#' @seealso
#' \link[ragg]{agg_png}, \link[base64enc]{dataURI}
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


