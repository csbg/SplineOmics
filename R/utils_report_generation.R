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
#' @param data A dataframe or a list of dataframes, containing data that should
#'             be directly embedded in the HTML report for downloading.
#' @param meta A dataframe, containing metadata that should
#'             be directly embedded in the HTML report for downloading.
#' @param report_type A character string specifying the report type 
#'                    ('screen_limma_hyperparams' or 'cluster_hits').
#' @param mode A character string specifying the mode 
#'            ('isolated' or 'integrated').
#' @param filename A character string specifying the filename for the report.
#' @param timestamp A timestamp to include in the report filename.
#' @param report_dir A character string specifying the report directory.
#'
#' @return No return value, called for side effects.
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
                                 data = NA,
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
  } else if (report_type == "screen_limma_hyperparams") {
    title <- paste("hyperparams screen |", filename)
  } else if (report_type == "create_limma_report") {
    title <- "limma report"
  } else if (report_type == "cluster_hits") {                         
    title <- "clustered hits | within level"
  } else if (report_type == "create_gsea_report") {                         
    title <- "gsea"
    
  } else {
    stop(paste("report_type must be explore_hits, screen_limma_hyperparams,", 
               "create_limma_report, or cluster_hits"),
         call. = FALSE)
  }
  
  fields_to_format <- c("data_description",
                        "method_description",
                        "results_summary",
                        "conclusions")
  
  for (field in fields_to_format) {
    if (field %in% names(report_info)) {
      report_info[[field]] <- format_text(report_info[[field]])
    }
  }
  
  header_text <- paste(title, 
                       report_info$omics_data_type, 
                       timestamp, sep = " | ")
  
  header_section <- get_header_section(title = title,
                                       header_text = header_text)

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
    
    if (field == "data" && !any(is.na(data)))  {
      
      if (report_type == "create_limma_report") {
        base64_df <- sprintf('<a href="%s" download="top_tables.xlsx">
                            <button>Download top_tables.xlsx</button></a>', 
                            encode_df_to_base64(data))
      } else {
        base64_df <- sprintf('<a href="%s" download="data.xlsx">
                             <button>Download data.xlsx</button></a>', 
                             encode_df_to_base64(data))
      }

      
    } else if (field == "meta" &&
               !is.null(meta) &&
               is.data.frame(meta) &&
               !any(is.na(meta))) {
      
      base64_df <- sprintf('<a href="%s" download="meta.xlsx">
                           <button>Download meta.xlsx</button></a>', 
                           encode_df_to_base64(meta))
    } else {
      base64_df <- ifelse(is.null(report_info[[field]]),
                          "NA", report_info[[field]])
    }
    field_display <- sprintf("%-*s", max_field_length, gsub("_", " ", field))
    header_section <- paste(header_section,
                            sprintf('<tr><td style="text-align: right;
                                    color:blue; padding-right: 5px;">%s
                                    :</td><td>%s</td></tr>',
                                    field_display, base64_df),
                            sep = "\n")
  }
  
  # Close the table
  header_section <- paste(header_section, "</table>", sep = "\n")
  
  if (report_type == "create_gsea_report") {
    databases_text <- paste(report_info$databases, collapse = ", ")
    header_section <- paste(header_section, 
                            "<p style='font-size: 20px;'>Databases used: ", 
                            databases_text, 
                            "</p>", 
                            sep = "\n")
  }
  
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
    
  } else if (report_type == "screen_limma_hyperparams") {
    
    build_hyperparams_screen_report(header_section = header_section, 
                                    plots = plots, 
                                    plots_sizes = plots_sizes, 
                                    output_file_path = output_file_path)
    
  } else if (report_type == "create_limma_report") {
    
    build_create_limma_report(header_section = header_section,
                       plots = plots,
                       plots_sizes = plots_sizes,
                       level_headers_info = level_headers_info,
                       output_file_path = output_file_path)
    
  } else if (report_type == "create_gsea_report") {
    
    build_create_gsea_report(header_section = header_section,
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


#' Format text
#'
#' @description
#' This function takes a character vector `text` and splits it into individual
#' characters. It then iterates over the characters and builds lines not exceeding
#' a specified character limit (default 70). Newlines are inserted between lines
#' using the `<br>` tag, suitable for HTML display.
#' 
#' @param text A character vector to be formatted.
#' 
#' @return A character vector with formatted text containing line breaks.
#' 
format_text <- function(text) {

  letters <- strsplit(text, "")[[1]]
  formatted_lines <- vector(mode = "character", length = 0)  
  current_line <- ""
  for (char in letters) {
    if (nchar(current_line) + nchar(char) <= 70) {
      current_line <- paste(current_line, char, sep = "")
    } else {
      formatted_lines <- c(formatted_lines, current_line)
      current_line <- char
    }
  }
  formatted_lines <- c(formatted_lines, current_line) 
  formatted_text <- paste(formatted_lines, collapse = "<br>")
}


#' Get Header Section
#'
#' @description
#' Generates the HTML header section for a report, including the title, header 
#' text, and logo. This section also includes the styling for the table and 
#' other HTML elements.
#'
#' @param title A string specifying the title of the HTML document.
#' @param header_text A string specifying the text to be displayed in the 
#' header of the report.
#'
#' @return A string containing the HTML header section.
#'
#' @details
#' The function checks the `DEVTOOLS_LOAD` environment variable to determine 
#' the path to the logo image. The logo image is then converted to a base64 
#' data URI and included in the HTML. The header section includes styles for 
#' tables, table cells, and header elements to ensure proper formatting and 
#' alignment.
#'
#' @importFrom base64enc dataURI
#' 
get_header_section <- function(title,
                               header_text) {
  
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
    sep = ""
  )
}


#' Encode DataFrame to Base64 for HTML Embedding
#'
#' @description
#' This function takes a dataframe as input and returns a base64 encoded 
#' CSV object. The encoded object can be embedded into an HTML document 
#' directly, with a button to download the file without pointing to a 
#' local file.
#'
#' @param df A dataframe to be encoded.
#' @param report_type (Optional) A string specifying for which report generation
#'                    this function is called. Generates different Excel sheet
#'                    names based on the report_type.
#' 
#' @return A character string containing the base64 encoded CSV data.
#' 
#' @importFrom openxlsx createWorkbook addWorksheet writeData saveWorkbook
#' 
encode_df_to_base64 <- function(df,
                                report_type = NA) {
  
  temp_file <- tempfile(fileext = ".xlsx")
  wb <- openxlsx::createWorkbook()
  
  if (is.data.frame(df)) {
    # Convert single dataframe to Excel with one sheet
    sheet_name <- "Sheet1"
    openxlsx::addWorksheet(wb, sheet_name)
    openxlsx::writeData(wb, sheet = sheet_name, df)
  } else if (is.list(df) && all(sapply(df, is.data.frame))) {
    # Convert list of dataframes to Excel with multiple sheets
    
    if (!is.na(report_type)) {
      if (report_type == "create_gsea_report") {
        all_names <- names(df)
        sheet_names <- sapply(all_names, extract_and_combine)
      }
    } else {
      sheet_names <- make.unique(names(df))
    }
    
    
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


#' Create Table of Contents
#'
#' @description
#' Creates the HTML content for the Table of Contents.
#'
#' @return A string containing the HTML for the Table of Contents.
#' 
create_toc <- function() {
  
  toc <- "<div id='toc' style='text-align: center; display: block; margin: auto;
          width: 80%;'>
        <h2 style='font-size: 40px;'>Table of Contents</h2>
        <ul style='display: inline-block; text-align: left;'>"
}


#' Define HTML Styles
#'
#' @description
#' Defines the CSS styles for section headers and Table of Contents (TOC) 
#' entries used in the GSEA report generation.
#'
#' @return A list containing the styles for section headers and TOC entries.
#' 
define_html_styles <- function() {
  
  section_header_style <- "font-size: 70px; color: #001F3F; text-align: center;"
  toc_style <- "font-size: 30px;"
  
  styles <- list(
    section_header_style = section_header_style,
    toc_style = toc_style
  )
}


#' Process Plots
#'
#' @description
#' Converts plots to base64 and appends them to the HTML content.
#'
#' @param plots A list of plots to be processed.
#' @param plots_sizes A list of sizes for the plots.
#' @param index An integer, specifying which element index to take of plots and
#'              plots_sizes.
#' @param html_content The current state of the HTML content.
#'
#' @return Updated HTML content with the plots included.
#' 
process_plots <- function(plots,
                          plots_sizes,
                          index,
                          html_content) {

  img_tag <- plot2base64(plots[[index]], plots_sizes[[index]])
  html_content <- paste(html_content, img_tag, sep = "\n")
}



# Level 3 internal functions ---------------------------------------------------


extract_and_combine <- function(input) {
  
  # Extract substring after 'cluster:' and before the next ','
  cluster_match <- str_extract(input, "(?<=cluster: )[^,]+")
  # Extract substring after 'database:'
  database_match <- str_extract(input, "(?<=database: )[^,]+")
  
  # Combine the substrings with a whitespace in between
  combined <- paste(cluster_match, database_match, sep = " ")
  
  # Truncate the combined string to 30 characters if necessary
  if (nchar(combined) > 30) {
    combined <- substr(combined, 1, 30)
  }
  
  return(combined)
}


