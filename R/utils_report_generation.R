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
generate_report_html <- function(
    plots, 
    plots_sizes, 
    report_info,
    data = NULL,
    meta = NA,
    topTables = NA,
    enrichr_format = NA,
    level_headers_info = NA,
    spline_params = NA,
    report_type = "explore_data",
    analysis_type = NA,  # only for cluster_hits()
    mode = NA,
    filename = "report",
    timestamp = format(Sys.time(), 
                      "%d_%m_%Y-%H_%M_%S"),
    report_dir = here::here()
    ) {
  
  feature_names_formula <- NA
  
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
    title <- paste("clustered hits |", analysis_type)
    feature_names_formula <- generate_feature_name_template(data)
    
  } else if (report_type == "create_gsea_report") {                         
    title <- "gsea"
    
  } else {
    stop(paste("report_type must be explore_hits, screen_limma_hyperparams,", 
               "create_limma_report, or cluster_hits"),
         call. = FALSE)
  }
  
  fields_to_format <- c(
    "data_description",
    "method_description",
    "results_summary",
    "conclusions"
    )
  
  for (field in fields_to_format) {
    if (field %in% names(report_info)) {
      report_info[[field]] <- format_text(report_info[[field]])
    }
  }

  header_text <- paste(
    title, 
    paste("Omics-Datatype:", report_info$omics_data_type), 
    paste("Date-Time:", timestamp), sep = " | "
    )
  header_text <- paste(header_text, "<br><br><br>")
  
  header_section <- get_header_section(
    title = title,
    header_text = header_text,
    report_type = report_type,
    feature_names_formula = feature_names_formula
    )
  
  report_info_fields <- c(
    "omics_data_type",
    "data_description", 
    "data_collection_date",
    "meta_condition",
    "meta_batch",
    "limma_design",
    "analyst_name", 
    "contact_info",
    "project_name",
    "method_description",
    "results_summary", 
    "conclusions"
    )

  download_fields <- c()
  if (!is.null(data)) download_fields <- c(
    download_fields,
    "data_with_annotation"
    )
  
  if (!all(is.na(meta))) download_fields <- c(
    download_fields,
    "meta"
    )
  
  if (!all(is.na(topTables))) download_fields <- c(
    download_fields,
    "limma_topTables"
    )
  
  if (!all(is.na(enrichr_format))) {
    download_fields <- c(
      download_fields,
      "Enrichr_clustered_genes",
      "Enrichr_background"
      )
  }
  
  max_field_length <- max(nchar(gsub("_", " ", report_info_fields)))
  
  report_info_section <- paste(
    '<hr style="border: none; height: 3px;', 
    'background-color: #333; margin: 40px 0;', 
    'width: 75%;"> <h2 style="font-size: 48px;">', 
    'Report Info ℹ️</h2><table>', sep = ""
    )
  
  downloads_section <- paste(
    '<hr><h2 style="font-size: 48px;">', 
    'Downloads 📥</h2><table>'
    )

  for (field in report_info_fields) {

    base64_df <- process_field(
      field,
      data,
      meta,
      topTables,
      report_info,
      encode_df_to_base64,
      report_type, 
      enrichr_format
      )
    
    field_display <- sprintf("%-*s", max_field_length, gsub("_", " ", field))
    report_info_section <- paste(
      report_info_section,
      sprintf('<tr><td style="text-align: right;
           color:blue; padding-right: 5px;">%s
           :</td><td>%s</td></tr>',
             field_display, base64_df),
      sep = "\n"
      )
  }
  
  # Close the Report Info table
  report_info_section <- paste(report_info_section, "</table>", sep = "\n")
  
  for (field in download_fields) {

    base64_df <- process_field(
      field,
      data,
      meta,
      topTables,
      report_info,
      encode_df_to_base64,
      report_type,
      enrichr_format
      )
    
    field_display <- sprintf("%-*s", max_field_length, gsub("_", " ", field))
    downloads_section <- paste(
      downloads_section,
      sprintf('<tr><td style="text-align: right;
           color:blue; padding-right: 5px;">%s
           :</td><td>%s</td></tr>',
             field_display, base64_df),
      sep = "\n"
      )
  }
  
  # Close the Downloads table
  downloads_section <- paste(downloads_section, "</table>", sep = "\n")
  
  # Preserve initial header_section content
  header_section <- paste(header_section,
                          report_info_section,
                          downloads_section, sep = "\n")
  
  
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
    
    build_explore_data_report(
      header_section = header_section, 
      plots = plots, 
      plots_sizes = plots_sizes, 
      output_file_path = output_file_path
      )
    
  } else if (report_type == "screen_limma_hyperparams") {
    
    build_hyperparams_screen_report(
      header_section = header_section, 
      plots = plots, 
      plots_sizes = plots_sizes, 
      output_file_path = output_file_path
      )
    
  } else if (report_type == "create_limma_report") {
    
    build_create_limma_report(
      header_section = header_section,
      plots = plots,
      plots_sizes = plots_sizes,
      level_headers_info = level_headers_info,
      output_file_path = output_file_path
      )
    
  } else if (report_type == "create_gsea_report") {
    
    build_create_gsea_report(
      header_section = header_section,
      plots = plots,
      plots_sizes = plots_sizes,
      level_headers_info = level_headers_info,
      output_file_path = output_file_path
      )
    
  } else {           # report_type == "cluster_hits"
    build_cluster_hits_report(
      header_section = header_section, 
      plots = plots, 
      plots_sizes = plots_sizes,
      level_headers_info = level_headers_info,
      spline_params = spline_params,
      mode = mode,
      output_file_path = output_file_path
      )
  }
}


#' Generate and Write HTML Content with Table of Contents
#'
#' @description
#' This function generates HTML content by inserting a table of contents (TOC) 
#' at a specified placeholder and writes the final HTML content to an output file.
#'
#' @param toc The HTML string representing the table of contents.
#' @param html_content The initial HTML content containing a placeholder for the TOC.
#' @param output_file_path The file path where the final HTML content will be written.
#'
#' @return NULL
#'
generate_and_write_html <- function(
    toc,
    html_content,
    output_file_path
    ) {
  
  output_file_path <- normalizePath(output_file_path, mustWork = FALSE)
  
  # Close the Table of Contents
  toc <- paste(toc, "</ul></div>", sep = "\n")
  
  # Insert the Table of Contents at the placeholder
  html_content <- gsub("<!--TOC-->", toc, html_content)
  
  # Append a horizontal line after the TOC
  html_content <- gsub("</ul></div>", "</ul></div>\n<hr>", html_content)
  
  # Append the final closing tags for the HTML body and document
  html_content <- paste(html_content, "</body></html>", sep = "\n")
  
  # Ensure the directory exists
  dir_path <- dirname(output_file_path)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  writeLines(html_content, output_file_path)
}



# Level 2 internal functions ---------------------------------------------------


#' Generate Column Headers String from Feature Name
#'
#' @description
#' This function takes a dataframe as input and returns a string that replaces the values
#' in the `feature_name` column of the first row with their respective column headers.
#'
#' @param df A dataframe containing a column called `feature_name`.
#'
#' @return A string in the format of the value in `feature_name` but containing
#' the column headers instead of the values.
#'
generate_feature_name_template <- function(df) {

  # Ensure the dataframe contains the feature_name column
  if (!"feature_name" %in% colnames(df)) {
    stop(
      "The dataframe must contain a column named 'feature_name'",
      call. = FALSE
    )
  }
  
  # Extract the feature_name from the first row
  feature_name_value <- df$feature_name[1]
  
  # Split the feature_name value into its components
  feature_name_parts <- strsplit(feature_name_value, "_")[[1]]
  
  # Initialize an empty vector to store the column headers
  column_headers <- vector()
  
  # Function to find the matching column header
  find_matching_column <- function(part) {
    for (i in 1:length(part)) {
      substring <- paste(part[1:i], collapse = "_")
      matching_columns <- which(sapply(df[1, ],
                                       function(x) identical(x, substring)) 
                                & colnames(df) != "feature_name")
      if (length(matching_columns) > 0) {
        return(colnames(df)[matching_columns])
      }
    }
    return(NULL)
  }
  
  # Iterate over the parts and find corresponding column headers
  i <- 1
  while (i <= length(feature_name_parts)) {
    j <- i
    while (j <= length(feature_name_parts)) {
      substring <- paste(feature_name_parts[i:j], collapse = "_")
      matching_column <- find_matching_column(feature_name_parts[i:j])
      if (!is.null(matching_column)) {
        column_headers <- c(column_headers, matching_column)
        i <- j + 1
        break
      }
      j <- j + 1
    }
    if (i == j) {
      column_headers <- c(column_headers, feature_name_parts[i])
      i <- i + 1
    }
  }
  
  # Combine the column headers into a single string with underscores
  result_string <- paste(column_headers, collapse = "_")
}


#' Format text
#'
#' @description
#' This function takes a character vector `text` and splits it into individual
#' characters. It then iterates over the characters and builds lines not 
#' exceeding
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
#' @param report_type A character specifying the type of HTML report.
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
get_header_section <- function(
    title,
    header_text,
    report_type,
    feature_names_formula
    ) {
  
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
    "  margin-top: 0; margin-bottom: 0;",  # Reset margin to ensure no extra space
    "}",
    ".logo {",
    "  position: absolute;",  # Position the logo absolutely
    "  top: 60px;",           # Adjust the top position to move the logo down
    "  right: 0;",            # Align with the right of the h1
    "  width: 400px;",        # Adjust the width to make the logo smaller
    "  height: auto;",        # Maintain aspect ratio
    "}",
    "hr {",
    "  margin-top: 20px;",    # Space above the line
    "  margin-bottom: 20px;", # Space below the line
    "}",
    "</style>",
    "</head><body>",
    "<h1>", header_text, "<img src='",
    logo_base64, "' alt='Logo' class='logo'></h1>",
    "<table>",
    sep = ""
  )
  
  sentence <- switch(
    report_type,
    "explore_data" = paste(
      '<div style="border: 2px solid #f00; padding: 15px;', 
      'position: relative; margin-bottom: 20px;', 
      'background-color: #fee; font-family: Arial,', 
      'sans-serif; width: 65%;">',
      '<div style="position: absolute; top: -25px;', 
      'right: -25px; transform: rotate(45deg);', 
      'background-color: #f00; color: #fff;', 
      'padding: 10px 15px; font-size: 2em;', 
      'font-weight: bold; z-index: 1;">Note!</div>',
      '<p style="font-size: 2em;">',
      "This HTML report contains the exploratory",
      "data analysis plots, such as density and", 
      "PCA plots. <br> Right-click on ", 
      "any plot in this report to save it as a .svg (vector graphic) file!</p>",
      '</div>'
      ),
    "screen_limma_hyperparams" = '<p style="font-size: 2em;"></p>',
    "create_limma_report" = paste(
      '<div style="border: 2px solid #f00; padding: 15px; position: relative;', 
      'margin-bottom: 20px; background-color: #fee; font-family: Arial,', 
      'sans-serif; width: 65%;">',
      '<div style="position: absolute; top: -25px; right: -25px; transform:', 
      'rotate(45deg); background-color: #f00; color: #fff; padding: 10px 15px;', 
      'font-size: 2em; font-weight: bold; z-index: 1;">Note!</div>',
      '<p style="font-size: 2em;">',
      "This HTML report contains all the plots visualizing the results from", 
      "the limma topTables. <br> Right-click on", 
      "any plot in this report to save it as a .svg (vector graphic) file!</p>",
      '</div>'
      ),
    "cluster_hits" = paste(
      '<div style="border: 2px solid #f00; padding: 15px;', 
      'position: relative; margin-bottom: 20px;', 
      'background-color: #fee; font-family: Arial, sans-serif; width: 65%;">',
      '<div style="position: absolute; top: -20px;', 
      'right: -27px; transform: rotate(45deg);', 
      'background-color: #f00; color: #fff; padding:', 
      '10px 15px; font-size: 2em; font-weight: bold;', 
      'z-index: 1;">Note!</div>',
      '<p style="font-size: 2em;">',
      "Clustering of features that show", 
      "significant changes over time (= hits).<br>", 
      "Clustering was done based on the min-max normalized shape",
      "of the spline. <br> Right-click on any plot", 
      "in this report to save it as a", 
      ".svg (vector graphic) file! <br> If one level of ", 
      "the experiment is not shown, it means it has < 2 hits!</p>",
      paste(
        '<span style="font-size:1.3em;">',
        'feature_name "formula": ', 
        '{annotation-column-x}_{annotation-column-y}_ ... :', 
        '<br><b>', 
        feature_names_formula, 
        '</b></span>'
      ),
      '</div>'
      ),
    "create_gsea_report" = '<p style="font-size: 2em;"></p>')
  
  header_section <- paste(
    header_section,
    "<p>", sentence, "</p>",
    "</body></html>",
    sep = ""
    )
  
  return(header_section)
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
  mime_type <- 
    "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
  
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
#' @seealso
#'   \link[base64enc]{dataURI}
#' 
#' @importFrom ggplot2 ggsave
#' @importFrom base64enc dataURI
#' @importFrom svglite svglite
#' 
plot2base64 <- function(
    plot, 
    height, 
    width = 7, 
    base_height_per_row = 2.5, 
    units = "in", 
    html_img_width = "100%"
    ) {
  
  additional_height_per_row <- 2.1
  height <- base_height_per_row + (height - 1) * additional_height_per_row
  
  # Create a temporary file for the SVG. SVG does not specify the quality 
  # already, but later, after exporting the figures from the HTML, you can
  # specify the quality.
  img_file <- tempfile(fileext = ".svg")

  svglite::svglite(file = img_file, width = width, height = height)
  
  # Draw the plot
  print(plot)

  # Turn off the device
  dev.off()
  
  # Read the SVG file content
  svg_content <- readLines(img_file, warn = FALSE)
  
  # Convert the SVG content to a single string
  svg_string <- paste(svg_content, collapse = "\n")
  
  # Encode the SVG content as base64
  svg_base64 <- base64enc::dataURI(charToRaw(svg_string),
                                   mime = "image/svg+xml")
  
  # Delete the temporary SVG file
  unlink(img_file)
  
  # Return the HTML img tag with the base64 string and a fixed width
  return(
    sprintf(
      '<img src="%s" alt="Plot" style="width:%s;">', 
      svg_base64, html_img_width
      )
    )
}



#' Create Table of Contents
#'
#' @description
#' Creates the HTML content for the Table of Contents.
#'
#' @return A string containing the HTML for the Table of Contents.
#' 
create_toc <- function() {
  
  toc <- paste("<hr style='border: none; height: 3px; background-color:
               #333; margin: 40px 0;'>",
               "<div id='toc' style='text-align: center; display: block; margin:
               auto; width: 80%;'>",
               "<h2 style='font-size: 60px;'>Table of Contents</h2>",
               "<ul style='display: inline-block; text-align: left;'>",
               sep = ""
               )
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
  toc_style <- "font-size: 40px;"
  
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
# process_plots <- function(
#     plots_element,
#     plots_size,
#     html_content
#     ) {
#   
#   
#   img_tag <- plot2base64(
#     plots_element,
#     plots_size
#     )
#   html_content <- paste(
#     html_content,
#     img_tag,
#     sep = "\n"
#     )
# }
process_plots <- function(
    plots_element,
    element_name,
    plots_size,
    html_content
    ) {
  
  if (startsWith(element_name, "individual_spline_plots")) {
    # Convert each individual plot to base64 and store in a list
    base64_list <- lapply(plots_element, function(plot) {
      # Extract the title from the plot
      title <- plot$labels$title
      plot <- plot + ggplot2::labs(title = NULL)  # Remove the title from the plot
      
      # Convert the plot to base64
      base64_plot <- plot2base64(plot, height = plots_size)
      
      # Return a list containing the title and the base64 plot
      list(title = title, plot = base64_plot)
    })
    
    # Arrange the base64 images and titles in a single-column layout
    grid_content <- ""
    for (i in seq_along(base64_list)) {
      # Add the title as HTML text
      grid_content <- paste0(grid_content, '<div style="padding: 5px; text-align: center; font-size: 32px;">', base64_list[[i]]$title, '</div>')
      # Start a new row for each plot
      grid_content <- paste0(grid_content, '<div style="display: flex;">')
      grid_content <- paste0(grid_content, '<div style="flex: 1; padding: 5px;">', base64_list[[i]]$plot, '</div>')
      grid_content <- paste0(grid_content, '</div>')
    }
    
    # Add the grid content to the HTML content
    html_content <- paste(
      html_content,
      grid_content,
      sep = "\n"
    )
  }
  else {
    # Process a single plot
    img_tag <- plot2base64(plots_element, height = plots_size)
    html_content <- paste(
      html_content,
      img_tag,
      sep = "\n"
    )
  }
  
  return(html_content)
}


#' Process and Encode Data Field for Report
#'
#' @description
#' This function processes a given field, encodes the associated data as base64, 
#' and generates a download link for the report. It handles different types of 
#' fields including data, meta, top tables, and Enrichr formatted gene lists.
#'
#' @param field A string specifying the field to process.
#' @param data A dataframe containing the main data.
#' @param meta A dataframe containing meta information.
#' @param topTables A dataframe containing the results of differential 
#'                  expression analysis.
#' @param report_info A list containing additional report information.
#' @param encode_df_to_base64 A function to encode a dataframe to base64.
#' @param report_type A string specifying the type of report.
#' @param enrichr_format A list with the formatted gene lists and background 
#'                       gene list.
#'
#' @return A string containing the HTML link for downloading the processed 
#'         field.
#'
#' @importFrom base64enc base64encode
#' 
process_field <- function(
    field,
    data,
    meta,
    topTables,
    report_info,
    encode_df_to_base64,
    report_type,
    enrichr_format
    ) {
  
  if (field == "data_with_annotation")  {
    base64_df <- sprintf('<a href="%s" download="data.xlsx">
                         <button>Download data_with_annotation.xlsx
                         </button></a>', 
                         encode_df_to_base64(data))
    
  } else if (field == "meta" &&
             !is.null(meta) &&
             is.data.frame(meta) &&
             !any(is.na(meta))) {
    base64_df <- sprintf('<a href="%s" download="meta.xlsx">
                         <button>Download meta.xlsx</button></a>', 
                         encode_df_to_base64(meta))

  } else if (field == "limma_topTables" && !any(is.na(topTables)))  {
    base64_df <- sprintf('<a href="%s" download="limma_topTables.xlsx">
                         <button>Download limma_topTables.xlsx</button></a>', 
                         encode_df_to_base64(topTables))

  } else if (field == "Enrichr_clustered_genes" && 
             !any(is.na(enrichr_format)) &&
             !is.null(enrichr_format$gene_lists)) {
    
    # Create ZIP file for Enrichr_clustered_genes
    zip_base64 <- create_enrichr_zip(enrichr_format)
    base64_df <- sprintf(
      paste('<a href="data:application/zip;base64,%s" download=',
            '"Enrichr_clustered_genes.zip">',
            '<button>Download Enrichr_clustered_genes.zip</button></a>'),
      zip_base64
    )
    
  } else if (field == "Enrichr_background" && 
             !any(is.na(enrichr_format)) &&
             !is.null(enrichr_format$background)) {
    
    base64_df <- sprintf(
      paste(
        '<a href="data:text/plain;base64,%s"', 
        'download="Enrichr_background.txt">',
        '<button>Download Enrichr_background.txt</button></a>'
      ),
      base64enc::base64encode(charToRaw(enrichr_format$background))
    )
    
  } else {
    base64_df <- ifelse(is.null(report_info[[field]]),
                        "NA", report_info[[field]])
  }
  return(base64_df)
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


#' Create a ZIP File for Enrichr Gene Lists
#'
#' @description
#' This function creates a ZIP file containing directories for each level of 
#' gene lists. Each directory contains text files for each cluster. The ZIP file 
#' is then encoded to base64 for easy download.
#'
#' @param enrichr_format A list with the formatted gene lists and background 
#' gene list, typically the output of `prepare_gene_lists_for_enrichr`.
#'
#' @return A base64-encoded string representing the ZIP file.
#'
#' @details
#' The function creates a temporary directory to store the files. For each level 
#' in the `enrichr_format$gene_lists`, it creates a directory named after the 
#' level. Within each level directory, it creates a text file for each cluster, 
#' containing the genes in that cluster. The directories and files are added 
#' to a ZIP file, which is then encoded to base64.
#'
#' @importFrom zip zip
#' @importFrom base64enc base64encode
#' 
create_enrichr_zip <- function(enrichr_format) {
  
  temp_dir <- tempfile(pattern = "enrichr")
  dir.create(temp_dir)
  zip_file <- tempfile(fileext = ".zip")
  
  for (level in names(enrichr_format$gene_lists)) {
    
    level_dir <- file.path(temp_dir, level)
    dir.create(level_dir, recursive = TRUE)
    
    for (cluster in names(enrichr_format$gene_lists[[level]])) {
      
      cluster_file <- file.path(level_dir, paste0(cluster, ".txt"))
      writeLines(enrichr_format$gene_lists[[level]][[cluster]], cluster_file)
    }
  }
  
  # Create the ZIP file using relative paths
  original_wd <- getwd()
  setwd(temp_dir)
  files_to_zip <- list.files(temp_dir, recursive = TRUE)
  zip::zip(zipfile = zip_file, files = files_to_zip)
  setwd(original_wd)
  
  # Read the ZIP file and encode it to base64
  zip_base64 <- base64enc::base64encode(zip_file)
  
  # Clean up temporary directory and files
  unlink(temp_dir, recursive = TRUE)
  unlink(zip_file)
  
  return(zip_base64)
}



