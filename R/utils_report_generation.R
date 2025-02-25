#' utils scripts contains shared functions that are used by at least two package
#' functions of the SplineOmics package. The level separation is only valid
#' internally in this script, and has no connection to the script level of the
#' respective exported functions scripts.

# Level 1 internal functions ---------------------------------------------------


#' Generate Report HTML
#' 
#' @noRd
#'
#' @description
#' Generates an HTML report with the provided plots, spline parameters, and
#' report information.
#'
#' @param plots A list of ggplot2 plot objects.
#' @param plots_sizes A list of integers specifying the size of each plot.
#' @param report_info A named list containing report information.
#' @param limma_result_2_and_3_plots List containing the list of lists with all
#' the plots for all the pairwise comparisons of the condition in terms of
#' average spline diff and interaction condition time, and another list of lists
#' where the respective names of each plot are stored.
#' @param data A dataframe or a list of dataframes, containing data that should
#'             be directly embedded in the HTML report for downloading.
#' @param meta A dataframe, containing metadata that should
#'             be directly embedded in the HTML report for downloading.
#' @param topTables List of limma topTables.
#' @param category_2_and_3_hits List of dataframes, where each df is the part
#' of the toptable that contains the significant features of the respective 
#' limma result category (2 or 3).
#' @param enrichr_format List, containing two lists: The gene list and the list
#'                       of background genes.
#' @param level_headers_info A list of header information for each level.
#' @param spline_params A list of spline parameters, such as dof and type.
#' @param adj_pthresholds Numeric vector with the values for the adj.p.tresholds
#'                        for each level.
#' @param adj_pthresh_avrg_diff_conditions Float, only for cluster_hits()
#' @param adj_pthresh_interaction_condition_time Float, only for cluster_hits()
#' @param report_type A character string specifying the report type
#'                    ('screen_limma_hyperparams' or 'cluster_hits').
#' @param feature_name_columns Character vector with the column names of the
#'                             annotation information, such as the columns
#'                             containing the gene names. These column names
#'                             are used to put the info in the HTML reports on
#'                             how the descriptions above the individual spline
#'                             plots where created. This is because those
#'                             descriptions can be made up of several column
#'                             values, and the specific columns are then stated
#'                             in the HTML report on top (e.g gene_uniprotID).
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
#'
generate_report_html <- function(
    plots,
    plots_sizes,
    report_info,
    limma_result_2_and_3_plots = NULL, # only for build_cluster_hits_report
    data = NULL,
    meta = NA,
    topTables = NA,
    category_2_and_3_hits = NA,  # only for cluster_hits()
    enrichr_format = NA,
    level_headers_info = NA,
    spline_params = NA,
    adj_pthresholds = NA,
    adj_pthresh_avrg_diff_conditions = NA, # only for cluster_hits()
    adj_pthresh_interaction_condition_time = NA, # only for cluster_hits()
    report_type = "explore_data",
    feature_name_columns = NA, # only for cluster_hits()
    mode = NA,
    filename = "report",
    timestamp = format(
      Sys.time(),
      "%d_%m_%Y-%H_%M_%S"
    ),
    report_dir = here::here()
    ) {
  
  feature_names_formula <- ""

  if (report_type == "explore_data") {
    if (filename == "explore_data") {
      title <- "explore data"
    } else {
      title <- "explore batch-corrected data"
    }
  } else if (report_type == "screen_limma_hyperparams") {
    title <- paste("hyperparams screen |", filename)
  } else if (report_type == "create_limma_report") {
    title <- "l(m)m report"
  } else if (report_type == "cluster_hits") {
    title <- "clustered hits"
    feature_names_formula <- paste(
      feature_name_columns,
      collapse = "_"
    )
  } else if (report_type == "create_gsea_report") {
    title <- "gsea report"
  } else {
    stop(
      paste(
        "report_type must be explore_hits, screen_limma_hyperparams,",
        "create_limma_report, or cluster_hits"
      ),
      call. = FALSE
    )
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

  splineomics_version <- utils::packageVersion("SplineOmics")
  header_text <- paste(
    title,
    paste("Omics-Datatype:", report_info$omics_data_type),
    paste("Date-Time:", timestamp),
    paste("SplineOmics Version:", splineomics_version),
    sep = " | "
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
    "plot_data_batch_correction",
    "Fixed effects",
    "Random effects",
    "mode",
    "analyst_name",
    "contact_info",
    "project_name",
    "method_description",
    "results_summary",
    "conclusions"
  )

  download_fields <- c()
  if (!is.null(data)) {
    download_fields <- c(
      download_fields,
      "data_with_annotation"
    )
  }

  if (!all(is.na(meta))) {
    download_fields <- c(
      download_fields,
      "meta"
    )
  }

  if (!all(is.na(topTables))) {
    download_fields <- c(
      download_fields,
      if (report_type == "cluster_hits") {
        "limma_topTables_clustered_time_effect_hits"
      } else {
        "limma_topTables"
      }
    )
  }
  
  if (!all(is.na(category_2_and_3_hits))) {
    download_fields <- c(
      download_fields,
      "limma_topTables_avrg_diff_conditions_hits",
      "limma_topTables_interaction_condition_time_hits"
    )
  }

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
    "background-color: #333; margin: 40px 0;",
    'width: 75%;"> <h2 style="font-size: 48px;">',
    "Report Info \u2139</h2><table>",
    sep = ""
  )

  downloads_section <- paste(
    '<hr><h2 style="font-size: 48px;">',
    "Downloads \U0001F4E5</h2><table>"
  )
  
  # Because these are not part of report_info already.
  report_info[["mode"]] <- mode 

  for (field in report_info_fields) {
    base64_df <- process_field(
      field = field,
      data = data,
      meta = meta,
      topTables = topTables,
      report_info = report_info,
      encode_df_to_base64 = encode_df_to_base64,
      enrichr_format = enrichr_format
    )

    field_display <- sprintf(
      "%-*s",
      max_field_length,
      gsub("_", " ", field)
    )

    report_info_section <- paste(
      report_info_section,
      sprintf(
        '<tr><td style="text-align: right;
           color:blue; padding-right: 5px;">%s
           :</td><td>%s</td></tr>',
        field_display, base64_df
      ),
      sep = "\n"
    )
  }

  # Close the Report Info table
  report_info_section <- paste(report_info_section, "</table>", sep = "\n")

  download_fields <- c(download_fields, "session_info")
  download_fields <- c(download_fields, "analysis_script")
  
  category_2_and_3_hit_counts <- count_hits(category_2_and_3_hits)
  
  # Truncate all df names so that they are <= 30 chars (required for sheet 
  # names of Excel files) --> every df becomes a sheet.
  category_2_and_3_hits <- truncate_nested_list_names(category_2_and_3_hits)
  
  for (field in download_fields) {
    base64_df <- process_field(
      field = field,
      data = data,
      meta = meta,
      topTables = topTables,
      category_2_and_3_hits = category_2_and_3_hits,
      report_info = report_info,
      encode_df_to_base64 = encode_df_to_base64,
      enrichr_format = enrichr_format
    )

    field_display <- sprintf("%-*s", max_field_length, gsub("_", " ", field))
    downloads_section <- paste(
      downloads_section,
      sprintf(
        '<tr><td style="text-align: right;
           color:blue; padding-right: 5px;">%s
           :</td><td>%s</td></tr>',
        field_display, base64_df
      ),
      sep = "\n"
    )
  }

  # Close the Downloads table
  downloads_section <- paste(downloads_section, "</table>", sep = "\n")

  # Preserve initial header_section content
  header_section <- paste(
    header_section,
    report_info_section,
    downloads_section,
    sep = "\n"
  )


  if (report_type == "create_gsea_report") {
    databases_text <- paste(report_info$databases, collapse = ", ")
    header_section <- paste(
      header_section,
      "<p style='font-size: 20px;'>Databases used: ",
      databases_text,
      "</p>",
      sep = "\n"
    )
  }

  file_name <- sprintf(
    "%s_%s_%s.html",
    filename,
    report_info$omics_data_type,
    timestamp
  )

  output_file_path <- here::here(report_dir, file_name)

  if (report_type == "explore_data") {
    build_explore_data_report(
      header_section = header_section,
      plots = plots,
      plots_sizes = plots_sizes,
      report_info = report_info,
      output_file_path = output_file_path
    )
  } else if (report_type == "screen_limma_hyperparams") {
    build_hyperparams_screen_report(
      header_section = header_section,
      plots = plots,
      plots_sizes = plots_sizes,
      report_info = report_info,
      output_file_path = output_file_path
    )
  } else if (report_type == "create_limma_report") {
    build_create_limma_report(
      header_section = header_section,
      plots = plots,
      plots_sizes = plots_sizes,
      level_headers_info = level_headers_info,
      report_info = report_info,
      output_file_path = output_file_path
    )
  } else if (report_type == "create_gsea_report") {
    build_create_gsea_report(
      header_section = header_section,
      plots = plots,
      plots_sizes = plots_sizes,
      level_headers_info = level_headers_info,
      report_info = report_info,
      output_file_path = output_file_path
    )
  } else { # report_type == "cluster_hits"
    build_cluster_hits_report(
      header_section = header_section,
      plots = plots,
      limma_result_2_and_3_plots = limma_result_2_and_3_plots,
      plots_sizes = plots_sizes,
      level_headers_info = level_headers_info,
      spline_params = spline_params,
      adj_pthresholds = adj_pthresholds,
      adj_pthresh_avrg_diff_conditions = adj_pthresh_avrg_diff_conditions,
      adj_pthresh_interaction_condition_time =
        adj_pthresh_interaction_condition_time,
      category_2_and_3_hit_counts = category_2_and_3_hit_counts,
      mode = mode,
      report_info = report_info,
      output_file_path = output_file_path
    )
  }
}


#' Generate and Write HTML Report
#' 
#' @noRd
#'
#' @description
#' This function generates an HTML report by inserting a table of contents,
#' embedding necessary JavaScript files, and writing the final HTML content
#' to a specified output file.
#'
#' @param toc A string containing the table of contents in HTML format.
#' @param html_content A string containing the main HTML content with a
#' placeholder for the table of contents.
#' @param report_info A list containing report information such as
#' `contact_info` and `analyst_name`.
#' @param output_file_path A string specifying the path where the final
#' HTML file will be written.
#'
generate_and_write_html <- function(
    toc,
    html_content,
    report_info,
    output_file_path) {
  output_file_path <- normalizePath(
    output_file_path,
    mustWork = FALSE
  )

  # Close the Table of Contents
  toc <- paste(
    toc,
    "</ul></div>",
    sep = "\n"
  )

  # Insert the Table of Contents at the placeholder
  html_content <- gsub(
    "<!--TOC-->",
    toc,
    html_content
  )

  # Append a horizontal line after the TOC
  html_content <- gsub(
    "</ul></div>",
    "</ul></div>\n<hr>",
    html_content
  )

  # Path to the external JavaScript file within the package
  js_file_path <- normalizePath(
    system.file(
      "www/hotkeys.js",
      package = "SplineOmics"
    ),
    mustWork = FALSE
  )
  if (js_file_path == "") {
    stop("JavaScript file not found.")
  }

  # Read the JavaScript file and replace placeholders with actual values
  js_content <- readLines(
    js_file_path,
    encoding = "UTF-8"
  )
  js_content <- gsub(
    "\\{\\{email\\}\\}",
    report_info$contact_info,
    js_content
  )
  js_content <- gsub(
    "\\{\\{name\\}\\}",
    report_info$analyst_name,
    js_content
  )

  # Read the content of JSZip and FileSaver JavaScript files as text
  jszip_path <- normalizePath(
    system.file(
      "www/jszip.min.js",
      package = "SplineOmics"
    ),
    mustWork = FALSE
  )
  filesaver_path <- normalizePath(
    system.file(
      "www/FileSaver.min.js",
      package = "SplineOmics"
    ),
    mustWork = FALSE
  )

  if (!file.exists(jszip_path)) {
    stop("JSZip file not found at: ", jszip_path)
  }
  if (!file.exists(filesaver_path)) {
    stop("FileSaver.js file not found at: ", filesaver_path)
  }

  jszip_content <- readLines(
    jszip_path,
    encoding = "UTF-8",
    warn = FALSE
  )
  filesaver_content <- readLines(
    filesaver_path,
    encoding = "UTF-8",
    warn = FALSE
  )

  # Combine all JavaScript content
  combined_js_content <- c(
    "<script>",
    jszip_content,
    filesaver_content,
    js_content,
    "</script>"
  )

  # Properly escape special characters in JavaScript content
  combined_js_content <- paste(
    combined_js_content,
    collapse = "\n"
  )
  combined_js_content <- gsub(
    "\\\\",
    "\\\\\\\\",
    combined_js_content
  ) # Escape backslashes
  combined_js_content <- gsub( # Escape double quotes
    "\"",
    "\\\"",
    combined_js_content
  )

  # Embed the combined JavaScript content before the closing body tag
  script_tag <- paste(
    combined_js_content,
    collapse = "\n"
  )
  html_content <- gsub(
    "</body>",
    paste(
      script_tag,
      "</body>",
      sep = "\n"
    ),
    html_content
  )

  # Append the final closing tags for the HTML body and document
  html_content <- paste(
    html_content,
    "</body></html>",
    sep = "\n"
  )

  # Ensure the directory exists
  dir_path <- dirname(output_file_path)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }

  con <- file(
    output_file_path,
    "w",
    encoding = "UTF-8"
  )
  writeLines(
    html_content,
    con,
    useBytes = TRUE
  )
  close(con)
}


#' Read and split section texts from a file
#' 
#' @noRd
#'
#' @description
#' This internal function reads the contents of a text file located in the
#' `inst/descriptions` directory of the package and splits it into individual
#' sections based on a specified delimiter.
#'
#' @param filename A character string specifying the name of the file
#' containing the section texts. The file should be located in the
#' `inst/descriptions` directory of the package.
#'
#' @return A character vector where each element is a section of the text
#' split by the delimiter `|`.
#'
read_section_texts <- function(filename) {
  file_path <- system.file(
    "descriptions",
    filename,
    package = "SplineOmics"
  )
  content <- readLines(
    file_path,
    warn = FALSE
  ) |> paste(collapse = " ")
  # Split the content by the delimiter
  strsplit(content, "\\|")[[1]]
}


# Level 2 internal functions ---------------------------------------------------


#' Format text
#' 
#' @noRd
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
#' @noRd
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
#' @param feature_names_formula String describing which columns of the
#'                              annotation info, such as gene and uniprotID,
#'                              where used to construct the description above
#'                              the individual spline plots. This is placed in
#'                              the beginning of the output HTML reports.
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
    feature_names_formula) {
  if (feature_names_formula == "") {
    feature_names_formula <- "No feature name columns provided!"
  }

  if (Sys.getenv("DEVTOOLS_LOAD") == "true") {
    logo_path <- file.path(
      "inst",
      "logos",
      "SplineOmics_logo.png"
    )
  } else {
    logo_path <- system.file(
      "logos",
      "SplineOmics_logo.png",
      package = "SplineOmics"
    )
  }

  logo_base64 <- base64enc::dataURI(
    file = logo_path,
    mime = "image/png"
    )

  header_section <- paste(
    "<html><head><title>", title, "</title>",
    "<meta charset=\"UTF-8\">", # Ensure UTF-8 encoding (JavaScript issues)
    "<style>",
    
    # Existing Table Styles
    "table { font-size: 30px; }",
    "td { padding: 8px; text-align: left; }",
    "td:first-child { text-align: right; color: blue; }",
    
    # Header and Logo Styles
    "h1 {",
    "  color: #333333;",
    "  display: flex;",
    "  align-items: center;",
    "  justify-content: space-between;",
    "  margin-top: 0;",
    "  margin-bottom: 0;",
    "}",
    ".logo {",
    "  position: absolute;",
    "  top: 60px;",
    "  right: 0;",
    "  width: 400px;",
    "  height: auto;",
    "}",
    
    # Horizontal Rule
    "hr {",
    "  margin-top: 20px;",
    "  margin-bottom: 20px;",
    "}",
    
    # Zoom Container
    ".zoom-container {",
    "  display: block;",
    "  text-align: center;",
    "  margin: 10px auto;",
    "  position: relative;",
    "}",
    
    # Zoomable Image
    ".zoomable-image {",
    "  width: 100%;",
    "  max-width: 800px;",
    "  cursor: zoom-in;",
    "  transition: transform 0.3s ease;",
    "}",
    
    # Zoom Button
    ".zoom-button {",
    "  position: absolute;",
    "  top: 10px;",
    "  right: 10px;",
    "  background: rgba(0, 0, 0, 0.6);",
    "  color: white;",
    "  border: none;",
    "  padding: 8px 12px;",
    "  border-radius: 50%;",
    "  cursor: pointer;",
    "  font-size: 18px;",
    "}",
    
    # Zoom Button Hover Effect
    ".zoom-button:hover {",
    "  background: rgba(0, 0, 0, 0.8);",
    "}",

    ".zoom-container.zoom-active .zoomable-image { cursor: zoom-out; }",
    
    "</style>",
    
    "<script>",
    "document.addEventListener('DOMContentLoaded', function() {",
    "  console.log('JavaScript is running!');",
    "  document.querySelectorAll('.zoom-container').forEach(container => {",
    "    console.log('Zoom container found');",
    "    const img = container.querySelector('.zoomable-image');",
    "    let scale = 1, zoomActive = false, isDragging = false;",
    "    let startX = 0, startY = 0, translateX = 0, translateY = 0;",
    "    let imgOffsetX = 0, imgOffsetY = 0;",
    "",
    "    img.addEventListener('click', function(event) {",
    "      zoomActive = !zoomActive;",
    "      console.log('Zoom toggled:', zoomActive);",
    "",
    "      if (!zoomActive) {",
    "        scale = 1; imgOffsetX = 0; imgOffsetY = 0;", 
    "        img.style.transform = 'scale(1) translate(0px, 0px)';",  
    "        img.style.cursor = 'zoom-in'; img.style.zIndex = '1';",
    "      } else {",
    "        img.style.cursor = 'grab'; img.style.zIndex = '1000';",
    "        img.style.position = 'relative';",
    "      }",
    "    });",
    "",
    "    img.addEventListener('wheel', function(event) {",
    "      if (!zoomActive) return;",
    "      event.preventDefault();",
    "      let rect = img.getBoundingClientRect();",
    "      let mouseX = (event.clientX - rect.left) / rect.width;",
    "      let mouseY = (event.clientY - rect.top) / rect.height;",
    "",
    "      scale += event.deltaY * -0.01;",  
    "      scale = Math.min(Math.max(1, scale), 5);",  
    "      imgOffsetX = (0.5 - mouseX) * 100 * (scale - 1);", 
    "      imgOffsetY = (0.5 - mouseY) * 100 * (scale - 1);",
    "",
    "      img.style.transform = `scale(${scale}) translate(${imgOffsetX}px,", 
    "        ${imgOffsetY}px)`;",
    "      img.style.transformOrigin = 'center center';",
    "      console.log('Zooming:', scale, 'Offset:', imgOffsetX, imgOffsetY);",
    "    });",
    "",
    "    img.addEventListener('mousedown', function(event) {",
    "      if (!zoomActive || scale === 1) return;",
    "      isDragging = true;",
    "      startX = event.clientX; startY = event.clientY;",
    "      img.style.cursor = 'grabbing';",
    "    });",
    "",
    "    document.addEventListener('mousemove', function(event) {",
    "      if (!isDragging) return;",
    "      event.preventDefault();",
    "      let dx = (event.clientX - startX) / scale;",
    "      let dy = (event.clientY - startY) / scale;",
    "      imgOffsetX += dx; imgOffsetY += dy;",
    "",
    "      img.style.transform = `scale(${scale}) translate(${imgOffsetX}px,", 
    "        ${imgOffsetY}px)`;",
    "      startX = event.clientX; startY = event.clientY;",
    "    });",
    "",
    "    document.addEventListener('mouseup', function() {",
    "      isDragging = false;",
    "      if (zoomActive) img.style.cursor = 'grab';",
    "    });",
    "  });",
    "});",
    "</script>",

    "</head><body>",
    
    "</style>",
    "</head><body>",
    
    # Add the Title with Logo
    "<h1>", 
    header_text,
    "<img src='",
    logo_base64,
    "' alt='Logo' class='logo'></h1>",
    "<table>",
    
    sep = ""
  )
  

  note <- switch(report_type,
    "explore_data" = paste(
      '<div style="border: 2px solid #f00; padding: 15px;',
      "position: relative; margin-bottom: 20px;",
      "background-color: #fee; font-family: Arial,",
      'sans-serif; width: 65%;">',
      '<div style="position: absolute; top: -25px;',
      "right: -25px; transform: rotate(45deg);",
      "background-color: #f00; color: #fff;",
      "padding: 10px 15px; font-size: 2em;",
      'font-weight: bold; z-index: 1;">Note!</div>',
      '<p style="font-size: 2em;">',
      "This HTML report contains the exploratory",
      "data analysis plots, (e.g. density plots) <br> Right-click on",
      "any plot in this report to save it as a .svg (vector graphic) file!</p>",
      "</div>"
    ),
    "screen_limma_hyperparams" = '<p style="font-size: 2em;"></p>',
    "create_limma_report" = paste(
      '<div style="border: 2px solid #f00; padding: 15px; position: relative;',
      "margin-bottom: 20px; background-color: #fee; font-family: Arial,",
      'sans-serif; width: 65%;">',
      '<div style="position: absolute; top: -25px; right: -25px; transform:',
      "rotate(45deg); background-color: #f00; color: #fff; padding: 10px 15px;",
      'font-size: 2em; font-weight: bold; z-index: 1;">Note!</div>',
      '<p style="font-size: 2em;">',
      "This HTML report contains plots visualizing the results from",
      "the limma topTables. <br> Right-click on",
      "any plot in this report to save it as a .svg (vector graphic) file!",
      "<br><br>To understand the three limma result categories shown in this ",
      'report, please <a href="',
      file.path(system.file(
        "descriptions",
        "limma_result_categories.pdf",
        package = "SplineOmics"
      )),
      '" download>download and review this PDF document</a><br><br>
      The grey shaded areas of the plots in this report cover the
      non-significant features!</p>',
      "</div>"
    ),
    "cluster_hits" = paste(
      '<div style="border: 2px solid #f00; padding: 15px;',
      "position: relative; margin-bottom: 20px;",
      'background-color: #fee; font-family: Arial, sans-serif; width: 65%;">',
      '<div style="position: absolute; top: -20px;',
      "right: -27px; transform: rotate(45deg);",
      "background-color: #f00; color: #fff; padding:",
      "10px 15px; font-size: 2em; font-weight: bold;",
      'z-index: 1;">Note!</div>',
      '<p style="font-size: 1em;">',
      '<ul style="font-size: 2em; padding-left: 20px;">',
      '<li style="margin-bottom: 15px;">Clustering of features that show
      significant changes over time (= hits).</li>',
      '<li style="margin-bottom: 15px;">Clustering was done based on the
      min-max normalized shape of the spline.</li>',
      '<li style="margin-bottom: 15px;">For this, the spline parameters in
      the limma topTable are used to generate 1000 curve datapoints between
      0 and 1.</li>',
      '<li style="margin-bottom: 15px;">These datapoints are used for
      hierarchical clustering.</li>',
      '<li style="margin-bottom: 15px;">Right-click on any plot in this
      report to save it as a .svg (vector graphic) file!</li>',
      '<li style="margin-bottom: 15px;">If one level of the experiment is
      not shown, it means it has < 2 hits!</li>',
      '<li style="margin-bottom: 15px;">The "avg CV" value on top of all
      individual spline plots is the average coefficient of variation across
      all timepoints. For example, a value of 10% means that the timepoints,
      on average, have a standard deviation of 10% of the mean.</li>',
      "</ul>",
      "</p>",
      "</div>",
      paste(
        '<span style="font-size:1.3em;">',
        'feature_name "formula": ',
        "{annotation-column-x}_{annotation-column-y}_ ... :",
        "<br><b>",
        feature_names_formula,
        "</b></span>"
      ),
      "</div>"
    ),
    "create_gsea_report" = '<p style="font-size: 2em;"></p>'
  )

  hotkeys_box <- paste(
    '<div style="border: 2px solid #00f; padding: 15px;',
    "position: relative; margin-bottom: 20px;",
    "background-color: #eef; font-family: Arial,",
    'sans-serif; width: 65%;">',
    '<div style="position: absolute; top: -5px;',
    "right: -65px; transform: rotate(45deg);",
    "background-color: #00f; color: #fff;",
    "padding: 10px 15px; font-size: 2em;",
    'font-weight: bold; z-index: 1;">Hotkeys</div>',
    '<p style="font-size: 2em;">',
    "Press:<br>",
    "<b>t</b>  --> Jump to <b>Table of Contents</b> and save current scroll",
    "position \U0001F4D1<br>",
    "<b>s</b> --> <b>Save</b> current scroll position \U0001F4CC<br>",
    "<b>b</b> --> Jump <b>back</b> to saved position \U0001F519<br>",
    "<b>d</b> --> <b>Download</b> all embedded files as zip \U0001F4E5<br>",
    "<b>e</b> --> Write an <b>email</b> to contact info \u2709<br>",
    "</p>",
    "</div>"
  )


  header_section <- paste(
    header_section,
    "<p>", note, "</p>",
    "<br><br>",
    hotkeys_box,
    "</body></html>",
    sep = ""
  )

  return(header_section)
}


#' Encode DataFrame to Base64 for HTML Embedding
#' 
#' @noRd
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
encode_df_to_base64 <- function(
    df,
    report_type = NA
    ) {
  
  temp_file <- tempfile(fileext = ".xlsx")
  wb <- openxlsx::createWorkbook()

  if (is.data.frame(df)) {
    # Convert single dataframe to Excel with one sheet
    sheet_name <- "Sheet1"
    openxlsx::addWorksheet(
      wb,
      sheet_name
    )
    openxlsx::writeData(
      wb,
      sheet = sheet_name,
      df
    )
  } else if (is.list(df) && all(vapply(df, is.data.frame, logical(1)))) {
    # Convert list of dataframes to Excel with multiple sheets

    if (!is.na(report_type)) {
      if (report_type == "create_gsea_report") {
        all_names <- names(df)
        sheet_names <- vapply(
          all_names,
          extract_and_combine,
          character(1)
        )
      }
    } else {
      sheet_names <- make.unique(names(df))
    }


    for (i in seq_along(sheet_names)) {
      openxlsx::addWorksheet(
        wb,
        sheet_names[i]
      )
      openxlsx::writeData(
        wb,
        sheet = sheet_names[i],
        df[[i]]
      )
    }
  } else {
    stop("Input must be a dataframe or a list of dataframes.")
  }

  openxlsx::saveWorkbook(
    wb,
    temp_file,
    overwrite = TRUE
  )

  # Read the file and encode to base64
  file_content <- readBin(
    temp_file,
    "raw",
    file.info(temp_file)$size
  )
  base64_file <- base64enc::base64encode(file_content)

  # Determine MIME type
  mime_type <-
    "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"

  # Create the data URI scheme
  data_uri <- paste0(
    "data:",
    mime_type,
    ";base64,",
    base64_file
  )

  # Remove the temporary file
  unlink(temp_file)

  return(data_uri)
}


#' Count Rows in Dataframes and Create Nested Dictionary
#'
#' @noRd
#'
#' @description
#' This function processes the nested list of `category_2_and_3_hits` to count
#' the rows in each dataframe. It organizes the counts into a nested list where:
#' - The outer keys are `category_2` and `category_3`.
#' - The inner keys are substrings of dataframe names after the second
#'   underscore.
#' - The values are the row counts of the corresponding dataframe.
#'
#' @param nested_list A nested list (e.g., category_2_and_3_hits) containing
#'   lists of named dataframes.
#'
#' @return A nested dictionary (list) with row counts.
#'
#' @examples
#' category_2_and_3_hits <- list(
#'   category_2_hits = list(
#'     "very_long_name_with_substring_ABC_vs_XYZ" = data.frame(A = 1:5),
#'     "another_long_name_DEF_vs_GHI" = data.frame(B = 6:10)
#'   ),
#'   category_3_hits = list(
#'     "extremely_long_name_with_substring_PQR_vs_STU" = data.frame(C = 11:15)
#'   )
#' )
#' row_counts_dict <- count_hits(category_2_and_3_hits)
count_hits <- function(nested_list) {
  # Initialize the output dictionary
  row_counts_dict <- list(
    category_2 = list(),
    category_3 = list()
  )
  
  # Helper function to extract substring after the second underscore
  extract_after_second_underscore <- function(name) {
    parts <- unlist(strsplit(name, "_"))
    if (length(parts) > 2) {
      paste(parts[3:length(parts)], collapse = "_")
    } else {
      name  # Return the full name if there aren't enough parts
    }
  }
  
  # Iterate over the outer list
  for (outer_key in names(nested_list)) {
    # Determine if we are processing category_2 or category_3
    category_key <- ifelse(
      outer_key == "category_2_hits", "category_2", "category_3"
      )
    inner_list <- nested_list[[outer_key]]
    
    # Iterate over the inner list of dataframes
    for (df_name in names(inner_list)) {
      # Get the substring after the second underscore
      extracted_name <- extract_after_second_underscore(df_name)
      
      # Count the rows in the dataframe
      row_count <- nrow(inner_list[[df_name]])
      
      # Store the result in the dictionary
      row_counts_dict[[category_key]][[extracted_name]] <- row_count
    }
  }
  
  return(row_counts_dict)
}


#' Truncate Dataframe Names in Nested Lists
#'
#' @noRd
#'
#' @description
#' This function processes a nested list (like category_2_and_3_hits) where each
#' element is a list of named dataframes. It ensures all dataframe names are
#' ≤ 30 characters, truncating them evenly for `_substring_vs_substring`.
#'
#' @param nested_list A list of lists containing named dataframes.
#'
#' @return The same nested list, but with dataframe names adjusted to meet the
#'         naming logic.
#'
#' @examples
#' category_2_and_3_hits <- list(
#'   category_2_hits = list(
#'     "very_long_name_with_substring_ABC_vs_XYZ" = data.frame(A = 1:5),
#'     "short_name" = data.frame(B = 6:10)
#'   ),
#'   category_3_hits = list(
#'     "extremely_long_name_with_substring_PQR_vs_STU" = data.frame(C = 11:15)
#'   )
#' )
#' result <- truncate_nested_list_names(category_2_and_3_hits)
truncate_nested_list_names <- function(nested_list) {
  
  # Helper function to truncate a single dataframe name
  truncate_name <- function(name) {
    max_length <- 30
    
    # Check if the name matches the special pattern
    if (grepl("_.*_vs_.*$", name)) {
      # Extract the base name and substrings
      base <- sub("(.*)_(.*)_vs_(.*)$", "\\1", name)
      substring1 <- sub("(.*)_(.*)_vs_(.*)$", "\\2", name)
      substring2 <- sub("(.*)_(.*)_vs_(.*)$", "\\3", name)
      
      # Calculate available length for substrings
      base_length <- nchar(base)
      remaining_length <- max_length - base_length - 5  # 5 for "_vs_"
      
      # Ensure there's space for both substrings
      if (remaining_length > 0) {
        trunc_length <- floor(remaining_length / 2)
        substring1 <- substr(substring1, 1, trunc_length)
        substring2 <- substr(substring2, 1, trunc_length)
        return(paste0(base, "_", substring1, "_vs_", substring2))
      } else {
        # If no space for substrings, truncate the entire name
        return(substr(name, 1, max_length))
      }
    } else {
      # For names that don't match the pattern, truncate normally
      return(substr(name, 1, max_length))
    }
  }
  
  # Process each inner list and truncate the names of its dataframes
  lapply(nested_list, function(inner_list) {
    # Update the names of the dataframes in the inner list
    names(inner_list) <- vapply(
      names(inner_list),
      truncate_name,
      FUN.VALUE = character(1)
      )
    inner_list
  })
}


#' Convert Plot to Base64
#' 
#' @noRd
#'
#' @description
#' Converts a ggplot2 plot to a Base64-encoded PNG image and returns an HTML
#' img tag for embedding in a report.
#'
#' @param plot A ggplot2 plot object.
#' @param height An integer specifying the height of the plot for correct
#'               representation in the HTML.
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
    html_img_width = "100%") {
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
  svg_base64 <- base64enc::dataURI(
    charToRaw(svg_string),
    mime = "image/svg+xml"
  )

  # Delete the temporary SVG file
  unlink(img_file)

  # # Return the HTML img tag with the base64 string and a fixed width
  # return(
  #   sprintf(
  #     '<img src="%s" alt="Plot" style="width:%s;">',
  #     svg_base64, html_img_width
  #   )
  # Return HTML with zoom functionality
  return(
    sprintf(
      '<div class="zoom-container">
         <img src="%s" alt="Plot" class="zoomable-image" style="width:%s;">
       </div>',
      svg_base64, html_img_width
    )
  )
}


#' Create Table of Contents
#' 
#' @noRd
#'
#' @description
#' Creates the HTML content for the Table of Contents.
#'
#' @return A string containing the HTML for the Table of Contents.
#'
create_toc <- function() {
  toc <- paste(
    "<hr style='border: none; height: 3px; background-color:
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
#' @noRd
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
#' @noRd
#'
#' @description
#' Converts plots to base64 and appends them to the HTML content.
#'
#' @param plots_element A list of plots to be processed.
#' @param element_name A character string specifying the name of the element.
#' @param plots_size A list of sizes for the plots.
#' @param html_content The current state of the HTML content.
#' @param toc The current state of the table of contents (TOC).
#' @param header_index An index to uniquely identify each section
#' for anchoring.
#'
#' @return Updated HTML content with the plots included.
#'
process_plots <- function(
    plots_element,
    plots_size,
    html_content,
    toc,
    header_index,
    element_name = NA) {
  if (
    !is.na(element_name) &&
      startsWith(
        element_name,
        "individual_spline_plots"
      )
  ) {
    spline_plots <- plots_element$spline_plots
    main_title <- plots_element$cluster_main_title

    # Create a TOC entry for the main title
    toc_entry <- paste0(
      "<li style='margin-left: 50px; font-size: 25px;'>",
      "<a href='#section", header_index, "'>",
      main_title,
      "</a></li>"
    )
    toc <- paste(
      toc,
      toc_entry,
      sep = "\n"
    )

    # Add the main title as a section title with an anchor before the first plot
    grid_content <- paste0(
      '<h2 id="section',
      header_index,
      '" style="text-align: center; margin-bottom: 20px; font-size: 50px;">',
      main_title, "</h2>"
    )

    # Convert each individual plot to base64 and store in a list
    base64_list <- lapply(spline_plots, function(plot) {
      title <- plot$labels$title
      
      # Remove the title from the plot (so that it is not there twice)
      plot <- plot + ggplot2::labs(title = NULL)
      
      # Remove HTML elements of the title, to add it back after removing
      plain_text_title <- gsub("<[^>]+>", "", title)
      plain_text_title <- sub("^\\s+", "", plain_text_title)
      plain_text_title <- sub(" ", "  | ", plain_text_title, fixed = TRUE)
      plot <- plot + ggplot2::labs(title = plain_text_title)
      

      base64_plot <- plot2base64(
        plot,
        height = plots_size
        )

      list(
        title = title,
        plot = base64_plot
      )
    })

    # Arrange the base64 images and titles in a single-column layout
    for (i in seq_along(base64_list)) {
      # Add the title as HTML text
      grid_content <- paste0(
        grid_content,
        '<div style="padding: 5px; text-align: center; font-size: 32px;">',
        base64_list[[i]]$title, "</div>"
      )
      # Start a new row for each plot
      grid_content <- paste0(
        grid_content,
        '<div style="display: flex;">'
      )
      grid_content <- paste0(
        grid_content,
        '<div style="flex: 1; padding: 5px;">', base64_list[[i]]$plot, "</div>"
      )

      grid_content <- paste0(grid_content, "</div>")

      # Add a horizontal line after each plot
      grid_content <- paste0(
        grid_content,
        '<hr style="border: 0; border-top: 1px solid #ccc; margin: 20px 0;">'
      )
    }

    # Add the grid content to the HTML content
    html_content <- paste(
      html_content,
      grid_content,
      sep = "\n"
    )
  } else if (
    !is.na(element_name) &&
      element_name == "cluster_mean_splines"
  ) {
    # One plot for each cluster
    for (i in seq_along(plots_element)) {
      html_content <- add_plot_to_html(
        html_content, plots_element[[i]], plots_size, header_index + i - 1
      )
    }
  } else {
    # Process a single plot
    html_content <- add_plot_to_html(
      html_content,
      plots_element,
      plots_size,
      header_index
    )
  }

  return(
    list(
      html_content = html_content,
      toc = toc
    )
  )
}



#' Process and Encode Data Field for Report
#' 
#' @noRd
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
#' @param enrichr_format A list with the formatted gene lists and background
#'                       gene list.
#' @param category_2_and_3_hits List of dataframes, where each df is the part
#' of the toptable that contains the significant features of the respective 
#' limma result category (2 or 3). Default is NA.
#'
#' @return A string containing the HTML link for downloading the processed
#'         field.
#'
#' @importFrom base64enc base64encode
#' @importFrom rstudioapi isAvailable getSourceEditorContext
#'
process_field <- function(
    field,
    data,
    meta,
    topTables,
    report_info,
    encode_df_to_base64,
    enrichr_format,
    category_2_and_3_hits = NA
    ) {
  
  if (field == "data_with_annotation") {
    base64_df <- sprintf(
      '<a href="%s" download="data.xlsx" class="embedded-file">
       <button>Download data_with_annotation.xlsx</button></a>',
      encode_df_to_base64(data)
    )
  } else if (field == "meta" &&
    !is.null(meta) &&
    is.data.frame(meta) &&
    !any(is.na(meta))) {
    base64_df <- sprintf(
      '<a href="%s" download="meta.xlsx" class="embedded-file">
       <button>Download meta.xlsx</button></a>',
      encode_df_to_base64(meta)
    )
  } else if (field == "limma_topTables_clustered_time_effect_hits" &&
    !any(is.na(topTables))) {
    base64_df <- sprintf(
      '<a href="%s" download="limma_topTables_clustered_time_effect_hits.xlsx"
      class="embedded-file">
       <button>Download limma_topTables_clustered_time_effect_hits.xlsx</button>
      </a>',
      encode_df_to_base64(topTables)
    )
  } else if (field == "limma_topTables" && !any(is.na(topTables))) {
    base64_df <- sprintf(
      '<a href="%s" download="limma_topTables.xlsx" class="embedded-file">
       <button>Download limma_topTables.xlsx</button></a>',
      encode_df_to_base64(topTables)
    )
  } else if (field == "limma_topTables_avrg_diff_conditions_hits") {
    base64_df <- sprintf(
      '<a href="%s" download="limma_topTables_avrg_diff_conditions_hits.xlsx"
      class="embedded-file">
       <button>Download limma_topTables_avrg_diff_conditions_hits.xlsx
      </button></a>',
      encode_df_to_base64(category_2_and_3_hits[[1]])   # category 2
    )
  } else if (field == "limma_topTables_interaction_condition_time_hits") {
    base64_df <- sprintf(
      '<a href="%s" 
      download="limma_topTables_interaction_condition_time_hits.xlsx"
      class="embedded-file">
       <button>Download limma_topTables_interaction_condition_time_hits.xlsx
      </button></a>',
      encode_df_to_base64(category_2_and_3_hits[[2]])   # category 3
    )
  } else if (field == "Enrichr_clustered_genes" &&
    !any(is.na(enrichr_format)) &&
    !is.null(enrichr_format$gene_lists)) {
    # Create ZIP file for Enrichr_clustered_genes
    zip_base64 <- create_enrichr_zip(enrichr_format)
    base64_df <- sprintf(
      '<a href="data:application/zip;base64,%s"
      download="Enrichr_clustered_genes.zip" class="embedded-file">
       <button>Download Enrichr_clustered_genes.zip</button></a>',
      zip_base64
    )
  } else if (field == "Enrichr_background" &&
    !any(is.na(enrichr_format)) &&
    !is.null(enrichr_format$background)) {
    base64_df <- sprintf(
      '<a href="data:text/plain;base64,%s" download="Enrichr_background.txt"
      class="embedded-file">
       <button>Download Enrichr_background.txt</button></a>',
      base64enc::base64encode(charToRaw(enrichr_format$background))
    )
  } else if (field == "session_info") {
    # Capture session info and encode to base64 on-the-fly
    session_details <- sessionInfo()
    session_info <- paste(capture.output(session_details), collapse = "\n")
    base64_session_info <- base64enc::base64encode(charToRaw(session_info))

    base64_df <- sprintf(
      '<a href="data:text/plain;base64,%s" download="session_info.txt"
      class="embedded-file">
       <button>Download R Session Info</button></a>',
      base64_session_info
    )
  } else if (field == "analysis_script") {
    # Capture the analysis script content if available (in RStudio)
    if (rstudioapi::isAvailable()) {
      script_content <- rstudioapi::getSourceEditorContext()$contents
      base64_script <- base64enc::base64encode(
        charToRaw(
          paste(
            script_content,
            collapse = "\n"
            )
          )
        )
      
      base64_df <- sprintf(
        '<a href="data:text/plain;base64,%s" download="analysis_script.txt"
        class="embedded-file">
         <button>Download Analysis Script</button></a>',
        base64_script
      )
    } else {
      base64_df <- "Analysis script is unavailable."
    }
    
  } else {
    base64_df <- ifelse(
      is.null(report_info[[field]]),
      "NA", report_info[[field]]
      )
  }
  return(base64_df)
}



# Level 3 internal functions ---------------------------------------------------


#' Extract and Combine Substrings from a String
#' 
#' @noRd
#'
#' @description
#' The `extract_and_combine()` function extracts specific substrings from an 
#' input string based on predefined patterns, combines them with a whitespace 
#' in between, and optionally truncates the result to 30 characters.
#'
#' @param input A character string containing the substrings to be extracted. 
#'   The string should follow a specific format where the desired substrings 
#'   are marked by the patterns `cluster:` and `database:`.
#'
#' @details
#' The function performs the following steps:
#' 1. Extracts the substring that appears after `cluster:` and before the next 
#'    comma (`,`).
#' 2. Extracts the substring that appears after `database:` and before the next 
#'    comma (`,`).
#' 3. Combines the two extracted substrings with a whitespace between them.
#' 4. If the combined string exceeds 30 characters, it is truncated to the 
#'    first 30 characters.
#'
#' This is useful for processing structured input strings and formatting the 
#' results for display or further analysis.
#'
#' @return
#' A character string combining the extracted substrings, truncated to 
#' 30 characters if necessary.
#'
extract_and_combine <- function(input) {
  # Extract substring after 'cluster:' and before the next ','
  cluster_match <- regmatches(
    input,
    regexpr(
      "(?<=cluster: )[^,]+",
      input,
      perl = TRUE
    )
  )

  # Extract substring after 'database:'
  database_match <- regmatches(
    input,
    regexpr(
      "(?<=database: )[^,]+",
      input,
      perl = TRUE
    )
  )

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
#' @noRd
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


#' Add Plot to HTML Content
#' 
#' @noRd
#'
#' @description
#' This function converts a plot to a base64 image and adds it to the
#' HTML content.
#'
#' @param html_content The current HTML content as a character string.
#' @param plot_element The plot element to be converted to base64.
#' @param plots_size An integer specifying the height of the plot.
#' @param section_index An integer specifying the section index.
#'
#' @return The updated HTML content as a character string.
#'
add_plot_to_html <- function(
    html_content,
    plot_element,
    plots_size,
    section_index) {
  img_tag <- plot2base64(
    plot_element,
    height = plots_size
  )
  paste(
    html_content,
    '<div id="section', section_index, '">',
    img_tag,
    "</div>",
    '<hr style="border: 0; border-top: 1px solid #ccc; margin: 20px 0;">',
    sep = "\n"
  )
}
