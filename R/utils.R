# Level 1 internal functions ---------------------------------------------------


check_spline_params <- function(spline_params, mode) {
  
  check_spline_params_generally(spline_params)
  check_spline_params_mode_dependent(spline_params, mode)
}


validate_report_info <- function(report_info) {
  
  mandatory_fields <- c("omics_data_type", "data_description", 
                        "data_collection_date", "analyst_name", 
                        "project_name")
  
  # Check if report_info is a named list
  if (!is.list(report_info) || is.null(names(report_info))) {
    stop("report_info must be a named list.")
  }
  
  # Check if all values in report_info are strings
  non_string_fields <- sapply(report_info, function(x) !is.character(x))
  if (any(non_string_fields)) {
    invalid_fields <- names(report_info)[non_string_fields]
    stop(paste("The following fields must be strings:", paste(invalid_fields, 
                                                              collapse = ", ")))
  }
  
  # Check if all mandatory fields are present
  missing_fields <- setdiff(mandatory_fields, names(report_info))
  if (length(missing_fields) > 0) {
    stop(paste("Missing mandatory fields:", paste(missing_fields, 
                                                  collapse = ", ")))
  }
  
  if (!grepl("^[a-zA-Z_]+$", report_info[["omics_data_type"]])) {
    stop("The 'omics_data_type' field must contain only alphabetic letters 
         and underscores.")
  }
  
  long_fields <- sapply(report_info, function(x) any(nchar(x) > 110))
  if (any(long_fields)) {
    too_long_fields <- names(report_info)[long_fields]
    stop(paste("The following fields have strings exceeding 110 characters:", 
               paste(too_long_fields, collapse = ", ")))
  }

  return(TRUE)
}


#' Generate HTML Report with Plot List
#'
#' This function generates an HTML report for a given list of plots. It 
#' includes
#' headers and design details specific to the data type being analyzed, and 
#' saves
#' the report in the specified directory. The function is internal and not 
#' exported.
#'
#' @param plots A list of ggplot objects to be included in the report.
#' @param plots_sizes The number of rows to arrange the plots in the report.
#' @param report_dir The directory where the report HTML file will be saved.
#' @importFrom here here
#' @importFrom tools file_path_sans_ext
#' @importFrom grDevices dev.off
#' @importFrom htmltools save_html
#' 
generate_report_html <- function(plots, 
                                 plots_sizes, 
                                 level_headers_info,
                                 spline_params,
                                 report_info,
                                 mode = NA,
                                 filename = "report",
                                 timestamp = format(Sys.time(), 
                                                    "%d_%m_%Y-%H_%M_%S"),
                                 report_dir = here::here()) {
  
  header_text <- paste("clustered hits", 
                       report_info$omics_data_type, 
                       timestamp, sep=" | ")
  
  header_section <- paste(
    "<html><head><title>clustered_hits</title>",
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
  
  build_plot_report_html(header_section = header_section, 
                         plots = plots, 
                         plots_sizes = plots_sizes,
                         level_headers_info = level_headers_info,
                         spline_params = spline_params,
                         mode = mode,
                         output_file_path = output_file_path)
}


# Level 2 internal functions ---------------------------------------------------


check_spline_params_generally <- function(spline_params) {

  if ("spline_type" %in% names(spline_params)) {
    if (!all(spline_params$spline_type %in% c("b", "n"))) {
      stop("Elements of spline_type must be either 'b' for B-splines or 'n'. for
             natural cubic splines, in spline_params")
    }
  } else {
    stop("spline_type is missing in spline_params.")
  }
  
  # Check if degrees exists and is an integer vector
  if ("degrees" %in% names(spline_params)) {
    if (!all(spline_params$degrees == as.integer(spline_params$degrees))) {
      stop("degrees must be an integer vector in spline_params.")
    }
  } else if (!all(spline_params$spline_type %in% c("n"))) {
    stop("degrees is missing in spline_params.")
  }
  
  # Check if DoFs exists and is an integer vector
  if ("DoFs" %in% names(spline_params)) {
    if (!all(spline_params$DoFs == as.integer(spline_params$DoFs))) {
      stop("DoFs must be an integer vector in spline_params.")
    }
    
    
    for (i in seq_along(spline_params$spline_type)) {
      if (spline_params$spline_type[i] == "b" && spline_params$DoFs[i] < 3) {
        stop(paste0("B-splines require DoF > 2, for spline_params spline_type ", 
                    "index ", i))
      }
    }
  }
  
  # Check if knots exists and is a list of numeric vectors
  if ("knots" %in% names(spline_params)) {
    if (!is.list(spline_params$knots) || 
        any(sapply(spline_params$knots, function(x) !is.numeric(x)))) {
      stop("knots must be a list of numeric vectors in spline_params.")
    }
  }
  
  if (("DoFs" %in% names(spline_params)) && 
      ("knots" %in% names(spline_params))) {
    stop("Either DoFs or knots must be present, but not both,in spline_params.")
  } else if (!("DoFs" %in% names(spline_params)) && 
             !("knots" %in% names(spline_params))) {
    stop("At least one of DoFs or knots must be present, in spline_params.")
  }
  
  # Check if bknots exists and is a list of numeric vectors
  if ("bknots" %in% names(spline_params)) {
    if (!is.list(spline_params$bknots) || 
        any(sapply(spline_params$bknots, function(x) !is.numeric(x)))) {
      stop("bknots must be a list of numeric vectors, in spline_params.")
    }
  }
}


check_spline_params_mode_dependent <- function(spline_params, mode) {
  
  if (mode == "integrated") {
    # Check that each vector in the main parameters has exactly one element
    if (any(sapply(spline_params, function(x) length(x) != 1))) {
      stop("All parameters in spline_params must have exactly one element when
            mode is 'integrated'. Different spline parameters for the different 
            levels is not supported for this mode")
    }
    
    # Additional check for 'knots' and 'bknots' if they exist
    if ("knots" %in% names(spline_params)) {
      if (length(spline_params$knots) != 1) {
        stop("All elements in 'knots' in spline_params must have length 1 when 
              mode is 'integrated'. Different spline parameters for the 
              different levels is not supported for this mode")
      }
    }
    if ("bknots" %in% names(spline_params)) {
      if (length(spline_params$bknots) != 1) {
        stop("All elements in 'bknots' in spline_params must have length 1 
              when mode is 'integrated'. Different spline parameters for the 
              different levels is not supported for this mode")
      }
    }
  } else if (mode == "isolated") {
    num_levels <- length(unique(meta[[condition]]))
    if (any(sapply(spline_params, length) != num_levels)) {
      stop("Each vector or list in spline_params must have as many elements as 
           there are unique elements in the ",
           condition, " column of meta when mode is 'isolated'.")
    }
    if ("knots" %in% names(spline_params)) {
      if (length(spline_params$knots) != num_levels) {
        stop("'knots' in spline_params must have the same number of elements as 
             there are unique elements in the ",
             condition, " column of meta when mode is 'isolated'.")
      }
    }
    if ("bknots" %in% names(spline_params)) {
      if (length(spline_params$bknots) != num_levels) {
        stop("'bknots' in spline_params must have the same number of elements as
             there are unique elements in the ",
             condition, " column of meta when mode is 'isolated'.")
      }
    }
  }
}


build_plot_report_html <- function(header_section, 
                                   plots, 
                                   plots_sizes, 
                                   level_headers_info = level_headers_info,
                                   spline_params = spline_params,
                                   mode,
                                   output_file_path) {  
  
  content_with_plots <- paste(header_section, "<!--TOC-->", sep="\n")
  
  # toc <- "<div id='toc'><h2 style='font-size: 40px;'>Table of Contents</h2><ul>"
  toc <- "<div id='toc' style='text-align: center; display: block; margin: auto; 
          width: 80%;'> <h2 style='font-size: 40px;'>Table of Contents</h2><ul>"
  
  
  section_header_style <- "font-size: 70px; color: #001F3F; text-align: center;"
  toc_style <- "font-size: 30px;"
  
  current_header_index <- 1
  
  # Generate the sections and plots
  for (index in seq_along(plots)) {
    if (current_header_index <= length(level_headers_info)) {
      header_info <- level_headers_info[[current_header_index]]
      header_placement <- header_info[[2]]
      
      if (index - 1 == header_placement) {
        section_header <- sprintf("<h2 style='%s' id='section%d'>%s</h2>", 
                                  section_header_style, index, header_info[[1]])
        content_with_plots <- paste(content_with_plots, section_header, 
                                    sep="\n")
        
        if (mode == "integrated") {
          j <- 1
        } else {      # mode == "isolated" or mode == NA
          j <- index
        }
        
        if (!is.null(spline_params$spline_type) && 
            length(spline_params$spline_type) >= j) {
          spline_params$spline_type[j] <- spline_params$spline_type[j]
        } else {
          spline_params$spline_type[j] <- NA
        }
        
        if (!is.null(spline_params$degrees) && 
            length(spline_params$degrees) >= j) {
          spline_params$degrees[j] <- spline_params$degrees[j]
        } else {
          spline_params$degrees[j] <- NA
        }
        
        if (!is.null(spline_params$DoFs) && 
            length(spline_params$DoFs) >= j) {
          spline_params$DoFs[j] <- spline_params$DoFs[j]
        } else {
          spline_params$DoFs[j] <- NA
        }
        
        if (!is.null(spline_params$knots) && 
            length(spline_params$knots) >= j) {
          spline_params$knots[j] <- spline_params$knots[j]
        } else {
          spline_params$knots[j] <- NA
        }
        
        if (!is.null(spline_params$bknots) && 
            length(spline_params$bknots) >= j) {
          spline_params$bknots[j] <- spline_params$bknots[j]
        } else {
          spline_params$bknots[j] <- NA
        }
        
        
        if (spline_params$spline_type[j] == "b") {
          spline_params_info <- 
            sprintf("<p style='text-align: center; font-size: 30px;'>
                    <span style='color: blue;'>Spline-type:</span> B-spline<br>
                    <span style='color: blue;'>Degree:</span> %s<br>
                    <span style='color: blue;'>DoF:</span> %s<br>
                    <span style='color: blue;'>Knots:</span> %s<br>
                    <span style='color: blue;'>Boundary-knots:</span> %s</p>", 
                    spline_params$degrees[j], spline_params$DoFs[j], 
                    spline_params$knots[j], spline_params$bknots[j])
        } else {    # == "n"
          spline_params_info <- 
            sprintf("<p style='text-align: center; font-size: 30px;'>
                    <span style='color: blue;'>Spline-type:</span> Natural cubic
                    spline<br>
                    <span style='color: blue;'>DoF:</span> %s<br>
                    <span style='color: blue;'>Knots:</span> %s<br>
                    <span style='color: blue;'>Boundary-knots:</span> %s</p>", 
                    spline_params$DoFs[j], spline_params$knots[j], 
                    spline_params$bknots[j])
          
        }
        
        content_with_plots <- paste(content_with_plots, spline_params_info, 
                                    sep="\n")
        
        toc_entry <- sprintf("<li style='%s'><a href='#section%d'>%s</a></li>", 
                             toc_style, index, header_info[[1]])
        toc <- paste(toc, toc_entry, sep="\n")
        
        current_header_index <- current_header_index + 1
      }
    }
    
    # Process each plot
    plot <- plots[[index]]
    plot_size <- plots_sizes[[index]]
    img_tag <- plot2base64(plot, plot_size)
    content_with_plots <- paste(content_with_plots, img_tag, sep="\n")
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




# Level 2 internal functions ---------------------------------------------------


#' @importFrom ggplot2 ggsave
#' @importFrom base64enc dataURI
#' @importFrom ragg agg_png
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
