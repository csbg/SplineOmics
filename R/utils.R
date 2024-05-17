#' utils.R contains shared functions that are used by at least two package 
#' functions of the splinetime package.

# Level 1 internal functions ---------------------------------------------------


check_mode <- function(mode) {
  if (!((mode == "integrated") || (mode == "isolated"))) {
    stop("mode must be either integrated or isolated. This is dependent on the
         used limma design formula. For example, this formula: 
         ~ 1 + Phase*X + Reactor would require mode = integrated, whereas this
         formula: ~ 1 + X + Reactor would require mode = isolated.")
  } else {
    sprintf("Mode %s choosen.", mode)
  }
}


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
                                 report_type = "limma_hyperparams_screen",
                                 mode = NA,
                                 filename = "report",
                                 timestamp = format(Sys.time(), 
                                                    "%d_%m_%Y-%H_%M_%S"),
                                 report_dir = here::here()) {
  
  if (report_type == "limma_hyperparams_screen") {
    title <- "hyperparams screen"
  } else if (report_type == "cluster_hits") {                         
    title <- "clustered hits"
  } else {
    stop("report_type must be either limma_hyperparams_screen or cluster_hits")
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
  
  if (report_type == "limma_hyperparams_screen") {
    build_hyperparams_screen_report(header_section = header_section, 
                                    plots = plots, 
                                    plots_sizes = plots_sizes, 
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


check_spline_params_generally <- function(spline_params) {

  if ("spline_type" %in% names(spline_params)) {
    if (!all(spline_params$spline_type %in% c("b", "n"))) {
      stop("Elements of spline_type must be either 'b' for B-splines or 'n'. for
             natural cubic splines, in spline_params")
    }
  } else {
    stop("spline_type is missing in spline_params.")
  }
  
  # Check if degree exists and is an integer vector
  if ("degree" %in% names(spline_params)) {
    if (!all(spline_params$degree == as.integer(spline_params$degree))) {
      stop("degree must be an integer vector in spline_params.")
    }
  } else if (!all(spline_params$spline_type %in% c("n"))) {
    stop("degree is missing in spline_params.")
  }
  
  if (("dof" %in% names(spline_params)) && 
      ("knots" %in% names(spline_params))) {
    stop("Either dof or knots must be present, but not both,in spline_params.")
  } else if (!("dof" %in% names(spline_params)) && 
             !("knots" %in% names(spline_params))) {
    stop("At least one of dof or knots must be present, in spline_params.")
  }
  
  # Check if dof exists and is an integer vector
  if ("dof" %in% names(spline_params)) {
    if (!all(spline_params$dof == as.integer(spline_params$dof))) {
      stop("dof must be an integer vector in spline_params.")
    }
    for (i in seq_along(spline_params$spline_type)) {
      if (spline_params$spline_type[i] == "b" && spline_params$dof[i] < 3) {
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


create_progress_bar <- function(iterable) {
  library(progress)
  
  # Create and return the progress bar
  pb <- progress_bar$new(
    format = "  Processing [:bar] :percent :elapsed",
    total = length(iterable),
    width = 60
  )
  
  return(pb)
}
