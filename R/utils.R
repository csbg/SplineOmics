# Level 1 internal functions ---------------------------------------------------


check_spline_params <- function(spline_params, mode) {
  
  check_spline_params_generally(spline_params)
  check_spline_params_mode_dependent(spline_params, mode)
}


validate_report_info <- function(report_info) {
  
  mandatory_fields <- c("omics_data_type", "data_description", 
                        "data_collection_date", "analyst_name", 
                        "project_name", "dataset_name", "limma_design")
  
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


#' Make dynamic HTML text descriptions, such as header- and design_text.
#' Future task
#' 
#' @importFrom here here
#' 
generate_report_html_new <- function(plots, 
                                 plot_len, 
                                 report_dir,
                                 filename, 
                                 timestamp) {
  header_text <- paste("limma hyperparams screen", timestamp, sep=" | ")
  
  design_text <- ""
  
  html_content <- paste(
    "<html><head><title>My Plots</title></head><body>",
    "<h1 style='color:red;'>", header_text, "</h1>",
    "<h2></h2>",
    "<p>", design_text, "</p>",
    "</body></html>",
    sep=""
  )
  
  file_name <- sprintf("%s_%s.html", filename, timestamp)
  output_file_path <- here::here(report_dir, file_name)
  
  build_plot_report_html(html_content, 
                         plots, 
                         plot_len, 
                         output_file_path)
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
                                   output_file_path) {

  header_section <- paste(header_section,
                          "<div style='margin-top: 30px; margin-bottom: 30px;
                          '><hr></div>",
                          sep="\n")

  for (index in seq_along(plots)) {
    plot <- plots[[index]]
    plot_size <- plots_sizes[[index]]
    img_tag <- plot2base64(plot, plot_size)
    header_section <- paste(header_section, img_tag, sep="\n")
  }

  html_content <- paste(header_section, "</body></html>", sep="\n")
  
  dir_path <- dirname(output_file_path)
  
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  writeLines(html_content, output_file_path)
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
