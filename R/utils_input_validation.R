#' utils scripts contains shared functions that are used by at least two package 
#' functions of the SplineOmics package.

# Level 1 internal functions ---------------------------------------------------


#' Check Mode
#'
#' @description
#' Validates that the mode is either 'integrated' or 'isolated', which depends 
#' on the design formula used in limma.
#'
#' @param mode A character string specifying the mode.
#'
#' @return A message indicating the chosen mode if valid; otherwise, an error 
#' is thrown.
#'
#' @examples
#' check_mode("integrated")
#' check_mode("isolated")
#'
#' @seealso
#' \code{\link{limma}}
#' 
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


#' Check Spline Parameters
#'
#' @description
#' Validates the spline parameters both generally and depending on the 
#' specified mode.
#'
#' @param spline_params A list of spline parameters.
#' @param mode A character string specifying the mode
#'             ('integrated' or 'isolated').
#'
#' @return No return value, called for side effects.
#'
#' @examples
#' spline_params <- list(spline_type = c("n"), dof = list(3))
#' mode <- "integrated"
#' check_spline_params(spline_params, mode)
#'
#' @seealso
#' \code{\link{check_spline_params_generally}}, 
#' \code{\link{check_spline_params_mode_dependent}}
#' 
check_spline_params <- function(spline_params, 
                                mode) {
  
  check_spline_params_generally(spline_params)
  check_spline_params_mode_dependent(spline_params, mode)
}


check_condition <- function(condition, 
                            meta) {
  
  # Check if condition is a single character
  if (!is.character(condition) || length(condition) != 1) {
    stop("'condition' must be a single character")
  }
  
  # Check if condition is a column in the meta dataframe
  if (!condition %in% colnames(meta)) {
    stop(sprintf("The specified condition '%s' must be a column name in the 
                 meta dataframe", condition))
  }
  
  return(TRUE)
}


check_adj_pthresholds <- function(adj_pthresholds) {
  
  if (!is.numeric(adj_pthresholds)) {
    stop("'adj_pthresholds' must be a numeric vector.")
  }
  
  # Check for elements <= 0
  if (any(adj_pthresholds <= 0)) {
    offending_elements <- which(adj_pthresholds <= 0)
    stop(paste0("'adj_pthresholds' must have all elements > 0. ",
                "Offending elements at indices: ",
                paste(offending_elements, collapse = ", "), 
                ". Values: ",
                paste(adj_pthresholds[offending_elements], collapse = ", "), 
                "."))
  }
  
  # Check for elements >= 1
  if (any(adj_pthresholds >= 1)) {
    offending_elements <- which(adj_pthresholds >= 1)
    stop(paste0("'adj_pthresholds' must have all elements < 1. ",
                "Offending elements at indices: ",
                paste(offending_elements, collapse = ", "), 
                ". Values: ",
                paste(adj_pthresholds[offending_elements], collapse = ", "), 
                "."))
  }
  
  return(TRUE)
}


check_time_unit <- function(time_unit) {
  
  if (!time_unit %in% c("s", "m", "h", "d")) {
    stop(paste0("time_unit must be one of s for seconds, m for minutes, h ", 
                "for hours, or d for days. Default is minutes."))
  }
}


check_and_create_report_dir <- function(report_dir) {
  
  # Attempt to create the directory if it does not exist
  if (!file.exists(report_dir)) {
    tryCatch({
      dir.create(report_dir, recursive = TRUE)
    }, warning = function(w) {
      stop(sprintf("Warning occurred while creating the directory: %s", w$message))
    }, error = function(e) {
      stop(sprintf("Error occurred while creating the directory: %s", e$message))
    })
  }
  
  # Verify that the directory exists and is a directory
  if (!file.exists(report_dir) || !file.info(report_dir)$isdir) {
    stop(sprintf("The specified path is not a valid directory: %s", report_dir))
  }
  
  return(TRUE)
}


check_padjust_method <- function(padjust_method) {
  
  supported_methods <- c("holm", "hochberg", "hommel", "bonferroni", "BH", 
                         "BY", "fdr", "none")
  if (!(is.character(padjust_method) && 
        padjust_method %in% supported_methods)) 
  {
    stop(sprintf("padjust_method must be a character and one of the 
                   supported methods (%s).",
                 paste(supported_methods, collapse = ", ")))
  }
}


#' Check Report Information
#'
#' @description
#' Validates the report information to ensure it contains all mandatory fields 
#' and adheres to the required formats.
#'
#' @param report_info A named list containing report information.
#'
#' @return TRUE if the report information is valid; otherwise, an error is 
#' thrown.
#'
#' @examples
#' report_info <- list(
#'   omics_data_type = "genomics",
#'   data_description = "Sample description",
#'   data_collection_date = "2023-01-01",
#'   analyst_name = "John Doe",
#'   project_name = "Project XYZ"
#' )
#' check_report_info(report_info)
#' 
check_report_info <- function(report_info) {
  
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


# Level 2 internal functions ---------------------------------------------------


#' Check Spline Parameters Generally
#'
#' @description
#' Validates the general structure and contents of spline parameters.
#'
#' @param spline_params A list of spline parameters.
#'
#' @return No return value, called for side effects.
#'
#' @examples
#' spline_params <- list(
#'   spline_type = c("b", "n"),
#'   degree = c(3, NA),
#'   dof = c(4, NA),
#'   knots = list(NA, c(1, 2, 3)),
#'   bknots = list(NA, c(0.5, 1.5))
#' )
#' check_spline_params_generally(spline_params)
#'
#' @seealso
#' \code{\link{check_spline_params_mode_dependent}}
#' 
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


#' Check Spline Parameters Mode Dependent
#'
#' @description
#' Validates the spline parameters depending on the specified mode.
#'
#' @param spline_params A list of spline parameters.
#' @param mode A character string specifying the mode 
#'            ('integrated' or 'isolated').
#'
#' @return No return value, called for side effects.
#'
#' @examples
#' spline_params <- list(
#'   spline_type = "b",
#'   degree = 3,
#'   dof = 4,
#'   knots = list(c(1, 2, 3)),
#'   bknots = list(c(0.5, 1.5))
#' )
#' check_spline_params_mode_dependent(spline_params, "integrated")
#'
#' @seealso
#' \code{\link{check_spline_params_generally}}
#' 
check_spline_params_mode_dependent <- function(spline_params, 
                                               mode) {
  
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
