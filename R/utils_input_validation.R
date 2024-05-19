#' utils scripts contains shared functions that are used by at least two package 
#' functions of the SplineOmics package.

# Level 1 internal functions ---------------------------------------------------


#' Check Data and Meta
#'
#' @description
#' This function checks the validity of the data and meta objects, ensuring that 
#' data is a matrix with numeric values and that meta is a dataframe containing 
#' the specified condition column. Additionally, it verifies that the number of 
#' columns in the data matrix matches the number of rows in the meta dataframe.
#'
#' @param data A matrix containing numeric values. 
#' @param meta A dataframe containing the metadata, including the 'Time' column 
#' and the specified condition column.
#' @param condition A single character string specifying the column name in 
#' the meta dataframe to be checked.
#' @param meta_batch_column An optional parameter specifying the column name 
#' in the meta dataframe used to remove the batch effect. Default is NA.
#' @param data_meta_index An optional parameter specifying the index of the 
#' data/meta pair for error messages. Default is NA.
#'
#' @return Returns TRUE if all checks pass. Stops execution and returns an 
#' error message if any check fails.
#'
#' @examples
#' # Example of how to use the function
#' data <- matrix(runif(9), nrow = 3, ncol = 3)
#' meta <- data.frame(Time = 1:3, Condition = c("A", "B", "C"))
#' check_data_and_meta(data, meta, "Condition")
#'
#' @export
#' 
check_data_and_meta <- function(data,
                                meta,
                                condition,
                                meta_batch_column = NA,
                                data_meta_index = NA) {
  
  check_data(data,
             data_meta_index)
  
  check_meta(meta,
             condition,
             meta_batch_column,
             data_meta_index)
  
  if (!(nrow(meta) == ncol(data))) {
    if (!is.na(index)) {
      stop(paste0("For index ", index, "data column number must be equal to ",
                  "meta row number"),
           call. = FALSE)
    } else {
      stop(paste0("data column number must be equal to meta row number"),
           call. = FALSE)
    }
  }
}


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
         formula: ~ 1 + X + Reactor would require mode = isolated.",
         call. = FALSE)
  } else {
    sprintf("Mode %s choosen.", mode)
  }
}


#' Check Design Formula
#'
#' @description
#' Validates the design formula ensuring it is a valid character string, 
#' contains allowed characters, includes the intercept term 'X', and references 
#' columns present in the metadata.
#'
#' @param formula A character string representing the design formula.
#' @param meta A data frame containing metadata.
#' @param meta_index An optional index for the data/meta pair.
#'
#' @return TRUE if the design formula is valid, otherwise an error is thrown.
#'
#' @examples
#' meta <- data.frame(Time = seq(1, 10), condition = rep(c("A", "B"), each = 5))
#' check_design_formula("~ Time + condition * X", meta)
#'
#' @seealso
#' \code{\link{stats::model.matrix}}
#' 
check_design_formula <- function(formula, 
                                 meta,
                                 meta_index = NA) {
  
  # Check if the formula is a valid character string
  if (!is.character(formula) || length(formula) != 1) {
    stop("The design formula must be a valid character string.",
         call. = FALSE)
  }
  
  # Ensure the formula contains allowed characters only
  allowed_chars <- "^[~ 1A-Za-z0-9_+*:()-]*$"
  if (!grepl(allowed_chars, formula)) {
    stop("The design formula contains invalid characters.",
         call. = FALSE)
  }
  
  # Ensure the formula contains the intercept term 'X'
  if (!grepl("\\bX\\b", formula)) {
    stop("The design formula must include the term 'X'.",
         call. = FALSE)
  }
  
  # Extract terms from the formula (removing interactions and functions)
  formula_terms <- unlist(strsplit(gsub("[~+*:()]", " ", formula), " "))
  formula_terms <- formula_terms[formula_terms != ""]
  
  # Remove '1' and 'X' from terms since they are not columns
  formula_terms <- setdiff(formula_terms, c("1", "X"))
  
  # Check if the terms are present in the dataframe
  missing_columns <- setdiff(formula_terms, names(meta))
  if (length(missing_columns) > 0) {
    if (!is.na(meta_index)) {
      stop(sprintf("%s (data/meta pair index: %s): %s",
                   "The following design columns are missing in meta",
                   meta_index, 
                   paste(missing_columns, collapse = ", ")),
           call. = FALSE)
      
    } else {
      stop(paste("The following design columns are missing in meta:", 
                 paste(missing_columns, collapse = ", ")),
           call. = FALSE)
    }
  }
  
  return(TRUE)
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


#' Check Adjusted p-Thresholds
#'
#' @description
#' This function checks the validity of the adjusted p-thresholds vector, 
#' ensuring that 
#' all elements are numeric, greater than 0, and less than 1. If any of these 
#' conditions 
#' are not met, the function stops execution and returns an error message 
#' indicating the 
#' offending elements.
#'
#' @param adj_pthresholds A numeric vector of adjusted p-thresholds.
#'
#' @return Returns TRUE if all checks pass. Stops execution and returns an 
#' error message if any check fails.
#'
#' @examples
#' # Example of how to use the function
#' valid_thresholds <- c(0.01, 0.05, 0.1)
#' check_adj_pthresholds(valid_thresholds)
#'
#' invalid_thresholds <- c(0.01, -0.05, 1.1)
#' try(check_adj_pthresholds(invalid_thresholds))  # Should raise an error
#'
#' @export
#' 
check_adj_pthresholds <- function(adj_pthresholds) {
  
  if (!is.numeric(adj_pthresholds)) {
    stop("'adj_pthresholds' must be a numeric vector.",
         call. = FALSE)
  }
  
  # Check for elements <= 0
  if (any(adj_pthresholds <= 0)) {
    offending_elements <- which(adj_pthresholds <= 0)
    stop(paste0("'adj_pthresholds' must have all elements > 0. ",
                "Offending elements at indices: ",
                paste(offending_elements, collapse = ", "), 
                ". Values: ",
                paste(adj_pthresholds[offending_elements], collapse = ", "), 
                "."),
         call. = FALSE)
  }
  
  # Check for elements >= 1
  if (any(adj_pthresholds >= 1)) {
    offending_elements <- which(adj_pthresholds >= 1)
    stop(paste0("'adj_pthresholds' must have all elements < 1. ",
                "Offending elements at indices: ",
                paste(offending_elements, collapse = ", "), 
                ". Values: ",
                paste(adj_pthresholds[offending_elements], collapse = ", "), 
                "."),
         call. = FALSE)
  }
  
  return(TRUE)
}


#' Check Time Unit
#'
#' @description
#' This function checks if the provided time unit is valid. The valid time 
#' units are:
#' 's' for seconds, 'm' for minutes, 'h' for hours, and 'd' for days. If the 
#' time unit 
#' is not one of these, the function stops execution and returns an error 
#' message.
#'
#' @param time_unit A character string specifying the time unit. Valid options 
#' are 
#' 's' for seconds, 'm' for minutes, 'h' for hours, and 'd' for days.
#'
#' @return Returns TRUE if the time unit is valid. Stops execution and returns 
#' an error message if the time unit is invalid.
#'
#' @examples
#' # Example of how to use the function
#' check_time_unit("m")  # Should pass without error
#' 
#' try(check_time_unit("x"))  # Should raise an error
#'
#' @export
#' 
check_time_unit <- function(time_unit) {
  
  if (!time_unit %in% c("s", "m", "h", "d")) {
    stop(paste0("time_unit must be one of s for seconds, m for minutes, h ", 
                "for hours, or d for days. Default is minutes."),
         call. = FALSE)
  }
}


#' Check and Create Report Directory
#'
#' @description
#' This function checks if the specified report directory exists and is a 
#' valid directory. 
#' If the directory does not exist, it attempts to create it. If there are any 
#' warnings or 
#' errors during directory creation, the function stops execution and returns 
#' an error message.
#'
#' @param report_dir A character string specifying the path to the report 
#' directory.
#'
#' @return Returns TRUE if the directory exists or is successfully created. 
#' Stops execution 
#' and returns an error message if the directory cannot be created or is not 
#' valid.
#'
#' @examples
#' # Example of how to use the function
#' report_dir <- "path/to/report_dir"
#' check_and_create_report_dir(report_dir)
#'
#' @export
#' 
check_and_create_report_dir <- function(report_dir) {
  
  # Attempt to create the directory if it does not exist
  if (!file.exists(report_dir)) {
    tryCatch({
      dir.create(report_dir, recursive = TRUE)
    }, warning = function(w) {
      stop(sprintf("Warning occurred while creating the directory: %s", 
                   w$message),
           call. = FALSE)
    }, error = function(e) {
      stop(sprintf("Error occurred while creating the directory: %s", 
                   e$message),
           call. = FALSE)
    })
  }
  
  # Verify that the directory exists and is a directory
  if (!file.exists(report_dir) || !file.info(report_dir)$isdir) {
    stop(sprintf("The specified path is not a valid directory: %s", 
                 report_dir),
         call. = FALSE)
  }
  
  return(TRUE)
}


#' Check p-Adjustment Method
#'
#' @description
#' This function checks if the provided p-adjustment method is valid. The valid 
#' methods are:
#' "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", and "none". 
#' If the method
#' is not one of these, the function stops execution and returns an error
#'  message.
#'
#' @param padjust_method A character string specifying the p-adjustment method. 
#' Valid options 
#' are "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", and 
#' "none".
#'
#' @return Returns TRUE if the p-adjustment method is valid. Stops execution and 
#' returns an error message if the method is invalid.
#'
#' @examples
#' # Example of how to use the function
#' check_padjust_method("BH")  # Should pass without error
#' 
#' try(check_padjust_method("invalid_method"))  # Should raise an error
#'
#' @export
#' 
check_padjust_method <- function(padjust_method) {
  
  supported_methods <- c("holm", "hochberg", "hommel", "bonferroni", "BH", 
                         "BY", "fdr", "none")
  if (!(is.character(padjust_method) && 
        padjust_method %in% supported_methods)) 
  {
    stop(sprintf("padjust_method must be a character and one of the 
                   supported methods (%s).",
                 paste(supported_methods, collapse = ", ")),
         call. = FALSE)
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
    stop("report_info must be a named list.",
         call. = FALSE)
  }
  
  # Check if all values in report_info are strings
  non_string_fields <- sapply(report_info, function(x) !is.character(x))
  if (any(non_string_fields)) {
    invalid_fields <- names(report_info)[non_string_fields]
    stop(paste("The following fields must be strings:", paste(invalid_fields, 
                                                              collapse = ", ")),
         call. = FALSE)
  }
  
  # Check if all mandatory fields are present
  missing_fields <- setdiff(mandatory_fields, names(report_info))
  if (length(missing_fields) > 0) {
    stop(paste("Missing mandatory fields:", paste(missing_fields, 
                                                  collapse = ", ")),
         call. = FALSE)
  }
  
  if (!grepl("^[a-zA-Z_]+$", report_info[["omics_data_type"]])) {
    stop(paste("The 'omics_data_type' field must contain only alphabetic",  
               "letters and underscores."),
         call. = FALSE)
  }
  
  long_fields <- sapply(report_info, function(x) any(nchar(x) > 110))
  if (any(long_fields)) {
    too_long_fields <- names(report_info)[long_fields]
    stop(paste("The following fields have strings exceeding 110 characters:", 
               paste(too_long_fields, collapse = ", ")),
         call. = FALSE)
  }
  
  return(TRUE)
}


# Level 2 internal functions ---------------------------------------------------


#' Check Data Matrix
#'
#' @description
#' This function checks the validity of the data matrix, ensuring that it is a 
#' matrix, contains only numeric values, 
#' has no missing values, and all elements are non-negative. Additionally, it
#'  verifies that no rows or columns are 
#' entirely zeros.
#'
#' @param data A matrix containing numeric values.
#' @param data_meta_index An optional parameter specifying the index of the data 
#' for error messages. Default is NA.
#'
#' @return Returns TRUE if all checks pass. Stops execution and returns an error
#'  message if any check fails.
#'
#' @examples
#' # Example of how to use the function
#' data <- matrix(runif(9), nrow = 3, ncol = 3)
#' check_data(data)
#' 
#' invalid_data <- matrix(c(1, 2, NA, 4, 5, 6, 7, 8, 9), nrow = 3, ncol = 3)
#' try(check_data(invalid_data))  # Should raise an error
#'
#' @export
#' 
check_data <- function(data,
                       data_meta_index = NA) {
  
  if (!is.matrix(data)) {
    stop(create_error_message("data must be a matrix.", 
                         data_meta_index),
         call. = FALSE)
  }
  
  if (!all(is.numeric(data))) {
    stop(create_error_message("data must contain only numeric values.", 
                              data_meta_index),
         call. = FALSE)
  }
  
  if (any(is.na(data))) {
    stop(create_error_message("data must not contain missing values.", 
                              data_meta_index),
         call. = FALSE)
  }
  
  if (any(data < 0)) {
    stop(create_error_message(paste("All elements of data must be non-negative.
                               The elements should represent concentrations, 
                               abudances, or intensities (which are inherently 
                               non-negative", 
                              data_meta_index)),
         call. = FALSE)
  }
  
  if (any(rowSums(data) == 0)) {
    stop(create_error_message("data must not contain rows with all zeros.", 
                              data_meta_index),
         call. = FALSE)
  }
  
  if (any(colSums(data) == 0)) {
    stop(create_error_message("data must not contain columns with all zeros.", 
                              data_meta_index),
         call. = FALSE)
  }
  
  return(TRUE)
}


#' Check Metadata
#'
#' @description
#' This function checks the validity of the metadata dataframe, ensuring it 
#' contains the 'Time' column, 
#' does not contain missing values, and that the specified condition column is 
#' valid and of the appropriate type.
#' Additionally, it checks for an optional batch effect column and prints 
#' messages regarding its use.
#'
#' @param meta A dataframe containing the metadata, including the 'Time' column.
#' @param condition A single character string specifying the column name in the 
#' meta dataframe to be checked.
#' @param meta_batch_column An optional parameter specifying the column name in
#'  the meta dataframe used to remove the batch effect. Default is NA.
#' @param data_meta_index An optional parameter specifying the index of the 
#' data/meta pair for error messages. Default is NA.
#'
#' @return Returns TRUE if all checks pass. Stops execution and returns an 
#' error message if any check fails.
#'
#' @examples
#' # Example of how to use the function
#' meta <- data.frame(Time = 1:3, Condition = c("A", "B", "C"))
#' check_meta(meta, "Condition")
#' 
#' invalid_meta <- data.frame(Time = c(1, NA, 3), Condition = c("A", "B", "C"))
#' try(check_meta(invalid_meta, "Condition"))  # Should raise an error
#'
#' @export
#' 
check_meta <- function(meta, 
                       condition,
                       meta_batch_column = NA,
                       data_meta_index = NA) {
  
  if (!is.data.frame(meta) || !"Time" %in% names(meta)) {
    stop(create_error_message("meta must be a dataframe with the column Time",
                              data_meta_index), 
         call. = FALSE)
  }
  
  if (any(is.na(meta))) {
    stop(create_error_message("meta must not contain missing values.",
                              data_meta_index), 
         call. = FALSE)
  }
  
  # Check if condition is a single character
  if (!is.character(condition) || length(condition) != 1) {
    stop("'condition' must be a single character",
         call. = FALSE)
  }
  
  # Check if condition is a column in the meta dataframe
  if (!condition %in% colnames(meta)) {
    stop(create_error_message(sprintf("The condition '%s' is not a %s", 
                                      condition, 
                                      paste("column in meta")),
                              data_meta_index), 
         call. = FALSE)
  }
  
  # Check if the factor column is of appropriate type
  if (!is.factor(meta[[condition]]) && 
      !is.character(meta[[condition]])) {
    
    stop(create_error_message(sprintf("The factor column '%s' must be of type 
                                      factor or character.", condition),
                              data_meta_index),
         call. = FALSE)
  }
  
  if (!is.na(meta_batch_column) && !(meta_batch_column %in% names(meta))) {
    stop(create_error_message(sprintf("Batch effect column '%s' %s", 
                                 meta_batch_column, "not found in meta"),
                              data_meta_index),
         call. = FALSE)
  } else if (!is.na(meta_batch_column)) {
    if (!is.na(data_meta_index)) {
      message(sprintf("Index: %s. %s", 
                      data_meta_index, 
                      paste("Column", meta_batch_column, "of meta will be used", 
                            "to remove the batch effect for the plotting")))
      
    } else {
      message(sprintf("Column '%s' of meta will be used to %s", 
                      meta_batch_column, 
                      paste("remove the batch effect for the plotting")))
    }
    
    

  } else {
    if (!is.na(data_meta_index)) {
      message(sprintf("Index: %s. Batch effect will not be removed for 
                      plotting!", 
                      data_meta_index))
    } else {
      message("Batch effect will not be removed for plotting!")
    }
  }
  
  return(TRUE)
}


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
      stop(paste("Elements of spline_type must be either 'b' for B-splines or",
                 "'n'. for natural cubic splines, in spline_params"),
           call. = FALSE)
    }
  } else {
    stop("spline_type is missing in spline_params.",
         call. = FALSE)
  }
  
  # Check if degree exists and is an integer vector
  if ("degree" %in% names(spline_params)) {
    if (!all(spline_params$degree == as.integer(spline_params$degree))) {
      stop("degree must be an integer vector in spline_params.",
           call. = FALSE)
    }
  } else if (!all(spline_params$spline_type %in% c("n"))) {
    stop("degree is missing in spline_params.",
         call. = FALSE)
  }
  
  if (("dof" %in% names(spline_params)) && 
      ("knots" %in% names(spline_params))) {
    stop("Either dof or knots must be present, but not both,in spline_params.",
         call. = FALSE)
  } else if (!("dof" %in% names(spline_params)) && 
             !("knots" %in% names(spline_params))) {
    stop("At least one of dof or knots must be present, in spline_params.",
         call. = FALSE)
  }
  
  # Check if dof exists and is an integer vector
  if ("dof" %in% names(spline_params)) {
    if (!all(spline_params$dof == as.integer(spline_params$dof))) {
      stop("dof must be an integer vector in spline_params.",
           call. = FALSE)
    }
    for (i in seq_along(spline_params$spline_type)) {
      if (spline_params$spline_type[i] == "b" && spline_params$dof[i] < 3) {
        stop(paste("B-splines require DoF > 2, for spline_params spline_type", 
                    "index", i),
             call. = FALSE)
      }
    }
  }
  
  # Check if knots exists and is a list of numeric vectors
  if ("knots" %in% names(spline_params)) {
    if (!is.list(spline_params$knots) || 
        any(sapply(spline_params$knots, function(x) !is.numeric(x)))) {
      stop("knots must be a list of numeric vectors in spline_params.",
           call. = FALSE)
    }
  }
  
  # Check if bknots exists and is a list of numeric vectors
  if ("bknots" %in% names(spline_params)) {
    if (!is.list(spline_params$bknots) || 
        any(sapply(spline_params$bknots, function(x) !is.numeric(x)))) {
      stop("bknots must be a list of numeric vectors, in spline_params.",
           call. = FALSE)
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
      stop(paste("All parameters in spline_params must have exactly one",
                 "element when mode is 'integrated'. Different spline", 
                 "parameters for the different levels is not supported for",
                 "this mode"),
           call. = FALSE)
    }
    
    # Additional check for 'knots' and 'bknots' if they exist
    if ("knots" %in% names(spline_params)) {
      if (length(spline_params$knots) != 1) {
        stop(paste("All elements in 'knots' in spline_params must have length",
                   "1 when mode is 'integrated'. Different spline parameters",
                   "for the different levels is not supported for this mode"),
             call. = FALSE)
      }
    }
    if ("bknots" %in% names(spline_params)) {
      if (length(spline_params$bknots) != 1) {
        stop(paste("All elements in 'bknots' in spline_params must have length",
                   "1 when mode is 'integrated'. Different spline parameters",
                   "for the different levels is not supported for this mode"),
             call. = FALSE)
      }
    }
  } else if (mode == "isolated") {
    num_levels <- length(unique(meta[[condition]]))
    if (any(sapply(spline_params, length) != num_levels)) {
      stop(paste("Each vector or list in spline_params must have as many", 
                 "elements as there are unique elements in the ", condition,
                 "column of meta when mode is 'isolated'."),
           call. = FALSE)
    }
    if ("knots" %in% names(spline_params)) {
      if (length(spline_params$knots) != num_levels) {
        stop(paste("'knots' in spline_params must have the same number of",
                   "elements as there are unique elements in the ", condition,
                   "column of meta when mode is 'isolated'."),
             call. = FALSE)
      }
    }
    if ("bknots" %in% names(spline_params)) {
      if (length(spline_params$bknots) != num_levels) {
        stop(paste("'bknots' in spline_params must have the same number of", 
                   "elements as there are unique elements in the ", condition,
                   "column of meta when mode is 'isolated'."),
             call. = FALSE)
      }
    }
  }
}


# Level 3 internal functions ---------------------------------------------------


#' Create Error Message
#'
#' @description
#' This function creates a formatted error message that includes the index of 
#' the data/meta pair if provided. 
#' If no index is provided, it returns the message as is.
#'
#' @param message A character string specifying the error message.
#' @param data_meta_index An optional parameter specifying the index of the 
#' data/meta pair for the error message. Default is NA.
#'
#' @return Returns a formatted error message string. If an index is provided, 
#' the message includes the index; otherwise, it returns the message as is.
#'
#' @examples
#' # Example of how to use the function
#' message <- "data must be a matrix."
#' create_error_message(message, 1)
#' create_error_message(message)
#'
#' @export
#' 
create_error_message <- function(message,
                                 data_meta_index = NA) {
  
  if (!is.na(data_meta_index)) {
    return(sprintf("data/meta pair index %d: %s", data_meta_index, message))
  } else {
    return(message)
  }
}
