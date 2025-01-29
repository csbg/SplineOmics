#' utils scripts contains shared functions that are used by at least two package
#' functions of the SplineOmics package.

# Level 1 internal functions ---------------------------------------------------


#' Create Progress Bar
#' 
#' @noRd
#'
#' @description
#' Creates a progress bar for tracking the progress of an iterable task.
#'
#' @param iterable An iterable object (e.g., list or vector) whose length
#' determines the total number of steps.
#' @param message A message to display with the progress bar
#' (default is "Processing").
#'
#' @return A progress bar object from the 'progress' package.
#'
#' @importFrom progress progress_bar
#'
#' @seealso
#' \code{\link{progress_bar}}
#'
create_progress_bar <- function(
    iterable,
    message = "Processing") {
  # Create and return the progress bar
  pb <- progress::progress_bar$new(
    format = paste(" ", message, " [:bar] :percent :elapsed"),
    total = length(iterable),
    width = 60,
    clear = FALSE
  )

  return(pb)
}



#' Create Design Matrix for Splines
#' 
#' @noRd
#'
#' @description
#' This function generates a design matrix using spline parameters and metadata.
#' It accommodates both B-splines and natural cubic splines based on the provided
#' spline type and parameters.
#'
#' @param meta A dataframe containing the metadata, including the time column.
#' @param spline_params A list containing the spline parameters. This list can
#' include `dof` (degrees of freedom), `knots`, `bknots` (boundary knots),
#' `spline_type`, and `degree`.
#' @param level_index An integer representing the current level index for which
#' the design matrix is being generated.
#' @param design A character string representing the design formula to be used
#' for generating the model matrix.
#'
#' @return A design matrix constructed using the specified spline parameters and
#' design formula.
#'
#' @importFrom splines bs ns
#' @importFrom stats model.matrix as.formula
#'
design2design_matrix <- function(
    meta,
    spline_params,
    level_index,
    design
    ) {
  
  args <- list(
    x = meta$Time,
    intercept = FALSE
  ) # Time column is mandatory

  if (!is.null(spline_params$dof)) {
    args$df <- spline_params$dof[level_index]
  } else {
    args$knots <- spline_params$knots[[level_index]]
  }

  if (!is.null(spline_params$bknots)) {
    args$Boundary.knots <- spline_params$bknots[[level_index]]
  }

  if (spline_params$spline_type[level_index] == "b") {
    args$degree <- spline_params$degree[level_index]
    meta$X <- do.call(splines::bs, args)
  } else { # natural cubic splines
    meta$X <- do.call(splines::ns, args)
  }

  design_matrix <- stats::model.matrix(
    stats::as.formula(design),
    data = meta
  )
  
  result <- list(
    design_matrix = design_matrix,
    meta = meta
  )
}


#' Merge Annotation with a Single Top Table
#' 
#' @noRd
#'
#' @description
#' This function merges annotation information into a single `top_table`
#' dataframe based on the `feature_nr` column.
#'
#' @param top_table A dataframe containing the `top_table` with a `feature_nr` column.
#' @param annotation A dataframe containing the annotation information.
#'
#' @return A dataframe with updated `top_table` containing merged annotation information.
#'
merge_top_table_with_annotation <- function(
    top_table,
    annotation) {
  top_table$feature_nr <- as.numeric(as.character(top_table$feature_nr))
  annotation_rows <- annotation[top_table$feature_nr, ]
  top_table <- cbind(top_table, annotation_rows)
}


#' Bind Data with Annotation
#' 
#' @noRd
#'
#' @description
#' This function converts a matrix to a dataframe, adds row names as the first
#' column,
#' and binds it with annotation data.
#'
#' @param data A matrix containing the numeric data.
#' @param annotation A dataframe containing the annotation information.
#'
#' @return A dataframe with `data` and `annotation` combined, and the row names
#'  of `data`
#' as the first column named `feature_names`.
#'
bind_data_with_annotation <- function(
    data,
    annotation = NULL) {
  data_df <- as.data.frame(data)

  # Add row names as the first column named feature_names
  combined_df <- cbind(
    feature_name = rownames(data_df),
    data_df
  )

  # If annotation is not NULL, check row count and bind with annotation
  if (!is.null(annotation)) {
    if (nrow(data_df) != nrow(annotation)) {
      stop("The number of rows in data and annotation must be the same.",
        call. = FALSE
      )
    }

    # Bind the annotation with the data
    combined_df <- cbind(combined_df, annotation)
  }

  return(combined_df)
}


#' Print Informational Message
#' 
#' @noRd
#'
#' @description
#' This function prints a nicely formatted informational message with a green
#'  "Info" label.
#'
#' @param message_prefix A custom message prefix to be displayed before the
#' success message.
#' @param report_dir The directory where the HTML reports are located.
#'
#' @return NULL
#'
print_info_message <- function(
    message_prefix,
    report_dir) {
  # Green color code for "Info"
  green_info <- "\033[32mInfo\033[0m"

  full_message <- paste(
    green_info, message_prefix, "completed successfully.\n",
    "Your HTML reports are located in the directory: ", report_dir, ".\n",
    "Please note that due to embedded files, the reports might be flagged as\n",
    "harmful by other software. Rest assured that they provide no harm.\n"
  )

  message(full_message)
}


#' Stop with custom message without call.
#' 
#' @noRd
#'
#' @description
#' A helper function that triggers an error with the specified message and
#' suppresses the function call in the error output. This function behaves
#' similarly to the base `stop()` function but automatically concatenates
#' multiple message strings if provided.
#'
#' @param ... One or more character strings specifying the error message.
#'            If multiple strings are provided, they will be concatenated
#'            with a space between them.
#'
#' @return This function does not return a value; it stops execution and
#'         throws an error.
#'
stop_call_false <- function(...) {
  # Concatenate all arguments into a single string
  message_text <- paste(..., sep = " ")

  # Call stop with the concatenated message and call. = FALSE
  stop(message_text, call. = FALSE)
}


#' Extract fixed and random effects from a model formula string.
#'
#' @noRd
#'
#' @description
#' This function processes a model formula string by separating it into 
#' fixed effects and random effects. Random effects are substrings in the
#' format "(...)", and fixed effects are the remaining part of the string.
#'
#' @param formula_string A character string representing the model formula.
#'
#' @return A list containing two components:
#'         - `fixed_effects`: The cleaned fixed effects string.
#'         - `random_effects`: The concatenated random effects string.
#'
#' @examples
#' extract_effects("~ 1 + Condition*Time + Plate + (1|Reactor)")
#' # Returns:
#' # $fixed_effects
#' # [1] "~ 1 + Condition*Time + Plate"
#' #
#' # $random_effects
#' # [1] "(1|Reactor)"
extract_effects <- function(formula_string) {
  # Match all substrings in the format "(...)"
  matches <- gregexpr("\\([^)]*\\)", formula_string, perl = TRUE)
  substrings <- regmatches(formula_string, matches)[[1]]
  
  # If there are random effects, clean them
  if (length(substrings) > 0) {
    # Remove one non-whitespace character before and after each match
    cleaned_random_effects <- sapply(substrings, function(substring) {
      substring <- sub("\\S\\(", "(", substring)
      substring <- sub("\\)\\S", ")", substring)
      return(substring)
    }, USE.NAMES = FALSE)
    
    # Concatenate random effects into a single string
    random_effects <- paste(cleaned_random_effects, collapse = " ")
  } else {
    random_effects <- ""
  }
  
  # Remove all random effects from the original string to get fixed effects
  fixed_effects <- formula_string
  if (length(substrings) > 0) {
    for (substring in substrings) {
      fixed_effects <- sub(substring, "", fixed_effects, fixed = TRUE)
    }
  }
  
  # Clean up any extra '+' or whitespace from fixed effects
  fixed_effects <- gsub("\\s*\\+\\s*$", "", fixed_effects)
  fixed_effects <- trimws(fixed_effects)
  
  # Return both fixed and random effects
  return(list(
    fixed_effects = fixed_effects,
    random_effects = random_effects
  ))
}

