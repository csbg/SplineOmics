#' utils scripts contains shared functions that are used by at least two package 
#' functions of the SplineOmics package.

# Level 1 internal functions ---------------------------------------------------


#' Create Progress Bar
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
    message = "Processing"
) {
  
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
    )   # Time column is mandatory
  
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
  } else {                                          # natural cubic splines
    meta$X <- do.call(splines::ns, args)
  }
  
  design_matrix <- stats::model.matrix(
    stats::as.formula(design),
    data = meta
    )
}


#' Determine Analysis Mode Based on Design Formula
#'
#' @description
#' This function determines whether each level should be analyzed in isolation 
#' or together based on the design formula. If the design formula includes 
#' interaction terms involving the factor of the experiment, the analysis 
#' mode is considered integrated (together). Otherwise, it is considered isolated.
#'
#' @param design A character string representing the design formula to be used 
#' for generating the model matrix.
#' @param factor_column A character string representing the column name of the 
#' factor of the experiment in the metadata.
#'
#' @return A character string indicating the analysis mode, either "integrated" 
#' if the design formula involves interaction terms with the factor of the 
#' experiment, or "isolated" otherwise.
#'
#' @keywords internal
#' 
determine_analysis_mode <- function(
    design, 
    factor_column
    ) {
  
  if (!is.character(design) || length(design) != 1) {
    stop(paste("design must be a character of length 1"),
         call. = FALSE)
  }
  
  if (!is.character(factor_column) || length(factor_column) != 1) {
    stop(paste("factor_column must be a character of length 1"),
         call. = FALSE)
  }
  
  # Extract terms from the design formula
  design_terms <- all.vars(stats::as.formula(design))
  
  # Check if the factor_column is part of any interaction term
  integrated <- any(grepl(paste0(factor_column, ":"), 
                          design) | grepl(paste0(":", factor_column), design))
  
  # If the factor_column is present without interactions, check if it's included
  if (!integrated && factor_column %in% design_terms) {
    # Check if the factor_column is part of any interaction term
    integrated <- any(grepl(paste0(factor_column, "*"), design))
  }
  
  # Return the mode
  if (integrated) {
    message(paste(
      "Interaction terms identified in limma design formula.",
      "Selecting mode == integrated. This means that all levels",
      "are analysed using the full data."
      ))
    return("integrated")
  } else {
    message(paste(
      "No interaction terms identified in limma design formula.",
      "Selecting mode == isolated. This means that every level",
      "is analysed using the level specific data."
      ))
    return("isolated")
  }
}


#' Merge Annotation with a Single Top Table
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
    annotation
) {
  
  top_table$feature_nr <- as.numeric(as.character(top_table$feature_nr))
  annotation_rows <- annotation[top_table$feature_nr, ]
  top_table <- cbind(top_table, annotation_rows)
}


#' Bind Data with Annotation
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
    annotation
) {
  
  data_df <- as.data.frame(data)
  
  # Check if the number of rows match (additional check, just for safety)
  if (nrow(data_df) != nrow(annotation)) {
    stop("The number of rows in data and annotation must be the same.",
         call. = FALSE)
  }
  
  # Add row names as the first column named feature_names
  combined_df <- cbind(
    feature_name = rownames(data_df),
    data_df, annotation
  )
}


#' Print Informational Message
#'
#' @description
#' This function prints a nicely formatted informational message with a green "Info" label.
#'
#' @param message_prefix A custom message prefix to be displayed before the success message.
#' @param report_dir The directory where the HTML reports are located.
#'
#' @return NULL
#'
print_info_message <- function(
    message_prefix,
    report_dir
    ) {
  
  # Green color code for "Info"
  green_info <- "\033[32mInfo\033[0m"
  
  full_message <- paste(
    green_info, message_prefix, "completed successfully.\n",
    "Your HTML reports are located in the directory: ", report_dir, ".\n",
    "Please note that due to embedded files, the reports might be flagged as\n",
    "harmful by other software. Rest assured that they provide no harm.\n"
  )
  
  cat(full_message)
}
