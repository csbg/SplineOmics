#' extract_data.R contains the exported package function extract_data. This 
#' function automatically recognises the data field in a table and returns the
#' data matrix, that serves as input for the other functions of this package. 
#' This is for convenience only.


# Exported function: extract_data ----------------------------------------------

#' Extract Numeric Matrix from Dataframe
#'
#' @description
#' This function takes a dataframe and identifies a rectangular or quadratic 
#' area containing numeric data, starting from the first occurrence of a 
#' 6x6 block of numeric values. It then extracts this area into a matrix, 
#' ensuring that each row contains only numeric values. Rows with any NA values 
#' are removed from the resulting matrix.
#'
#' @param data A dataframe loaded from a csv file, potentially containing a 
#' rectangular or quadratic area with numeric data amidst other values.
#' @param feature_name_columns (Optional) A character vector, specifying the
#'                             columns of the dataframe data, that should be 
#'                             used to construct the feature names. If ommited,
#'                             the feature names are just numbers (stored as
#'                             characters) starting from 1 (1, 2, 3, etc.)
#'
#' @return A matrix containing only the numeric data from the identified area, 
#' with rows containing NA values removed.
#'
#' @examples
#' \dontrun{
#' # Example of how to use the function
#' data <- read.csv("path_to_your_file.csv")
#' result_matrix <- extract_numeric_matrix(data)
#' print(result_matrix)}
#'
#' @importFrom stats complete.cases
#'
#' @export
#' 
extract_data <- function(data,
                         feature_name_columns = NA) {
  
  control_inputs_extract_data(data = data,
                              feature_name_columns = feature_name_columns)
  
  data <- as.data.frame(data)
  
  numeric_block_finder <- NumericBlockFinder$new(data)
  upper_left_cell <- numeric_block_finder$find_upper_left_cell()
  lower_right_cell <- numeric_block_finder$find_lower_right_cell()
  
  upper_left_row <- upper_left_cell$upper_left_row
  upper_left_col <- upper_left_cell$upper_left_col
  
  lower_right_row <- lower_right_cell$lower_right_row
  lower_right_col <- lower_right_cell$lower_right_col

  
  # Extract the numeric data block
  numeric_data <- data[upper_left_row:lower_right_row,
                       upper_left_col:lower_right_col]

  numeric_data[] <-
    lapply(numeric_data, function(col) suppressWarnings(
      as.numeric(as.character(col))))

  # Check if every element of numeric_data is numeric
  if (any(sapply(numeric_data, function(col) all(is.na(col))))) {
    stop(paste("All elements of the numeric data must be numeric. Please",
               "ensure there is an empty column between the numeric data and",
               "the annotation information, which, if present, must be on the",
               "right of the numeric data, not on the left."),
               call. = FALSE)
  }

  # Remove rows and columns that are entirely NA
  numeric_data <- numeric_data[rowSums(is.na(numeric_data)) !=
                                 ncol(numeric_data), ]
  numeric_data <- numeric_data[, colSums(is.na(numeric_data)) !=
                                 nrow(numeric_data)]

  # Remove rows with any NA values
  clean_data <- numeric_data[stats::complete.cases(numeric_data), ]
  
  # Extract headers for each column above the identified block
  headers <- sapply(upper_left_col:lower_right_col, function(col_idx) {
    header_values <- data[1:(upper_left_row-1), col_idx]
    header_values <- header_values[!is.na(header_values)]
    paste(header_values, collapse = "_")
  })
  
  colnames(clean_data) <- headers
  
  clean_data <- add_feature_names(data = data,
                                  clean_data = clean_data,
                                  feature_name_columns = feature_name_columns)
}



# Class NumericBlockFinder -----------------------------------------------------


#' NumericBlockFinder: A class for finding numeric blocks in data
#'
#' This class provides methods to identify the upper-left and lower-right 
#' cells of a numeric block within a dataframe.
#'
#' @field data A dataframe containing the input data.
#' @field upper_left_cell A list containing the row and column indices of the
#'                        upper-left cell.
#'
NumericBlockFinder <- R6::R6Class("NumericBlockFinder",
  public = list(
    data = NULL,
    upper_left_cell = NULL,
    
    
    #' Initialize a NumericBlockFinder object
    #'
    #' @param data A dataframe containing the input data.
    #' @return A new instance of the NumericBlockFinder class.
    #' 
    initialize = function(data) {
      
      self$data <- as.data.frame(data)
    },
    
    
    #' Find the upper-left cell of the first 6x6 block of numeric values
    #'
    #' This method identifies the upper-left cell of the first 6x6 block of 
    #' numeric values in the dataframe.
    #'
    #' @return A list containing the row and column indices of the upper-left 
    #'         cell.
    #'         
    find_upper_left_cell = function() {
      
      upper_left_row <- NA
      upper_left_col <- NA
      num_rows <- nrow(self$data)
      num_cols <- ncol(self$data)
      
      for (i in 1:(num_rows - 5)) {
        for (j in 1:(num_cols - 5)) {
          block <- self$data[i:(i+5), j:(j+5)]
          block_num <- suppressWarnings(as.numeric(as.matrix(block)))
          if (all(!is.na(block_num)) && (all(is.numeric(block_num)))) {
            upper_left_row <- i
            upper_left_col <- j
            break
          }
        }
        if (!is.na(upper_left_row)) break
      }
      
      if (is.na(upper_left_row) || is.na(upper_left_col)) {
        stop("No at least 6x6 block of numeric values found.",
             call. = FALSE)
      }
      
      self$upper_left_cell <- list(upper_left_row = upper_left_row,
                                   upper_left_col = upper_left_col)
      return(self$upper_left_cell)
    },
    
    
    #' Find the lower-right cell of a block of contiguous non-NA values
    #'
    #' This method identifies the lower-right cell of a block of contiguous
    #' non-NA values starting from a given upper-left cell in the dataframe.
    #'
    #' @return A list containing the row and column indices of the lower-right
    #'         cell.
    #' 
    find_lower_right_cell = function() {
      
      if (is.null(self$upper_left_cell)) {
        stop(paste("Upper-left cell has not been identified.", 
                   "Call find_upper_left_cell first."), call. = FALSE)
      }
      
      upper_left_row <- self$upper_left_cell$upper_left_row
      upper_left_col <- self$upper_left_cell$upper_left_col
      num_rows <- nrow(self$data)
      num_cols <- ncol(self$data)
      
      # Expand the block vertically
      lower_right_row <- upper_left_row
      for (i in (upper_left_row+1):num_rows) {
        if (is.na(self$data[i, upper_left_col])) {
          break
        }
        lower_right_row <- i
      }
      
      # Expand the block horizontally
      lower_right_col <- upper_left_col
      for (j in (upper_left_col+1):num_cols) {
        if (is.na(self$data[upper_left_row, j])) {
          break
        }
        lower_right_col <- j
      }
      
      list(lower_right_row = lower_right_row,
           lower_right_col = lower_right_col)
    }
  )
)


# Level 1 internal functions ---------------------------------------------------


#' Control Inputs for Extracting Data
#'
#' @description
#' This function checks the validity of input data and the feature name column.
#' It ensures that the input data is a dataframe, the feature name column is 
#' specified correctly, and contains valid data.
#'
#' @param data A dataframe containing the input data.
#' @param feature_name_columns A character vector specifying the names of the 
#'                             feature name columns. The columns must be present
#'                             in the dataframe data. If `NA`, no column is 
#'                             checked.
#'
#' @details
#' The function performs the following checks:
#' - Ensures the input data is a dataframe.
#' - Checks if the feature name column is a single string and exists in the 
#' data.
#' - Ensures the specified feature name column does not contain only `NA` 
#' values.
#' - Checks if the input dataframe is not empty.
#'
#' If any of these checks fail, the function stops with an error message.
#' 
control_inputs_extract_data <- function(data,
                                        feature_name_columns) {
  
  if (!is.data.frame(data)) {
    stop("Input data must be a dataframe.", call. = FALSE)
  }
  
  if (!any(is.na(feature_name_columns))) {
    
    if (!is.character(feature_name_columns)) {
      stop("feature_name_columns should be a character vector.", call. = FALSE)
    }
    
    if (!all(feature_name_columns %in% colnames(data))) {
      stop("Not all feature_name_columns are present in the data.", 
           call. = FALSE)
    }
    
    if (!any(is.na(feature_name_columns))) {
      if (all(is.na(data[feature_name_columns]))) {
        stop(paste("Columns '", paste(feature_name_columns, collapse = ", "),
                   "' contain only NA values.", sep = ""),
             call. = FALSE)
      }
    }
  }
  
  if (nrow(data) == 0) {
    stop("Input dataframe is empty.", call. = FALSE)
  }
}


#' Add Feature Names to Data
#'
#' @description
#' This function assigns feature names to the rows of a dataframe based on a 
#' specified column from another dataframe. If no column is specified, it 
#' assigns sequential numbers as feature names.
#'
#' @param data A dataframe containing the original data with feature names.
#' @param clean_data A dataframe to which the feature names will be added.
#' @param feature_name_columns A string specifying the name of the feature 
#'                             columns in `data`. If `NA`, sequential numbers 
#'                             will be used as feature names.
#'
#' @details
#' The function performs the following operations:
#' - Extracts feature names from the specified column in `data`, ignoring
#'  `NA` values.
#' - Ensures the feature names are unique and match the number of rows in 
#' `clean_data`.
#' - Assigns the feature names to the rows of `clean_data`.
#' - If `feature_name_column` is `NA`, assigns sequential numbers 
#' (1, 2, 3, etc.) 
#'   as feature names and issues a message.
#'
#' @return The `clean_data` dataframe with updated row names.
#'         
add_feature_names <- function(data, 
                              clean_data, 
                              feature_name_columns) {
  
  if (!any(is.na(feature_name_columns))) {
    
    non_na_index <- which(apply(data[feature_name_columns], 1, 
                                function(row) all(!is.na(row))))[1]
    data_filtered <- data[non_na_index:nrow(data), , drop = FALSE]
    clean_data_filtered <- clean_data[non_na_index:nrow(clean_data), ,
                                      drop = FALSE]
    
    # Extract and combine the feature names from specified columns
    feature_names <- apply(data_filtered[feature_name_columns], 1, 
                           function(row) paste(row, collapse = "_"))
    
    # Check for NA values in combined feature names
    feature_names <- as.character(feature_names)
    feature_names <- feature_names[!is.na(feature_names)]

    # Ensure unique feature names
    if (length(feature_names) != length(unique(feature_names))) {
      stop("Combined feature names must be unique, ignoring NA values.", 
           call. = FALSE)
    }
    
    # Ensure the length matches the number of rows in clean_data
    if (length(feature_names) != nrow(clean_data)) {
      stop(paste("Length of combined feature names does not match the number of",
                 "rows in clean_data."), 
           call. = FALSE)
    }
    
    # Assign combined feature names as row names
    rownames(clean_data) <- feature_names
    
  } else {
    
    rownames(clean_data) <- as.character(seq_len(nrow(clean_data)))
    message(paste("No feature_name column specified. Setting numbers 1, 2, 3,", 
                  "etc. as the feature names"))
  }
  return(clean_data)
}
