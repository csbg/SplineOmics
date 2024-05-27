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
                         feature_name_column = NA) {
  
  control_inputs_extract_data(data = data,
                              feature_name_column = feature_name_column)
  
  data <- as.data.frame(data)

  # Identify the starting point of the numeric data block
  upper_left_cell <- find_upper_left_cell(data = data)
  upper_left_row <- upper_left_cell$upper_left_row
  upper_left_col <- upper_left_cell$upper_left_col
  

  lower_right_cell <- find_lower_right_cell(data = data,
                                            upper_left_row = upper_left_row,
                                            upper_left_col = upper_left_col)
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
                                  feature_name_column = feature_name_column)
}



# Level 1 internal functions ---------------------------------------------------


#' Control Inputs for Extracting Data
#'
#' @description
#' This function checks the validity of input data and the feature name column.
#' It ensures that the input data is a dataframe, the feature name column is 
#' specified correctly, and contains valid data.
#'
#' @param data A dataframe containing the input data.
#' @param feature_name_column A string specifying the name of the feature 
#'                            column.
#'                            It must be a single string that matches a column 
#'                            name in the data. If `NA`, no column is checked.
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
                                        feature_name_column) {
  
  if (!is.data.frame(data)) {
    stop("Input data must be a dataframe.", call. = FALSE)
  }
  
  if (!is.na(feature_name_column)) {
    
    if (!is.character(feature_name_column) || 
        length(feature_name_column) != 1) {
      stop("Feature name column must be a single string.", call. = FALSE)
    }
    
    if (!(feature_name_column %in% colnames(data))) {
      stop(paste("Column", feature_name_column, "not found in the data."), 
           call. = FALSE)
    }
    
    if (all(is.na(data[[feature_name_column]]))) {
      stop(paste("Column '", feature_name_column,
                 "' contains only NA values.", sep = ""),
           call. = FALSE)
    }
  }
  
  if (nrow(data) == 0) {
    stop("Input dataframe is empty.", call. = FALSE)
  }
}


#' Find Upper Left Cell of Numeric Block
#'
#' @description
#' This function identifies the upper-left cell of the first 6x6 block 
#' of numeric values in a given dataframe.
#'
#' @param data A dataframe containing the input data.
#'
#' @details
#' The function searches for the first 6x6 block within the dataframe where 
#' all values are numeric. It iterates through the dataframe and checks each 
#' possible 6x6 block. If a valid block is found, it returns the row and 
#' column indices of the upper-left cell of the block.
#'
#' If no such block is found, the function stops with an error message.
#'
#' @return A list containing the row and column indices of the upper-left cell 
#'         of the identified numeric block. The list has two elements:
#'         \item{upper_left_row}{The row index of the upper-left cell.}
#'         \item{upper_left_col}{The column index of the upper-left cell.}
#'
#' @throws An error if no 6x6 block of numeric values is found.
#' 
find_upper_left_cell <- function(data) {
  
  upper_left_row <- NA
  upper_left_col <- NA
  num_rows <- nrow(data)
  num_cols <- ncol(data)
  
  for (i in 1:(num_rows - 5)) {
    for (j in 1:(num_cols - 5)) {
      block <- data[i:(i+5), j:(j+5)]
      block_num <- suppressWarnings(as.numeric(as.matrix(block)))
      if (all(!is.na(block_num)) && (all(is.numeric((block_num))))) {
        upper_left_row <- i
        upper_left_col <- j
        break
      }
    }
    if (!is.na(upper_left_row)) break
  }
  
  if (is.na(upper_left_row) || is.na(upper_left_col)) {
    stop("No at least 6x6 block of numeric values found.", call. = FALSE)
  }
  
  list(upper_left_row = upper_left_row,
       upper_left_col = upper_left_col)
}


#' Find Lower Right Cell of a Block
#'
#' @description
#' This function identifies the lower-right cell of a block of contiguous 
#' non-NA values starting from a given upper-left cell in a dataframe.
#'
#' @param data A dataframe containing the input data.
#' @param upper_left_row An integer specifying the row index of the upper-left 
#'                       cell of the block.
#' @param upper_left_col An integer specifying the column index of the 
#'                       upper-left 
#'                       cell of the block.
#'
#' @details
#' The function expands the block of contiguous non-NA values vertically and 
#' horizontally from the specified upper-left cell to identify the lower-right 
#' cell.
#' It iterates through the dataframe, expanding the block until it encounters NA 
#' values.
#'
#' @return A list containing the row and column indices of the lower-right cell 
#'         of the identified block. The list has two elements:
#'         \item{lower_right_row}{The row index of the lower-right cell.}
#'         \item{lower_right_col}{The column index of the lower-right cell.}
#'         
find_lower_right_cell <- function(data,
                                  upper_left_row,
                                  upper_left_col) {
  
  num_rows <- nrow(data)
  num_cols <- ncol(data)
  
  # Expand the block vertically
  lower_right_row <- upper_left_row
  for (i in (upper_left_row+1):num_rows) {
    if (is.na(data[i, upper_left_col])) {
      break
    }
    lower_right_row <- i
  }
  
  # Expand the block horizontally
  lower_right_col <- upper_left_col
  for (j in (upper_left_col+1):num_cols) {
    if (is.na(data[upper_left_row, j])) {
      break
    }
    lower_right_col <- j
  }
  
  list(lower_right_row = lower_right_row,
       lower_right_col = lower_right_col)
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
#' @param feature_name_column A string specifying the name of the feature column 
#'                            in `data`. If `NA`, sequential numbers will be
#'                             used 
#'                            as feature names.
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
#' @throws An error if:
#'         \item{feature names are not unique.}
#'         \item{the length of feature names does not match the number of rows 
#'         in `clean_data`.}
#'         
add_feature_names <- function(data,
                              clean_data,
                              feature_name_column) {
  
  if (!is.na(feature_name_column)) {
    
    feature_names <- as.character(data[[feature_name_column]])
    feature_names <- feature_names[!is.na(feature_names)]
    
    if (length(feature_names) != length(unique(feature_names))) {
      stop("Feature names in the column must be unique, ignoring NA values.",
           call. = FALSE)
    }
    
    if (length(feature_names) != nrow(clean_data)) {
      stop(paste("Length of feature names does not match the number of rows in",
                 "clean_data."), 
           call. = FALSE)
    }
    
    rownames(clean_data) <- feature_names
    
  } else {
    rownames(clean_data) <- as.character(seq_len(nrow(clean_data)))
    message(paste("No feature_name column specified. Setting numbers 1, 2, 3,", 
                  "etc. as the feature names"))
  }
  
  return(clean_data)
}
