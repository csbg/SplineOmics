#' extract_data.R contains the exported package function extract_data. This 
#' function automatically recognises the data field in a table and returns the
#' data matrix, that serves as input for the other functions of this package. 
#' This is for convenience only.


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
extract_data <- function(data, feature_name_column) {
  
  if (!is.data.frame(data)) {
    stop("Input data must be a dataframe.", call. = FALSE)
  }
  
  if (!is.character(feature_name_column) || length(feature_name_column) != 1) {
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
  
  if (nrow(data) == 0) {
    stop("Input dataframe is empty.", call. = FALSE)
  }
  
  # Convert data to character to handle mixed data types
  data_char <- as.data.frame(lapply(data, as.character), 
                             stringsAsFactors = FALSE)

  # Identify the starting point of the numeric data block
  upper_left_row <- NA
  upper_left_col <- NA
  num_rows <- nrow(data_char)
  num_cols <- ncol(data_char)
  
  for (i in 1:(num_rows - 5)) {
    for (j in 1:(num_cols - 5)) {
      block <- data_char[i:(i+5), j:(j+5)]
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
  
  # Expand the block vertically
  lower_right_row <- upper_left_row
  for (i in (upper_left_row+1):num_rows) {
    if (is.na(data_char[i, upper_left_col])) {
      break
    }
    lower_right_row <- i
  }
  
  # Expand the block horizontally
  lower_right_col <- upper_left_col
  for (j in (upper_left_col+1):num_cols) {
    if (is.na(data_char[upper_left_row, j])) {
      break
    }
    lower_right_col <- j
  }

  # Extract the numeric data block
  numeric_data <- data_char[upper_left_row:lower_right_row, 
                            upper_left_col:lower_right_col]
  
  numeric_data <- suppressWarnings(as.data.frame(lapply(numeric_data, 
                                                        as.numeric)))
  
  # Remove rows and columns that are entirely NA
  numeric_data <- numeric_data[rowSums(is.na(numeric_data)) != 
                                 ncol(numeric_data), ]
  numeric_data <- numeric_data[, colSums(is.na(numeric_data)) != 
                                 nrow(numeric_data)]
  
  # Remove rows with any NA values
  clean_data <- numeric_data[stats::complete.cases(numeric_data), ]
  
  # Extract headers for each column above the identified block
  headers <- sapply(upper_left_col:lower_right_col, function(col_idx) {
    header_values <- data_char[1:(upper_left_row-1), col_idx]
    header_values <- header_values[!is.na(header_values)]
    paste(header_values, collapse = "_")
  })
  
  colnames(clean_data) <- headers
  
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
  
  clean_data
}
