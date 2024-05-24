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
extract_data <- function(data,
                         feature_name_column) {
  
  if (!(feature_name_column %in% colnames(data))) {
    stop(paste("Column", feature_name_column, "not found in the data."), 
         call. = FALSE)
  }
  
  # Convert data to character to handle mixed data types
  data_char <- as.data.frame(lapply(data, as.character), 
                             stringsAsFactors = FALSE)
  
  # Identify the starting point of the numeric data block
  start_row <- NA
  start_col <- NA
  num_rows <- nrow(data_char)
  num_cols <- ncol(data_char)
  
  for (i in 1:(num_rows - 5)) {
    for (j in 1:(num_cols - 5)) {
      block <- data_char[i:(i+5), j:(j+5)]
      block_num <- suppressWarnings(as.numeric(as.matrix(block)))
      if (all(!is.na(block_num))) {
        start_row <- i
        start_col <- j
        break
      }
    }
    if (!is.na(start_row)) break
  }
  
  if (is.na(start_row) || is.na(start_col)) {
    stop("No at least6x6 block of numeric values found.", call. = FALSE)
  }
  
  # Extract the numeric data block
  numeric_data <- data_char[start_row:num_rows, start_col:num_cols]
  numeric_data <- suppressWarnings(as.data.frame(lapply(numeric_data, 
                                                        as.numeric)))
  
  # Remove rows and columns that are entirely NA
  numeric_data <- numeric_data[rowSums(is.na(numeric_data)) != 
                                 ncol(numeric_data), ]
  numeric_data <- numeric_data[, colSums(is.na(numeric_data)) != 
                                 nrow(numeric_data)]
  
  # Remove rows with any NA values
  clean_data <- numeric_data[stats::complete.cases(numeric_data), ]
  
  feature_names <- as.character(data[[feature_name_column]])
  feature_names <- feature_names[!is.na(feature_names)]
  
  if (length(feature_names) != nrow(clean_data)) {
    stop(paste("Length of feature names does not match the number of rows in",
               "clean_data."), 
         call. = FALSE)
  }
  
  rownames(clean_data) <- feature_names
  
  clean_data
}