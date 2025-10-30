#' Extract a rectangular block from a dataframe and make it to a matrix
#'
#' @description
#' This function extracts a rectangular block from a dataframe using
#' user-specified
#' top/bottom row indices and left/right column identifiers (numeric or
#' Excel-style letters).
#' It ensures the block contains only numeric values or NAs, and returns a
#' cleaned matrix.
#'
#' @param data `data.frame`: A dataframe containing the full input, including 
#' annotation columns and the numeric block to extract.
#'
#' @param top_row `integer(1)`: Specifies the first (top) row of the numeric 
#' data block. Row indexing is 1-based.
#'
#' @param bottom_row `integer(1)`: Specifies the last (bottom) row of the 
#' numeric data block. Must be >= `top_row`.
#'
#' @param left_col `integer(1)` | `character(1)`: Column specifier for the 
#' left-most column of the data block. Can be either:
#'   - An integer index (e.g., 2 for the second column), or
#'   - A character string using Excel-style letters (e.g., "A", "AB").
#'
#' Column names (e.g., "age") are **not** allowed here. Only letters or numeric 
#' indices are accepted.
#'
#' @param right_col `integer(1)` | `character(1)`: Same format as `left_col`. 
#' Specifies the right-most column of the numeric block. Must be >= `left_col` 
#' after conversion.
#'
#' @param feature_name_columns `character()` | `NULL`: Optional character 
#' vector specifying columns in `data` to be used as row (feature) names in 
#' the output. If `NA`, generic feature names are used. These row names are 
#' used everywhere to label the features, such as to label the plots in the 
#' `cluster_hits()` function report.
#'
#' @param use_row_index `logical(1)`: If \code{TRUE}, prepend the row index to 
#' each combined feature name to ensure uniqueness. Defaults to \code{FALSE}.
#'
#' @return A numeric matrix with cleaned values and appropriate column names.
#'
#' @examples
#' # Tiny demo table with two header rows, feature columns, and numeric block
#' df <- data.frame(
#'     feat_id = c(NA, NA, "g1", "g2", "g3"),
#'     feat_sym = c(NA, NA, "TP53", "EGFR", "BAX"),
#'     A = c("cond", "t0", 1, 2, 3),
#'     B = c("cond", "t1", 4, 5, 6),
#'     C = c("ctrl", "t0", 7, 8, 9),
#'     D = c("ctrl", "t1", 10, 11, 12),
#'     check.names = FALSE
#' )
#'
#' # Example 1: extract numeric block using Excel letters, build headers from
#' # the two rows above (they get collapsed like "cond_t0", "ctrl_t1", ...)
#' m1 <- extract_data(
#'     data = df,
#'     top_row = 3,
#'     bottom_row = 5,
#'     left_col = "A",
#'     right_col = "D",
#'     feature_name_columns = c("feat_id", "feat_sym"),
#'     use_row_index = FALSE
#' )
#' m1
#'
#' # Example 2: same extraction but with numeric column indices and row index
#' # prepended to ensure uniqueness of feature names
#' m2 <- extract_data(
#'     data = df,
#'     top_row = 3,
#'     bottom_row = 5,
#'     left_col = 3,
#'     right_col = 6,
#'     feature_name_columns = c("feat_id", "feat_sym"),
#'     use_row_index = TRUE
#' )
#' m2
#'
#' @export
#'
extract_data <- function(
    data,
    bottom_row,
    right_col,
    top_row = 1,
    left_col = 1,
    feature_name_columns = NA,
    use_row_index = FALSE) {
    control_inputs_extract_data(
        data = data,
        top_row = top_row,
        bottom_row = bottom_row,
        left_col = left_col,
        right_col = right_col,
        feature_name_columns = feature_name_columns,
        use_row_index = use_row_index
    )

    # Convert Excel-style letters to column indices if needed
    left_col_idx <- excel_col_to_index(left_col)
    right_col_idx <- excel_col_to_index(right_col)

    # Extract the data block
    data <- as.data.frame(data)
    numeric_data <- data[top_row:bottom_row, left_col_idx:right_col_idx]

    # Coerce to numeric
    numeric_data[] <- lapply(numeric_data, function(col) {
        as.numeric(as.character(col))
    })

    # Check that all values are numeric or NA
    valid_block <- all(vapply(numeric_data, function(col) {
        all(is.numeric(col) | is.na(col))
    }, logical(1)))

    if (!valid_block) {
        stop_call_false(
            "The selected block must contain only numeric values or NAs."
        )
    }

    headers <- vapply(left_col_idx:right_col_idx, function(col_idx) {
        if (top_row == 1) {
            # Just return existing column name
            return(colnames(data)[col_idx])
        }

        # Otherwise: build header from all non-NA entries above the data block
        header_values <- data[seq_len(top_row - 1), col_idx]
        header_values <- header_values[!is.na(header_values)]

        if (length(header_values) == 0) {
            return(paste0("X", col_idx)) # Fallback if nothing usable
        }

        paste(header_values, collapse = "_")
    }, character(1))

    colnames(numeric_data) <- headers

    numeric_data <- add_feature_names(
        data = data,
        data_matrix = numeric_data,
        feature_name_columns = feature_name_columns,
        top_row = top_row,
        bottom_row = bottom_row,
        use_row_index = use_row_index
    )

    clean_matrix <- as.matrix(numeric_data)
    rownames(clean_matrix) <- rownames(numeric_data)

    clean_matrix
}


# Level 1 internal functions ---------------------------------------------------


#' Control Inputs for Extracting Data
#'
#' @noRd
#'
#' @description
#' This function checks the validity of input data, user-defined row and column
#' boundaries, and the feature name columns. Ensures the selected block is valid
#' and the inputs are correctly specified.
#'
#' @param data A dataframe containing the input data.
#' @param top_row Integer specifying the top-most row of the numeric block.
#' @param bottom_row Integer specifying the bottom-most row of the numeric
#'  block.
#' @param left_col Column identifier (letter or numeric) for the left-most
#'  column.
#' @param right_col Column identifier (letter or numeric) for the right-most
#'  column.
#' @param feature_name_columns A character vector specifying the feature name
#' columns, or NA.
#' @param use_row_index Logical. If \code{TRUE}, prepend the row index to each
#'  combined feature name to ensure uniqueness. Defaults to \code{FALSE}.
#'
control_inputs_extract_data <- function(
    data,
    top_row,
    bottom_row,
    left_col,
    right_col,
    feature_name_columns,
    use_row_index) {
    if (!is.data.frame(data)) {
        stop_call_false(paste0(
            "Input 'data' must be a dataframe, but got: ",
            class(data)
        ))
    }

    if (nrow(data) == 0 || ncol(data) == 0) {
        stop_call_false(paste0(
            "Input dataframe is empty or has no columns (nrow = ",
            nrow(data),
            ", ncol = ",
            ncol(data), ")."
        ))
    }

    # Validate row inputs
    if (!is.numeric(top_row) || !is.numeric(bottom_row)) {
        stop_call_false(paste0(
            "'top_row' and 'bottom_row' must be numeric. Got: top_row = ",
            top_row,
            ", bottom_row = ",
            bottom_row
        ))
    }
    if (top_row <= 0 ||
        bottom_row <= 0 ||
        top_row != floor(top_row) ||
        bottom_row != floor(bottom_row)) {
        stop_call_false(paste0(
            "'top_row' and 'bottom_row' must be positive integers.",
            "Got: top_row = ",
            top_row,
            ", bottom_row = ",
            bottom_row
        ))
    }
    if (top_row > bottom_row) {
        stop_call_false(paste0(
            "'top_row' cannot be greater than 'bottom_row'. Got: top_row = ",
            top_row, ",
      bottom_row = ",
            bottom_row
        ))
    }
    if (bottom_row > nrow(data)) {
        stop_call_false(paste0(
            "'bottom_row' (",
            bottom_row,
            ") exceeds number of rows in data (nrow = ",
            nrow(data),
            ")."
        ))
    }

    left_col_idx <- excel_col_to_index(left_col)
    right_col_idx <- excel_col_to_index(right_col)

    if (any(is.na(c(left_col_idx, right_col_idx)))) {
        stop_call_false(paste0(
            "Invalid column reference(s). Got: left_col = '",
            left_col,
            "', right_col = '",
            right_col,
            "'"
        ))
    }

    if (left_col_idx <= 0 || right_col_idx <= 0) {
        stop_call_false(paste0(
            "Column indices must be positive integers. Got: left_col_idx = ",
            left_col_idx,
            ", right_col_idx = ",
            right_col_idx
        ))
    }

    if (left_col_idx > right_col_idx) {
        stop_call_false(paste0(
            "'left_col' index (",
            left_col_idx,
            ") cannot be after 'right_col' index (",
            right_col_idx,
            ")."
        ))
    }

    if (right_col_idx > ncol(data)) {
        stop_call_false(paste0(
            "'right_col' index (",
            right_col_idx,
            ") exceeds number of columns in data (ncol = ",
            ncol(data),
            ")."
        ))
    }

    # Validate feature name columns
    if (!any(is.na(feature_name_columns))) {
        if (!is.character(feature_name_columns)) {
            stop_call_false(paste0(
                "'feature_name_columns' must be a character vector. Got: ",
                class(feature_name_columns)
            ))
        }

        missing_columns <- setdiff(feature_name_columns, colnames(data))
        if (length(missing_columns) > 0) {
            stop_call_false(paste0(
                "The following 'feature_name_columns' are not present",
                "in the data: ",
                paste(missing_columns, collapse = ", ")
            ))
        }

        if (all(is.na(data[feature_name_columns]))) {
            stop_call_false(paste0(
                "Columns '", paste(feature_name_columns, collapse = ", "),
                "' contain only NA values."
            ))
        }
    }

    # Validate use_row_index argument
    if (!is.logical(use_row_index) ||
        length(use_row_index) != 1 ||
        is.na(use_row_index)) {
        stop_call_false(paste0(
            "'use_row_index' must be either TRUE or FALSE. Got: ",
            use_row_index
        ))
    }
}


#' Convert Excel-style column letters to numeric indices
#'
#' @description
#' Converts Excel-style column letters (e.g., "A", "Z", "AA", "AZ") into the
#' corresponding
#' numeric column indices. This is useful for allowing users to input column
#' references
#' in a more familiar format.
#'
#' @param col A character vector of column letters (e.g., "A", "BC") or numeric
#' column indices.
#'            If numeric, the input is returned unchanged (as integers).
#'
#' @return An integer vector of column indices corresponding to the input
#' letters.
#'
#' @examples
#' excel_col_to_index("A") # returns 1
#' excel_col_to_index("Z") # returns 26
#' excel_col_to_index("AA") # returns 27
#' excel_col_to_index("AZ") # returns 52
#' excel_col_to_index(c("A", "B", "Z", "AA", "AZ")) # returns
#' c(1, 2, 26, 27, 52)
#' excel_col_to_index(5) # returns 5 (numeric input is returned as-is)
#'
#' @noRd
#'
excel_col_to_index <- function(col) {
    if (is.numeric(col)) {
        return(as.integer(col))
    }
    col <- toupper(as.character(col))

    vapply(col, function(cname) {
        if (is.na(cname) || cname == "") {
            return(NA_integer_)
        }
        ch <- strsplit(cname, "", fixed = TRUE)[[1]]

        # Validate A–Z only
        if (any(ch < "A" | ch > "Z")) stop("Invalid column label: ", cname)

        nums <- vapply(ch, function(l) utf8ToInt(l) - 64L, integer(1))
        pows <- 26L^((length(nums) - 1L):0L)
        as.integer(sum(nums * pows))
    }, integer(1))
}


#' Add Feature Names to Data
#'
#' @noRd
#'
#' @description
#' This function assigns feature names to the rows of a dataframe based on a
#' specified column from another dataframe. If no column is specified, it
#' assigns sequential numbers as feature names.
#'
#' @param data A dataframe containing the original data with feature names.
#' @param data_matrix A dataframe to which the feature names will be added.
#' @param feature_name_columns A string specifying the name of the feature
#'                             columns in `data`. If `NA`, sequential numbers
#'                             will be used as feature names.
#' @param top_row Integer. Specifies the first (top) row of the numeric data
#' block. Row indexing is 1-based.
#' @param bottom_row Integer. Specifies the last (bottom) row of the numeric
#' data block. Must be >= `top_row`.
#' @param use_row_index Logical. If \code{TRUE}, prepend the row index to each
#'  combined feature name to ensure uniqueness. Defaults to \code{FALSE}.
#'
#' @details
#' The function performs the following operations:
#' - Extracts feature names from the specified column in `data`, ignoring
#'  `NA` values.
#' - Ensures the feature names are unique and match the number of rows in
#' `data_matrix`.
#' - Assigns the feature names to the rows of `data_matrix`.
#' - If `feature_name_column` is `NA`, assigns sequential numbers
#' (1, 2, 3, etc.)
#'   as feature names and issues a message.
#'
#' @return The `data_matrix` dataframe with updated row names.
#'
add_feature_names <- function(
    data,
    data_matrix,
    feature_name_columns,
    top_row,
    bottom_row,
    use_row_index = FALSE) {
    if (!any(is.na(feature_name_columns))) {
        rows <- seq.int(top_row, bottom_row)

        data_filtered <- data[rows, , drop = FALSE]

        feature_names <- apply(
            data_filtered[feature_name_columns],
            1,
            function(row) paste(row, collapse = "_")
        )

        feature_names <- as.character(feature_names)

        if (use_row_index) {
            row_indices <- seq_len(length(feature_names))
            feature_names <- paste(
                row_indices,
                feature_names,
                sep = "_"
            )
        }

        if (length(feature_names) != length(unique(feature_names))) {
            stop_call_false(
                "Combined feature names must be unique, ignoring NA values."
            )
        }

        if (length(feature_names) != nrow(data_matrix)) {
            stop_call_false(
                paste(
                    "Length of combined feature names does not match",
                    "the number of",
                    "rows in data_matrix"
                )
            )
        }

        rownames(data_matrix) <- feature_names
    } else {
        rownames(data_matrix) <- as.character(seq_len(nrow(data_matrix)))
        message(paste(
            "No feature_name column specified. Setting numbers 1, 2, 3,",
            "etc. as the feature names"
        ))
    }
    return(data_matrix)
}
