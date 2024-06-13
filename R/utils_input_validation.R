

# InputControl class -----------------------------------------------------------


#' InputControl: A class for controlling and validating inputs
#'
#' This class provides methods to validate the inputs of a function.
#'
#' @field args A list of arguments to be validated.
#'
InputControl <- R6::R6Class("InputControl",
  inherit = Level2Functions,
  
  public = list(
    args = NULL,
    
    #' Initialize an InputControl object
    #'
    #' @param args A list of arguments to be validated.
    #' @return A new instance of the InputControl class.
    #' 
    initialize = function(args) {
      self$args <- args
    },
    
    
    #' Automatically Validate All Arguments
    #'
    #' This method automatically validates all arguments by sequentially 
    #' calling
    #' various validation methods defined within the class. Each validation 
    #' method
    #' checks specific aspects of the input arguments and raises an error if the
    #' validation fails.
    #'
    #' The following validation methods are called in sequence:
    #' - \code{self$check_data_and_meta()}
    #' - \code{self$check_datas_and_metas()}
    #' - \code{self$check_datas_descr()}
    #' - \code{self$check_mode()}
    #' - \code{self$check_modes()}
    #' - \code{self$check_design_formula()}
    #' - \code{self$check_designs_and_metas()}
    #' - \code{self$check_spline_params()}
    #' - \code{self$check_spline_test_configs()}
    #' - \code{self$check_adj_pthresholds()}
    #' - \code{self$check_clusters()}
    #' - \code{self$check_time_unit()}
    #' - \code{self$check_report_dir()}
    #' - \code{self$check_padjust_method()}
    #' - \code{self$check_report_info()}
    #' - \code{self$check_report()}
    #'
    #' @return NULL. The function is used for its side effects of validating 
    #' input
    #' arguments and raising errors if any validation fails.
    #'
    auto_validate = function() {
      self$check_data_and_meta()
      self$check_datas_and_metas()
      self$check_datas_descr()
      self$check_top_tables()
      self$check_mode()
      self$check_modes()
      self$check_design_formula()
      self$check_designs_and_metas()
      self$check_spline_params()
      self$check_spline_test_configs()
      self$check_adj_pthresholds()
      self$check_clusters()
      self$check_time_unit()
      self$check_report_dir()
      self$check_padjust_method()
      self$check_report_info()
      self$check_report()
    },
    

    #' Check Data and Meta
    #'
    #' @description
    #' This function checks the validity of the data and meta objects, 
    #' ensuring that
    #' data is a matrix with numeric values and that meta is a dataframe 
    #' containing
    #' the specified condition column. Additionally, it verifies that the 
    #' number of
    #' columns in the data matrix matches the number of rows in the meta 
    #' dataframe.
    #'
    #' @param data A matrix containing numeric values.
    #' @param meta A dataframe containing the metadata, including the 'Time' 
    #' column
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
    check_data_and_meta = function() {
      
      data <- self$args[["data"]]
      meta <- self$args[["meta"]]
      condition <- self$args[["condition"]]
      meta_batch_column <- self$args[["meta_batch_column"]]
      meta_batch2_column <- self$args[["meta_batch2_column"]]
      data_meta_index <- self$args[["data_meta_index"]]

      required_args <- list(data, meta, condition)
      
      if (any(sapply(required_args, is.null))) {
        return(NULL)
      }
      
      self$check_data(data,
                      data_meta_index)
      
      self$check_meta(meta = meta,
                      condition = condition,
                      meta_batch_column = meta_batch_column,
                      meta_batch2_column = meta_batch2_column,
                      data_meta_index = data_meta_index)
      
      if (!(nrow(meta) == ncol(data))) {
        if (!is.null(data_meta_index)) {
          stop(paste0("For index ", data_meta_index,
                      "data column number must be equal to ",
                      "meta row number"),
               call. = FALSE)
        } else {
          stop(paste0("data column number must be equal to meta row number"),
               call. = FALSE)
        }
      }
    },
    
    
    #' Check Multiple Data and Meta Pairs
    #'
    #' @description
    #' Iterates over multiple data and meta pairs to validate each pair using 
    #' the `check_data_and_meta` function.
    #'
    #' @param datas A list of matrices containing numeric values.
    #' @param metas A list of data frames containing metadata.
    #' @param condition A character string specifying the column name in the 
    #' meta dataframe to be checked.
    #' @param meta_batch_column An optional parameter specifying the column name 
    #' in the meta dataframe used to remove the batch effect. Default is NA.
    #' @param meta_batch2_column An optional parameter specifying the column
    #'  name 
    #' in the meta dataframe used to remove the second batch effect. Default 
    #' is NA.
    #'
    #' @return NULL if any check fails, otherwise returns TRUE.
    #'
    check_datas_and_metas = function() {
      
      datas <- self$args$datas
      metas <- self$args$metas
      condition <- self$args$condition
      meta_batch_column <- self$args$meta_batch_column
      meta_batch2_column <- self$args$meta_batch2_column
      
      required_args <- list(datas, metas, condition)
      
      if (any(sapply(required_args, is.null))) {
        return(NULL)
      }
      
      if (length(datas) != length(metas)) {
        stop(paste0("datas and metas must have the same length."),
             call. = FALSE)
      }
      
      data_storage <- self$args$data
      meta_storage <- self$args$meta
      
      for (i in seq_along(datas)) {
        self$args$data <- datas[[i]]
        self$args$meta <- metas[[i]]
        self$args$data_meta_index <- i
        self$args$condition <- condition
        self$args$meta_batch_column <- meta_batch_column
        self$args$meta_batch2_column <- meta_batch2_column
        
        self$check_data_and_meta()
      }
      
      self$args$data <- data_storage
      self$args$meta <- meta_storage
      
      return(TRUE)
    },
    
    
    #' Check Data Descriptions
    #'
    #' @description
    #' Validates that the data descriptions are character vectors with each 
    #' element 
    #' not exceeding 80 characters in length.
    #'
    #' @param datas_descr A character vector of data descriptions.
    #'
    #' @return No return value, called for side effects.
    #'
    #' @seealso
    #' \code{\link{stop}} for error handling.
    #' 
    check_datas_descr = function() {
      
      datas_descr <- self$args$datas_descr
      
      required_args <- list(datas_descr)
      
      if (any(sapply(required_args, is.null))) {
        return(NULL)
      }
      
      if (!is.character(datas_descr) || any(nchar(datas_descr) > 80)) {
        long_elements_indices <- which(nchar(datas_descr) > 80)
        long_elements <- datas_descr[long_elements_indices]
        error_message <- sprintf(
          "'datas_descr' must be a character vector with no element over 80 
      characters. Offending element(s) at indices %s: '%s'. Please shorten
      the description.",
          paste(long_elements_indices, collapse = ", "),
          paste(long_elements, collapse = "', '")
        )
        stop(error_message)
      }
    },
    
    
    #' Check Top Tables
    #'
    #' @description
    #' Validates that the top tables are a list of dataframes and checks each
    #' dataframe using the `check_dataframe` function.
    #'
    #' @param top_tables A list of top tables from limma analysis.
    #'
    #' @return No return value, called for side effects.
    #'
    check_top_tables = function() {
      
      top_tables <- self$args$top_tables

      required_args <- list(top_tables)
      
      if (any(sapply(required_args, is.null))) {
        return(NULL)
      }
      
      # Helper function to check data frames in a list
      check_list_of_dataframes <- function(df_list) {
        
        if (!is.list(df_list) || !all(sapply(df_list, is.data.frame))) {
          
          stop("Expected a list of dataframes", call. = FALSE)
          
        } else {
          
          if (length(df_list) != 2) {
            stop(paste("top_tables must be a list of two lists if you want", 
                       "to cluster the hits of the limma results of", 
                       "avrg_diff_conditions or interaction_condition_time.", 
                       "The list must contain one of those two plus", 
                       "the time_effect limma result, in any order"),
                 call. = FALSE)
          }
          
          underscore_count <- sapply(names(df_list), 
                                     function(name) sum(grepl("_", name)))
          if (sum(underscore_count == 1) != 1 ||
              sum(underscore_count == 4) != 1) {
            stop(paste("top_tables must be a list of two lists if you want", 
                       "to cluster the hits of the limma results of", 
                       "avrg_diff_conditions or interaction_condition_time.", 
                       "The list must contain one of those two plus", 
                       "the time_effect limma result, in any order"),
                 call. = FALSE)
          }
          
          for (df in df_list) {
            self$check_dataframe(df)
          }
        }
      }
      
      # Check if top_tables is a list
      if (!is.list(top_tables)) {
        stop("top_tables must be a list", call. = FALSE)
      }
      
      for (i in seq_along(top_tables)) {
        
        element <- top_tables[[i]]
        element_name <- names(top_tables)[i]
        
        # Means it is a list of lists (so the user wants to cluster the hits of 
        # avrg_diff_conditions or interaction_condition_time with the help of 
        # the spline coeffs in time_effect)
        if (!tibble::is_tibble(element) && is.list(element)) {
          
          check_list_of_dataframes(element)
          
        } else if (tibble::is_tibble(element)) {

          matches <- gregexpr("_", string)
          underscore_count <- sum(sapply(matches,
                              function(x) if (x[1] == -1) 0 else length(x)))
          if (underscore_count != 1) {
            stop(paste("Very likely you did not pass the time_effect result", 
                 "of limma but rather one of avrg_diff_conditions or", 
                 "interaction_condition_time, which cannot be passed alone", 
                 "(to cluster the hits of one of those, put only one of them", 
                 "at a time into a list with the time_effect result. Further", 
                 "note: Please do not edit those list element names by hand."),
                 call. = FALSE)
          }
          
          self$check_dataframe(element)
          
        } else {
          stop(paste("top_tables must contain either data frames or lists of", 
                     "data frames"),
               call. = FALSE)
        }
      }
    },
    
    
    #' Check Mode
    #'
    #' @description
    #' Validates that the mode is either 'integrated' or 'isolated', which 
    #' depends
    #' on the design formula used in limma.
    #'
    #' @param mode A character string specifying the mode.
    #'
    #' @return A message indicating the chosen mode if valid; otherwise, an 
    #' error
    #' is thrown.
    #'
    #' @seealso
    #' \code{\link{limma}}
    #'
    check_mode = function() {

      mode <- self$args[["mode"]]
      
      required_args <- list(mode)
      
      if (any(sapply(required_args, is.null))) {
        return(NULL)
      }

      if (!((mode == "integrated") || (mode == "isolated"))) {
        stop(paste("mode must be either integrated or isolated. This is", 
                   "dependent on the used limma design formula. For example,", 
                   "this formula: ~ 1 + Phase*X + Reactor would require", 
                   "mode = integrated, whereas this formula:", 
                   "~ 1 + X + Reactor would require mode = isolated."),
             call. = FALSE)
      } else {
        sprintf("Mode %s chosen.", mode)
      }
    },
    
    
    #' Check Modes
    #'
    #' This method validates multiple modes by ensuring that the lengths of the 
    #' modes and designs are the same. It then iterates through each mode, 
    #' sets it 
    #' in the class arguments, and calls \code{self$check_mode()} to perform 
    #' the validation for each mode.
    #'
    #' @details
    #' - The method first checks if the lengths of \code{designs} and 
    #' \code{modes} 
    #'   are equal. If not, it raises an error.
    #' - It then checks if \code{modes} is \code{NULL}, and if so, the method
    #'  returns 
    #'   without performing any further checks.
    #' - For each mode in \code{modes}, the method sets \code{self$args$mode} 
    #' to the 
    #'   current mode and calls \code{self$check_mode()}.
    #'
    #' @return NULL. The function is used for its side effects of validating 
    #' each 
    #' mode and raising errors if any validation fails.
    #' 
    check_modes = function() {
      
      modes <- self$args$modes
      designs <- self$args$designs
      
      if (is.null(modes) || is.null(designs)) {
        return(NULL)
      }
      
      if (length(designs) != length(modes)) {
        stop(paste("designs and modes must have the same length."),
             call. = FALSE)
      }
      
      storage <- self$args$mode
      
      for (mode in modes) {
        self$args$mode <- mode
        self$check_mode()
      }
      self$args$mode <- storage
    },
    
    
    #' Check Design Formula
    #'
    #' @description
    #' Validates the design formula ensuring it is a valid character string,
    #' contains allowed characters, includes the intercept term 'X', and 
    #' references
    #' columns present in the metadata.
    #'
    #' @param formula A character string representing the design formula.
    #' @param meta A data frame containing metadata.
    #' @param meta_index An optional index for the data/meta pair.
    #'
    #' @return TRUE if the design formula is valid, otherwise an error is 
    #' thrown.
    #'
    #' @seealso \code{\link[stats]{model.matrix}}
    #'
    check_design_formula = function() {
      
      formula <- self$args[["design"]]
      meta <- self$args$meta
      meta_index <- self$args$meta_index
      
      required_args <- list(formula, meta)
      
      if (any(sapply(required_args, is.null))) {
        return(NULL)
      }

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
        if (!is.null(meta_index)) {
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
    },
    
    
    #' Check Multiple Designs and Metas
    #'
    #' @description
    #' Iterates over multiple design formulas and corresponding metadata
    #' to validate each pair using the `check_design_formula` function.
    #'
    #' @param designs A vector of character strings representing design 
    #' formulas.
    #' @param metas A list of data frames containing metadata.
    #' @param meta_indices A vector of optional indices for the data/meta pairs.
    #'
    #' @return NULL if any check fails, otherwise returns TRUE.
    #'
    check_designs_and_metas = function() {
      
      designs <- self$args$designs
      metas <- self$args$metas
      meta_indices <- self$args$meta_indices
      
      if (is.null(designs) || is.null(metas)) {
        return(NULL)
      }
      
      for (i in seq_along(designs)) {
        self$args$design <- designs[i]
        self$args$meta <- metas[[i]]
        self$args$meta_index <- 
          ifelse(!is.null(meta_indices), meta_indices[i], NA)
        self$check_design_formula()
      }
      
      return(TRUE)
    },
    
    
    #' Check Spline Parameters
    #'
    #' @description
    #' Validates the spline parameters both generally and depending on the
    #' specified mode.
    #'
    #' @param spline_params A list of spline parameters.
    #' @param mode A character string specifying the mode
    #'             ('integrated' or 'isolated').
    #' @param meta A dataframe containing metadata.
    #' @param condition A character string specifying the condition.
    #'
    #' @return No return value, called for side effects.
    #'
    check_spline_params = function() {
      
      spline_params <- self$args$spline_params
      mode <- self$args[["mode"]]
      meta <- self$args$meta
      condition <- self$args$condition
      
      required_args <- list(spline_params, mode, meta, condition)
      
      if (any(sapply(required_args, is.null))) {
        return(NULL)
      }
      
      self$check_spline_params_generally(spline_params)
      self$check_spline_params_mode_dependent(spline_params,
                                              mode,
                                              meta,
                                              condition)
    },
    
    
    #' Check Spline Test Configurations
    #'
    #' @description
    #' This function verifies the spline test configurations and associated 
    #' metadata
    #' within the object's arguments. It performs a series of checks on the 
    #' configurations, including column verification, spline type validation, 
    #' and ensuring that the degrees of freedom (dof) are within acceptable 
    #' ranges.
    #'
    #' @return
    #' Returns `NULL` if any required arguments 
    #' (`spline_test_configs` or `metas`) 
    #' are missing. Otherwise, it performs a series of validation checks.
    #'
    check_spline_test_configs = function() {
      
      spline_test_configs <- self$args$spline_test_configs
      metas <- self$args$metas
      
      required_args <- list(spline_test_configs, metas)
      
      if (any(sapply(required_args, is.null))) {
        return(NULL)
      }
      
      self$check_colums_spline_test_configs(spline_test_configs)
      
      self$check_spline_type_column(spline_test_configs)
      
      self$check_spline_type_params(spline_test_configs)
      
      self$check_max_and_min_dof(spline_test_configs, metas)
    },
    
    
    #' Check Limma Top Tables Structure
    #'
    #' This function checks if the provided limma top tables data structure 
    #' is correctly formatted. It ensures that the data structure contains 
    #' exactly three named elements ('time_effect', 'avrg_diff_conditions', 
    #' and 'interaction_condition_time') and that each element contains 
    #' dataframes with the correct columns and data types.
    #'
    #' @param self An object containing the data structure to check.
    #'
    #' @return This function does not return a value. It stops execution 
    #' if the data structure does not match the expected format.
    #' 
    check_limma_top_tables = function() {
      
      top_tables <- self$args$run_limma_splines_result
      
      required_names <- c("time_effect",
                          "avrg_diff_conditions", 
                          "interaction_condition_time")
      
      if (!all(names(top_tables) %in% required_names) ||
          length(top_tables) != 3) {
        stop("The list must contain exactly three named elements: 
         'time_effect', 'avrg_diff_conditions', 
         'interaction_condition_time'", call. = FALSE)
      }
      
      expected_cols1_3 <- c("X1", "X2", "AveExpr", "F", "P.Value", "adj.P.Val",
                            "feature_index", "feature_names", "intercept")
      
      expected_cols2 <- c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val",
                          "B", "feature_index", "feature_names", "intercept")
      
      for (df in top_tables$time_effect) {
        check_columns(df, expected_cols1_3)
      }
      
      for (df in top_tables$avrg_diff_conditions) {
        check_columns(df, expected_cols2)
      }
      
      for (df in top_tables$interaction_condition_time) {
        check_columns(df, expected_cols1_3)
      }
    },

    
    #' Check Adjusted p-Thresholds
    #'
    #' @description
    #' This function checks the validity of the adjusted p-thresholds vector,
    #' ensuring that
    #' all elements are numeric, greater than 0, and less than 1. If any of 
    #' these
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
    check_adj_pthresholds = function() {
      
      # Exploited argument slicing.
      adj_pthresholds <- self$args$adj_pthresh

      required_args <- list(adj_pthresholds)
      
      if (any(sapply(required_args, is.null))) {
        return(NULL)
      }
      
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
    },
    
    
    #' Check Clusters
    #'
    #' @description
    #' This function verifies the cluster configurations within the object's 
    #' arguments.
    #' It checks if the clusters argument is present and performs validation 
    #' on its 
    #' content. If no clusters are specified, it defaults to automatic cluster 
    #' estimation.
    #'
    #' @details
    #' The function performs the following checks:
    #' - Whether the `clusters` argument is present.
    #' - If `clusters` is a single character string "auto", it defaults to 
    #' automatic 
    #'   cluster estimation and prints a message.
    #' - If `clusters` is not a list or if the list contains non-character and 
    #'   non-numeric types, it throws an error.
    #'
    check_clusters = function() {
      
      clusters <- self$args$clusters
      
      required_args <- list(clusters)
      
      if (any(sapply(required_args, is.null))) {
        return(NULL)
      }
      
      if (is.character(clusters) && 
          length(clusters) == 1 && clusters[1] == "auto") {
        
        print(paste0("No cluster amounts were specified as arguments. Default ",
                     "argument is automatic cluster estimation."))
        
      } else if (!is.list(clusters) ||
                 !all(sapply(clusters, function(x) is.character(x) ||
                             is.numeric(x)))) {
        stop(paste("clusters must be a list containing only character", 
                   "or numeric types."),
             call. = FALSE)
      }
    },
    
    
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
    #' @param time_unit A character string specifying the time unit. Valid 
    #' options
    #' are
    #' 's' for seconds, 'm' for minutes, 'h' for hours, and 'd' for days.
    #'
    #' @return Returns TRUE if the time unit is valid. Stops execution and 
    #' returns
    #' an error message if the time unit is invalid.
    #'
    check_time_unit = function() {
      
      time_unit <- self$args$time_unit
      
      required_args <- list(time_unit)
      
      if (any(sapply(required_args, is.null))) {
        return(NULL)
      }
      
      if (!is.character(time_unit) ||
           length(time_unit) != 1 ||
           nchar(time_unit) > 15) {
        stop(paste0("time_unit must be a single character vector with", 
                    "< 16 letters"),
             call. = FALSE)
      }
    },
    
    
    #' Check and Create Report Directory
    #'
    #' @description
    #' This function checks if the specified report directory exists and is a
    #' valid directory.
    #' If the directory does not exist, it attempts to create it. If there are
    #'  any
    #' warnings or
    #' errors during directory creation, the function stops execution and 
    #' returns
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
    check_report_dir = function() {
      
      report_dir <- self$args$report_dir
      if (is.null(report_dir)) {
        # some functions just have a different name.
        report_dir <- self$args$output_dir  
      }
      
      required_args <- list(report_dir)
      
      if (any(sapply(required_args, is.null))) {
        return(NULL)
      }
      
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
    },
    
    
    #' Check p-Adjustment Method
    #'
    #' @description
    #' This function checks if the provided p-adjustment method is valid. The 
    #' valid
    #' methods are:
    #' "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", and 
    #' "none".
    #' If the method
    #' is not one of these, the function stops execution and returns an error
    #'  message.
    #'
    #' @param padjust_method A character string specifying the p-adjustment 
    #' method.
    #' Valid options
    #' are "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", and
    #' "none".
    #'
    #' @return Returns TRUE if the p-adjustment method is valid. Stops execution
    #'  and
    #' returns an error message if the method is invalid.
    #'
    check_padjust_method = function() {
      
      padjust_method <- self$args$padjust_method
      
      required_args <- list(padjust_method)
      
      if (any(sapply(required_args, is.null))) {
        return(NULL)
      }
      
      supported_methods <- c("holm", "hochberg", "hommel", "bonferroni", "BH",
                             "BY", "fdr", "none")
      if (!(is.character(padjust_method) &&
            padjust_method %in% supported_methods)) {
        stop(sprintf(paste("padjust_method must be a character and one of the",
                           "supported methods (%s).",
                           paste(supported_methods, collapse = ", "))),
             call. = FALSE)
      }
    },
    
    
    #' Check Report Information
    #'
    #' @description
    #' Validates the report information to ensure it contains all mandatory 
    #' fields
    #' and adheres to the required formats.
    #'
    #' @param report_info A named list containing report information.
    #'
    #' @return TRUE if the report information is valid; otherwise, an error is
    #' thrown.
    #'
    check_report_info = function() {
      
      report_info <- self$args$report_info
      
      required_args <- list(report_info)
      
      if (any(sapply(required_args, is.null))) {
        return(NULL)
      }
      
      mandatory_fields <- c("omics_data_type",
                            "data_description",
                            "data_collection_date",
                            "analyst_name",
                            "contact_info",
                            "project_name")
      
      all_fields <- c(mandatory_fields,
                      "method_description",
                      "results_summary",
                      "conclusions")
      
      # Check if report_info is a named list
      if (!is.list(report_info) || is.null(names(report_info))) {
        stop("report_info must be a named list.",
             call. = FALSE)
      }
      
      # Check if all values in report_info are strings
      non_string_fields <- sapply(report_info, function(x) !is.character(x))
      if (any(non_string_fields)) {
        invalid_fields <- names(report_info)[non_string_fields]
        stop(paste("The following fields in report_info must be strings:", 
                   paste(invalid_fields, collapse = ", ")),
             call. = FALSE)
      }
      
      # Check if all mandatory fields are present
      missing_fields <- setdiff(mandatory_fields, names(report_info))
      if (length(missing_fields) > 0) {
        stop(paste("Missing mandatory fields in report_info:",
                   paste(missing_fields, collapse = ", ")),
             call. = FALSE)
      }
      
      # Check if there are any extra fields not in all_fields
      extra_fields <- setdiff(names(report_info), all_fields)
      if (length(extra_fields) > 0) {
        stop(paste("The following fields in report_info are not recognized:",
                   paste(extra_fields, collapse = ", ")),
             call. = FALSE)
      }
      
      # Check omics_data_type format
      if (!grepl("^[a-zA-Z_]+$", report_info[["omics_data_type"]])) {
        stop(paste("The 'omics_data_type' field must contain only alphabetic",
                   "letters and underscores."),
             call. = FALSE)
      }
      
      excluded_fields <- c("data_description",
                           "method_description",
                           "results summary",
                           "conclusions")  
      excluded_limit <- 700  
      
      check_long_fields <- function(data,
                                    excluded_fields,
                                    excluded_limit) {
        
        long_fields <- sapply(data, function(x) {
          if (any(names(data) %in% excluded_fields)) {
            any(nchar(x) > excluded_limit)
          } else {
            any(nchar(x) > 70)
          }
        })
        return(long_fields)
      }
      
      # Check if any field exceeds 70 characters
      long_fields <- check_long_fields(report_info, 
                                       excluded_fields, 
                                       excluded_limit)
      
      if (any(long_fields)) {
        too_long_fields <- names(report_info)[long_fields]
        stop(paste("The following fields have strings exceeding 70 characters:",
                   paste(too_long_fields, collapse = ", "),
                   sep = "\n"), call. = FALSE)
      }
      
      return(TRUE)
    },
    
    
    #' Check Report
    #'
    #' @description
    #' This function verifies the `report` argument within the object's 
    #' arguments.
    #' It checks if the `report` argument is present and validates its 
    #' Boolean value.
    #'
    #' @details
    #' The function performs the following checks:
    #' - Whether the `report` argument is present.
    #' - If `report` is not a Boolean value (`TRUE` or `FALSE`), it throws 
    #' an error.
    #'
    check_report = function() {
      
      report <- self$args$report
      
      required_args <- list(report)
      
      if (any(sapply(required_args, is.null))) {
        return(NULL)
      }
      
      if (report != TRUE && report != FALSE) {
        stop("report must be either Boolean TRUE or FALSE", call. = FALSE)
      }
    }
  )
)



# Level2Functions class --------------------------------------------------------


#' Level2Functions: A class providing level 2 functionalities
#'
#' This class provides various level 2 functionalities, including
#' methods to check dataframes and spline parameters.
#'
#' @seealso \code{\link{InputControl}}
#'
Level2Functions <- R6::R6Class("Level2Functions",
 inherit = Level3Functions,
 
 public = list(
   
   #' Check Data Matrix
   #'
   #' @description
   #' This function checks the validity of the data matrix, ensuring that it 
   #' is a
   #' matrix, contains only numeric values,
   #' has no missing values, and all elements are non-negative. Additionally, it
   #' verifies that no rows or columns are
   #' entirely zeros.
   #'
   #' @param data A dataframe containing numeric values.
   #' @param data_meta_index An optional parameter specifying the index of the
   #'  data
   #' for error messages. Default is NA.
   #'
   #' @return Returns TRUE if all checks pass. Stops execution and returns an 
   #' error
   #' message if any check fails.
   #'
   check_data = function(data, 
                         data_meta_index = NULL) {
     
     # Check if the input is a dataframe
     if (!is.data.frame(data)) {
       stop(self$create_error_message("data must be a dataframe.", 
                                      data_meta_index),
            call. = FALSE)
     }
     
     # Check if all elements in the dataframe are numeric
     if (!all(sapply(data, is.numeric))) {
       stop(self$create_error_message("data must contain only numeric values.",
                                      data_meta_index),
            call. = FALSE)
     }
     
     # Check for missing values
     if (any(is.na(data))) {
       stop(self$create_error_message("data must not contain missing values.",
                                      data_meta_index),
            call. = FALSE)
     }
     
     # Check for non-negative values
     if (any(data < 0)) {
       stop(
         self$create_error_message(
           paste("All elements of data must be ",
                 "non-negative. The elements should",
                 "represent concentrations, abundances, or",
                 "intensities (which are inherently",
                 "non-negative)."),
            data_meta_index),
            call. = FALSE)
     }
     
     # Check for rows with all zeros
     if (any(rowSums(data) == 0)) {
       stop(
         self$create_error_message("data must not contain rows with all zeros.",
                                    data_meta_index),
            call. = FALSE)
     }
     
     # Check for columns with all zeros
     if (any(colSums(data) == 0)) {
       stop(self$create_error_message(paste("data must not contain columns", 
                                            "with all zeros."),
                                      data_meta_index),
            call. = FALSE)
     }
     
     # Check if row headers (rownames) are all strings
     row_headers <- rownames(data)
     if (is.null(row_headers) || !all(sapply(row_headers, is.character))) {
       stop(self$create_error_message("All row headers must be strings.",
                                      data_meta_index),
            call. = FALSE)
     }
     
     return(TRUE)
   },
   
   
   #' Check Metadata
   #'
   #' @description
   #' This function checks the validity of the metadata dataframe, ensuring it
   #' contains the 'Time' column,
   #' does not contain missing values, and that the specified condition column
   #'  is
   #' valid and of the appropriate type.
   #' Additionally, it checks for an optional batch effect column and prints
   #' messages regarding its use.
   #'
   #' @param meta A dataframe containing the metadata, including the 'Time' 
   #' column.
   #' @param condition A single character string specifying the column name 
   #' in the
   #' meta dataframe to be checked.
   #' @param meta_batch_column An optional parameter specifying the column 
   #' name in
   #'  the meta dataframe used to remove the batch effect. Default is NA.
   #' @param meta_batch2_column An optional parameter specifying the column
   #'  name in
   #'  the meta dataframe used to remove the batch effect. Default is NA.
   #' @param data_meta_index An optional parameter specifying the index of the
   #' data/meta pair for error messages. Default is NA.
   #'
   #' @return Returns TRUE if all checks pass. Stops execution and returns an
   #' error message if any check fails.
   #'
   check_meta = function(meta,
                         condition,
                         meta_batch_column = NULL,
                         meta_batch2_column = NULL,
                         data_meta_index = NULL) {
     
     if (!is.data.frame(meta) || 
         !"Time" %in% names(meta) ||
         !is.numeric(meta[["Time"]])) {
       stop(self$create_error_message(paste("meta must be a dataframe with", 
                                            "the numeric column Time"),
                                      data_meta_index),
            call. = FALSE)
     }
     
     self$check_time_column_pattern(meta)
     
     if (any(is.na(meta))) {
       stop(self$create_error_message("meta must not contain missing values.",
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
       stop(self$create_error_message(sprintf("The condition '%s' is not a %s",
                                              condition,
                                              paste("column in meta")),
                                      data_meta_index),
            call. = FALSE)
     }
     
     # Check if the factor column is of appropriate type
     if (!is.factor(meta[[condition]]) &&
         !is.character(meta[[condition]])) {
       stop(self$create_error_message(
       sprintf("The factor column '%s' must be of type
                 factor or character.", condition),
                                      data_meta_index),
            call. = FALSE)
     }
     
     # Check condition and time pattern consistency
     self$check_condition_time_consistency(meta, condition)
     
     if (!is.character(meta_batch2_column)) {
       meta_batch2_column <- NULL
     }
     
     if (is.null(meta_batch_column) && !is.null(meta_batch2_column)) {
       stop(paste("For removing the batch effect, batch2 can only be used when",
                  "batch is used!"), call. = FALSE)
     }
     
     if (!is.null(meta_batch_column)) {
       
       if (meta_batch_column == "Time" || meta_batch_column == condition) {
         stop(paste("meta_batch_column must not be == 'Time' or", condition),
              call. = FALSE)
       }

       self$check_batch_column(meta,
                               meta_batch_column,
                               data_meta_index)
     }
     
     if (!is.null(meta_batch2_column)) {
       
       if (meta_batch2_column == "Time" || meta_batch2_column == condition) {
         stop(paste("meta_batch2_column must not be == 'Time' or", condition),
              call. = FALSE)
       }
       
       if (meta_batch_column == meta_batch2_column) {
         stop(paste("meta_batch_column must not be equal to",
                    "meta_batch2_column"),
              call. = FALSE)
       }
       
       self$check_batch_column(meta,
                               meta_batch2_column,
                               data_meta_index)
     }
     
     return(TRUE)
   },
   
   
   #' Check Dataframe
   #'
   #' @description
   #' Validates that the dataframe contains all required columns with the 
   #' correct data types.
   #'
   #' @param df A dataframe to check.
   #'
   #' @return TRUE if the dataframe is valid, otherwise an error is thrown.
   #'
   check_dataframe = function(df) {
     
     # Define the required columns and their expected types
     required_columns <- list(
       AveExpr = "numeric",
       P.Value = "numeric",
       adj.P.Val = "numeric",
       feature_index = "integer",
       feature_names = "character",
       intercept = "numeric"
     )
     
     # Check if all required columns are present
     missing_columns <- setdiff(names(required_columns), names(df))
     if (length(missing_columns) > 0) {
       stop(paste("Missing columns in top_table:", 
                  paste(missing_columns, collapse = ", ")),
            call. = FALSE)
     }
     
     # Check if columns have the correct type
     for (col in names(required_columns)) {
       if (!inherits(df[[col]], required_columns[[col]])) {
         stop(paste("top_table column", col, "must be of type",
                    required_columns[[col]]),
              call. = FALSE)
       }
     }
     
     return(TRUE)
   },
   
   
   #' Check Spline Parameters Generally
   #'
   #' @description
   #' Validates the general structure and contents of spline parameters.
   #'
   #' @param spline_params A list of spline parameters.
   #'
   #' @return No return value, called for side effects.
   #'
   check_spline_params_generally = function(spline_params) {
     if ("spline_type" %in% names(spline_params)) {
       if (!all(spline_params$spline_type %in% c("b", "n"))) {
         stop(
           paste("Elements of spline_type must be either 'b' for B-splines or",
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
       stop(
         "Either dof or knots must be present, but not both,in spline_params.",
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
           stop(
             paste("B-splines require DoF > 2, for spline_params spline_type",
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
   },
   
   
   #' Check Spline Parameters Mode Dependent
   #'
   #' @description
   #' Validates the spline parameters depending on the specified mode.
   #'
   #' @param spline_params A list of spline parameters.
   #' @param mode A character string specifying the mode
   #'            ('integrated' or 'isolated').
   #' @param meta A dataframe containing metadata.
   #' @param condition A character string specifying the condition.
   #'
   #' @return No return value, called for side effects.
   #'
   #'
   check_spline_params_mode_dependent = function(spline_params,
                                                 mode,
                                                 meta,
                                                 condition) {
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
           stop(
             paste("All elements in 'knots' in spline_params must have length",
                   "1 when mode is 'integrated'. Different spline parameters",
                   "for the different levels is not supported for this mode"),
                call. = FALSE)
         }
       }
       if ("bknots" %in% names(spline_params)) {
         if (length(spline_params$bknots) != 1) {
           stop(
             paste("All elements in 'bknots' in spline_params must have length",
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
           stop(
             paste("'knots' in spline_params must have the same number of",
                   "elements as there are unique elements in the ", condition,
                   "column of meta when mode is 'isolated'."),
                call. = FALSE)
         }
       }
       if ("bknots" %in% names(spline_params)) {
         if (length(spline_params$bknots) != num_levels) {
           stop(
             paste("'bknots' in spline_params must have the same number of",
                   "elements as there are unique elements in the ", condition,
                   "column of meta when mode is 'isolated'."),
                call. = FALSE)
         }
       }
     }
   },
   
   
   #' Check Columns in Spline Test Configurations
   #'
   #' @description
   #' Validates that the spline test configurations contain the required columns 
   #' in the correct order.
   #'
   #' @param spline_test_configs A dataframe containing spline test
   #'  configurations.
   #'
   #' @return No return value, called for side effects.
   #' 
   check_colums_spline_test_configs = function(spline_test_configs) {
     
     required_columns <- c("spline_type", "degree", "dof", "knots", "bknots")
     
     # Check for exact match of column names and order
     if (!identical(names(spline_test_configs), required_columns)) {
       # Find the missing or extra columns
       missing_columns <- setdiff(required_columns, names(spline_test_configs))
       extra_columns <- setdiff(names(spline_test_configs), required_columns)
       error_message <- "Error: Incorrect columns in dataframe. "
       
       # Append specific issues to the error message
       if (length(missing_columns) > 0) {
         error_message <- paste0(error_message, "Missing columns: ",
                                 paste(missing_columns, collapse=", "), ". ")
       }
       if (length(extra_columns) > 0) {
         error_message <- paste0(error_message, "Extra columns: ",
                                 paste(extra_columns, collapse=", "), ". ")
       }
       error_message <- paste0(error_message,
                               "Expected columns in order: ",
                               paste(required_columns, collapse=", "), ".")
       
       stop(error_message)
     }
   },
   
   
   #' Check Spline Type Column
   #'
   #' @description
   #' Validates that the 'spline_type' column in the spline test configurations 
   #' contains only 'n' or 'b'.
   #'
   #' @param spline_test_configs A dataframe containing spline test 
   #' configurations.
   #'
   #' @return No return value, called for side effects.
   #' 
   check_spline_type_column = function(spline_test_configs) {
     
     if (!all(spline_test_configs$spline_type %in% c("n", "b"))) {
       # Identify invalid entries
       invalid_entries <- spline_test_configs$spline_type[
         !spline_test_configs$spline_type %in% c("n", "b")
       ]
       error_message <- sprintf(
         "Error: 'spline_type' contains invalid entries. 
        Only 'n' or 'b' are allowed. Invalid entries found: %s",
         paste(unique(invalid_entries), collapse=", ")
       )
       
       stop(error_message)
     }
   },
   
   
   #' Check Spline Type Parameters
   #'
   #' @description
   #' Validates the parameters for each row in the spline test configurations 
   #' based on the spline type.
   #'
   #' @param spline_test_configs A dataframe containing spline test 
   #' configurations.
   #'
   #' @return TRUE if all checks pass, otherwise an error is thrown.
   #' 
   check_spline_type_params = function(spline_test_configs) {
     
     for (i in seq_len(nrow(spline_test_configs))) {
       row <- spline_test_configs[i,]
       switch(as.character(row$spline_type),
              "n" = {
                if (!is.na(row$degree)) 
                  stop("degree must be NA for spline_type n")
                if (!((is.na(row$dof) && !is.na(row$knots)) ||
                      (!is.na(row$dof) && is.na(row$knots)))) {
                  stop("Either dof or knots must be NA, but not both, for 
                    spline_type n")
                }
                if (!is.na(row$dof) && !is.integer(row$dof)) {
                  stop("dof must be an integer when it is not NA for 
                    spline_type n")
                }
                if (!is.na(row$knots) && !is.numeric(row$knots)) {
                  stop("knots must be a numeric vector when it is not NA for 
                    spline_type n")
                }
                if (!is.na(row$bknots) && (!is.numeric(row$bknots) || 
                                           length(row$bknots) != 2)) {
                  stop("bknots must be a numeric vector of exactly two elements 
                    or NA for spline_type n")
                }
              },
              "b" = {
                if (!is.integer(row$degree)) stop("degree must be an integer for 
                                               spline_type b")
                if (!is.na(row$dof) && (!is.integer(row$dof) || 
                                        row$dof < row$degree)) {
                  stop("dof must be an integer at least as big as degree for 
                    spline_type b")
                }
                if (!is.na(row$knots) && !is.numeric(row$knots)) {
                  stop("knots must be a numeric vector when it is not NA for 
                    spline_type b")
                }
                if (!is.na(row$bknots) && (!is.numeric(row$bknots) || 
                                           length(row$bknots) != 2)) {
                  stop("bknots must be a numeric vector of exactly two elements 
                    or NA for spline_type b")
                }
              },
              stop("spline_type must be either 'n' or 'b'")
       )
     }
     return(TRUE)
   },
   
   
   #' Check Maximum and Minimum Degrees of Freedom
   #'
   #' @description
   #' Validates the degrees of freedom (DoF) for each row in the spline test 
   #' configurations based on the metadata.
   #'
   #' @param spline_test_configs A dataframe containing spline test 
   #' configurations.
   #' @param metas A list of metadata corresponding to the data matrices.
   #'
   #' @return No return value, called for side effects.
   #' 
   check_max_and_min_dof = function(spline_test_configs, 
                                     metas) {
     
     for (i in seq_len(nrow(spline_test_configs))) {
       row <- spline_test_configs[i, ]
       spline_type <- row[["spline_type"]]
       dof <- row[["dof"]]
       degree <- row[["degree"]]
       knots <- row[["knots"]]
       
       # Calculate k and DoF if dof is NA
       if (is.na(dof)) {
         k <- length(knots)
         if (spline_type == "b") {
           dof <- k + degree
         } else if (spline_type == "n") {
           dof <- k + 1
         } else {
           stop("Unknown spline type '", spline_type, "' at row ", i)
         }
       }
       
       # Check if calculated or provided DoF is valid
       if (dof < 2) {
         stop("DoF must be at least 2, found DoF of ", dof, " at row ", i)
       }
       
       for (j in seq_along(metas)) {
         meta <- metas[[j]]
         nr_timepoints <- length(unique(meta$Time))
         if (dof > nr_timepoints) {
           stop("DoF (", dof, ") cannot exceed the number of unique time points 
           (", nr_timepoints, ") in meta", j, " at row ", i)
         }
       }
       
     }
     
     invisible(NULL)
  },
  
  
  #' Check Dataframe Columns
  #'
  #' This function checks if the columns of a dataframe match the expected 
  #' column names and their respective data types.
  #'
  #' @param df A dataframe to check.
  #' @param expected_cols A character vector of expected column names.
  #'
  #' @return This function does not return a value. It stops execution if the
  #' dataframe columns or their classes do not match the expected structure.
  #' 
  check_columns = function(df, 
                           expected_cols) {
    
    actual_cols <- names(df)
    if (!all(expected_cols %in% actual_cols)) {
      stop("Dataframe columns do not match expected structure",
           call. = FALSE)
    }
    expected_classes <- c("numeric", "numeric", "numeric", "numeric",
                          "numeric", "numeric", "integer", "character",
                          "numeric")
    actual_classes <- sapply(df, class)
    if (!all(actual_classes == expected_classes)) {
      stop("Dataframe column classes do not match expected classes",
           call. = FALSE)
    }
  }
 )
)



# Level3Functions class --------------------------------------------------------


#' Level3Functions: A class for level 3 utility functions
#'
#' This class provides methods for creating error messages and checking
#' batch columns.
#' 
#' @seealso \code{\link{Level2Functions}}
#'
Level3Functions <- R6::R6Class("Level3Functions",

 inherit = Level4Functions,

 public = list(


   #' Check Batch Column
   #'
   #' @description
   #' This method checks the batch column in the metadata and provides 
   #' appropriate messages.
   #'
   #' @param meta A dataframe containing metadata.
   #' @param meta_batch_column A character string specifying the batch column
   #'  in the metadata.
   #' @param data_meta_index An optional parameter specifying the index of the
   #'  data/meta pair. Default is NA.
   #'
   #' @return NULL. The method is used for its side effects of throwing errors 
   #' or printing messages.
   #'
   check_batch_column = function(meta,
                                 meta_batch_column,
                                 data_meta_index) {

     if (!is.null(meta_batch_column) && !(meta_batch_column %in% names(meta))) {
       stop(self$create_error_message(sprintf("Batch effect column '%s' %s",
                                              meta_batch_column,
                                              "not found in meta"),
                                      data_meta_index),
            call. = FALSE)
     } else if (!is.null(meta_batch_column)) {
       if (!is.null(data_meta_index)) {
         message(sprintf("Index: %s. %s",
                         data_meta_index,
                         paste("Column", meta_batch_column,
                               "of meta will be used",
                               "to remove the batch effect for the plotting")))
       } else {
         message(sprintf("Column '%s' of meta will be used to %s",
                         meta_batch_column,
                         paste("remove the batch effect for the plotting")))
       }
     } else {
       if (!is.null(data_meta_index)) {
         message(sprintf(
           "Index: %s. Batch effect will NOT be removed for plotting!",
           data_meta_index))
       } else {
         message("Batch effect will NOT be removed for plotting!")
       }
     }
   },
   
   
   #' Check Time Column Pattern
   #'
   #' @description
   #' This function checks if the time values in the metadata follow a 
   #' consistent 
   #' pattern with a fixed number of replicates. It ensures that the time values 
   #' change consistently after each block of replicates and identifies the 
   #' repeating pattern.
   #'
   #' @param meta A data frame containing a column named "Time" with time 
   #'             values.
   #'
   #' @return Returns TRUE if the time values follow a consistent pattern. It 
   #'         also outputs messages if inconsistencies are found.
   #'
   check_time_column_pattern = function(meta) {
     time_values <- meta[["Time"]]
     
     # Find the number of replicates
     replicate_count <- 1
     for (i in 2:length(time_values)) {
       if (time_values[i] != time_values[1]) {
         replicate_count <- i - 1
         break
       }
     }

     # Check that the value changes after the replicate_count consistently
     for (i in seq(replicate_count + 1, length(time_values), 
                   by = replicate_count)) {
       if (i + replicate_count - 1 > length(time_values)) break
       current_block <- time_values[i:(i + replicate_count - 1)]
       previous_value <- time_values[i - 1]
       if (!all(current_block == current_block[1]) ||
           current_block[1] <= previous_value) {
         message(paste("The time values do not change consistently", 
                       "after the replicates."))
         break
       }
     }
     
     # Find the pattern
     pattern_end <- 1
     for (i in (replicate_count + 1):length(time_values)) {
       if (time_values[i] < time_values[i - 1]) {
         pattern_end <- i - 1
         break
       }
     }
     
     pattern <- time_values[1:pattern_end]
     pattern_length <- length(pattern)
     
     # Check if the pattern is fully repeated
     for (i in seq(pattern_length + 1, length(time_values),
                   by = pattern_length)) {
       if (i + pattern_length - 1 > length(time_values)) break
       if (!all(time_values[i:(i + pattern_length - 1)] == pattern)) {
         message("The time pattern is not fully repeated.")
         break
       }
     }
     
     # Check if there are leftover values
     if (length(time_values) %% pattern_length != 0) {
       message(paste("The time pattern is not fully repeated.", 
                     "There are leftover values."))
     }
     
     return(TRUE)
   },
   
   
   #' Check Condition Time Consistency
   #'
   #' @description
   #' This function checks whether the values in the `condition` column
   #' have unique values for each block of identical `Time` values in the 
   #' `meta` dataframe.
   #' Additionally, it ensures that every new block of a given time has a 
   #' new value
   #' in the `condition` column.
   #'
   #' @param meta A dataframe containing the metadata, including the `Time`
   #'  column.
   #' @param condition A character string specifying the column name in `meta`
   #'                  used to define groups for analysis.
   #'
   #' @return Logical TRUE if the condition values are consistent with the
   #'  time series pattern.
   #'
   check_condition_time_consistency = function(meta,
                                               condition) {
     
     # Get the unique times in the order they appear
     unique_times <- unique(meta[["Time"]])
     
     # Initialize a list to store previously seen conditions for each time
     # segment
     seen_conditions <- list()
     
     # Iterate through each block of unique times
     for (time in unique_times) {
       # Get the indices for the current time segment
       current_indices <- which(meta[["Time"]] == time)
       
       # Get the condition values for the current time segment
       current_conditions <- meta[current_indices, condition, drop = TRUE]
       
       # Ensure that all conditions in the current segment are the same
       if (length(unique(current_conditions)) != 2) {
         stop(paste("Every block of identical time values in meta must",
                    "have unique values in the column", condition),
              call. = FALSE)
       }
       
       # Check if the condition value has been seen before for this time segment
       if (!is.null(seen_conditions[[as.character(time)]])) {
         if (current_conditions[1] %in% seen_conditions[[as.character(time)]]) {
           stop(sprintf("Condition '%s' for time '%s' has been seen before.",
                        current_conditions[1], time), call. = FALSE)
         }
       }
       
       # Update the seen conditions for the current time segment
       if (is.null(seen_conditions[[as.character(time)]])) {
         seen_conditions[[as.character(time)]] <- character(0)
       }
       seen_conditions[[as.character(time)]] <- c(
         seen_conditions[[as.character(time)]],
         current_conditions[1]
       )
     }
     
     return(TRUE)
   }
 )
)



# Level4Functions class --------------------------------------------------------

#' Level4Functions: A class for level 3 utility functions
#'
#' This class provides methods for creating error messages and checking
#' batch columns.
#'
#' @seealso \code{\link{Level3Functions}}
#'
Level4Functions <- R6::R6Class("Level4Functions",
 public = list(

   #' Create Error Message
   #'
   #' @description
   #' This method creates a formatted error message that includes the index of
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
   create_error_message = function(message,
                                   data_meta_index = NULL) {

     if (!is.null(data_meta_index)) {
       return(sprintf("data/meta pair index %d: %s", data_meta_index, message))
     } else {
       return(message)
     }
   }
 )
)



# Utility input control functions ----------------------------------------------


#' Check for NULL Elements in Arguments
#'
#' This function checks if any elements in the provided list of arguments 
#' are `NULL`.
#' If any `NULL` elements are found, it stops the execution and returns 
#' an informative error message.
#'
#' @param args A list of arguments to check for `NULL` elements.
#'
#' @return This function does not return a value. It stops execution if 
#' any `NULL` elements are found in the input list.
#'
check_null_elements <- function(args) {
  
  if (any(sapply(args, is.null))) {
    stop("One or more function arguments are NULL.",
         call. = FALSE)
  }
}
