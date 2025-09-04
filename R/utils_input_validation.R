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
      if (!is.null(args$splineomics) &&
        inherits(args$splineomics, "SplineOmics")) {
        args <- c(args, args$splineomics)
        args$splineomics <- NULL
      }

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
    #' - \code{self$check_design_formula()}
    #' - \code{self$check_mode()}
    #' - \code{self$check_modes()}
    #' - \code{self$check_designs_and_metas()}
    #' - \code{self$check_spline_params()}
    #' - \code{self$check_spline_test_configs()}
    #' - \code{self$check_adj_pthresholds()}
    #' - \code{self$check_clusters()}
    #' - \code{self$check_time_unit()}
    #' - \code{self$check_raw_data()}
    #' - \code{self$check_report_dir()}
    #' - \code{self$check_padjust_method()}
    #' - \code{self$check_report_info()}
    #' - \code{self$check_report()}
    #' - \code{self$check_feature_name_columns()}
    #'
    #' @return NULL. The function is used for its side effects of validating
    #' input
    #' arguments and raising errors if any validation fails.
    #'
    auto_validate = function() {
      self$check_data_and_meta()
      self$check_annotation()
      self$check_alphas()
      self$check_batch_effects()
      self$check_datas_and_metas()
      self$check_datas_descr()
      self$check_top_tables()
      self$check_design_formula()
      self$check_robust_fit()
      self$check_dream_params()
      self$check_mode()
      self$check_modes()
      self$check_designs_and_metas()
      self$check_spline_params()
      self$check_spline_test_configs()
      self$check_adj_pthresholds()
      self$check_adj_pthresh_limma_category_2_3()
      self$check_clusters()
      self$check_genes()
      self$check_plot_info()
      self$check_plot_options()
      self$check_raw_data()
      self$check_report_dir()
      self$check_padjust_method()
      self$check_report_info()
      self$check_report()
      self$check_support()
      self$check_feature_name_columns()
    },


    #' Check Data and Meta
    #' 
    #' @noRd
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

      if (any(vapply(required_args, is.null, logical(1)))) {
        return(NULL)
      }

      self$check_data(
        data,
        data_meta_index
      )

      self$check_meta(
        meta = meta,
        condition = condition,
        meta_batch_column = meta_batch_column,
        meta_batch2_column = meta_batch2_column,
        data_meta_index = data_meta_index
      )

      if (!(nrow(meta) == ncol(data))) {
        if (!is.null(data_meta_index)) {
          stop(
            paste0(
              "For index ", data_meta_index,
              "data column number must be equal to ",
              "meta row number"
            ),
            call. = FALSE
          )
        } else {
          stop(paste0("data column number must be equal to meta row number"),
            call. = FALSE
          )
        }
      }
    },


    #' Check Annotation Consistency
    #' 
    #' @noRd
    #'
    #' @description
    #' This method checks the consistency of the annotation with the data.
    #' It ensures
    #' that the annotation is a dataframe and that it has the same number
    #' of rows as the data.
    #'
    #' @details
    #' The method performs the following checks:
    #'
    #' * Ensures that both `annotation` and `data` are provided.
    #' * Confirms that `annotation` is a dataframe.
    #' * Verifies that `annotation` and `data` have the same number of rows.
    #'
    #' If any of these checks fail, an informative error message is returned.
    #'
    #' @return NULL if any required arguments are missing. Otherwise, performs
    #' checks and potentially raises errors if checks fail.
    #'
    check_annotation = function() {
      annotation <- self$args[["annotation"]]
      data <- self$args[["data"]]

      required_args <- list(annotation, data)

      if (any(vapply(required_args, is.null, logical(1)))) {
        return(NULL)
      }

      if (!is.data.frame(annotation)) {
        stop(
          "annotation is not a dataframe but must be one!",
          call. = FALSE
        )
      }

      if (nrow(annotation) != nrow(data)) {
        stop(
          "annotation and data don't have the same nr. of rows but must have!",
          call. = FALSE
        )
      }
    },
    
    
    #' Check Alphas Argument
    #'
    #' @noRd
    #'
    #' @description
    #' This method validates the `alphas` argument, ensuring it is either a 
    #' single numeric threshold or a named list of numeric thresholds to be used
    #' for determining significance levels.
    #'
    #' @details
    #' The method performs the following checks:
    #'
    #' * Accepts either a numeric scalar or a named list of numeric values.
    #' * If a list is provided, it must be named with non-empty names.
    #' * All values must be between 0 and 1 and not `NA`.
    #'
    #' If any of these checks fail, an informative error is raised.
    #'
    #' @return NULL if `alphas` is not provided. Otherwise, performs
    #' validation checks and raises errors if checks fail.
    #' 
    check_alphas = function() {
      alphas <- self$args[["alphas"]]
      
      if (is.null(alphas)) {
        return(NULL)
      }
      
      if (is.numeric(alphas)) {
        if (length(alphas) != 1 || is.na(alphas) || alphas < 0 || alphas > 1) {
          stop_call_false(
            "alphas must be a single numeric value between 0 and 1."
            )
        }
      } else if (is.list(alphas)) {
        if (length(alphas) == 0 ||
            any(!vapply(alphas, is.numeric, logical(1)))) {
          stop_call_false(
            "alphas must be a named list of numeric values."
            )
        }
        
        if (any(is.na(unlist(alphas))) ||
            any(unlist(alphas) < 0) ||
            any(unlist(alphas) > 1)) {
          stop_call_false(
            "All values in the alphas list must be numeric & between 0 and 1."
            )
        }
        
        if (is.null(names(alphas)) || any(names(alphas) == "")) {
          stop_call_false(
            "alphas list must be named and contain no empty names."
            )
        }
      } else {
        stop_call_false(paste0(
          "alphas must be either a numeric scalar or a named list of numeric", 
          "values."
        ))
      }
    },
    
    
    #' Check Batch Effects Argument
    #'
    #' @noRd
    #'
    #' @description
    #' This method validates the `batch_effects` argument of the `find_pvc()`
    #' function. It ensures that the argument is either `FALSE` or a character
    #' vector containing up to two valid column names for batch correction.
    #'
    #' @details
    #' The method performs the following checks:
    #'
    #' * Accepts either `FALSE` (no batch correction), or a character vector.
    #' * If a character vector is provided, it must have length 1 or 2.
    #' * All entries must be non-empty strings.
    #'
    #' This validation does not check whether the provided names exist in
    #' the metadata, as that depends on runtime context.
    #'
    #' @return NULL if `batch_effects` is `FALSE`. Otherwise, performs
    #' validation checks and raises an error if checks fail.
    #'
    check_batch_effects = function() {
      batch_effects <- self$args[["batch_effects"]]

      if (identical(batch_effects, FALSE) || is.null(batch_effects)) {
        return(NULL)
      }
      
      if (!is.character(batch_effects)) {
        stop_call_false(
          "batch_effects must be either FALSE or a character vector."
          )
      }
      
      if (length(batch_effects) < 1 || length(batch_effects) > 2) {
        stop_call_false("batch_effects must have one or two elements.")
      }
      
      if (any(batch_effects == "") || any(is.na(batch_effects))) {
        stop_call_false("batch_effects must not contain empty or NA entries.")
      }
    },
    

    #' Check Multiple Data and Meta Pairs
    #' 
    #' @noRd
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

      if (any(vapply(required_args, is.null, logical(1)))) {
        return(NULL)
      }

      if (length(datas) != length(metas)) {
        stop(paste0("datas and metas must have the same length."),
          call. = FALSE
        )
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
    #' @noRd
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

      if (any(vapply(required_args, is.null, logical(1)))) {
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
    
    
    #' Check `bp_cfg`
    #'
    #' @noRd
    #'
    #' @description
    #' Validates that `bp_cfg` is either **`NULL`** *or* a length-2 numeric
    #' (integer) vector **named exactly**  
    #' `c("n_cores", "blas_threads")`, where each value is a positive whole
    #' number (>= 1). Any deviation triggers an informative `stop()`.
    #'
    #' @param bp_cfg A numeric named vector with elements
    #'   `"n_cores"` and `"blas_threads"`, or `NULL`.
    #'
    #' @return No return value. Called for its side-effect of throwing an error
    #'  on malformed input.
    #'
    #' @seealso
    #'   \code{\link[base]{stop}} for error signalling.
    #'
    check_bp_cfg = function() {
      
      bp_cfg <- self$args[["bp_cfg"]]
      
      required_args <- list(bp_cfg)
      
      if (any(vapply(required_args, is.null, logical(1)))) {
        return(NULL)
      }
      
      # Allow NULL (falls back to defaults elsewhere)
      if (is.null(bp_cfg))
        return(invisible(NULL))
      
      ## Basic structure checks
      valid_structure <-
        is.numeric(bp_cfg) &&
        length(bp_cfg) == 2 &&
        setequal(names(bp_cfg), c("n_cores", "blas_threads"))
      
      if (!valid_structure) {
        stop(
          "'bp_cfg' must be NULL or a numeric vector with *exactly* two named ",
          "elements: c(n_cores = , blas_threads = )."
        )
      }
      
      # Value checks
      if (any(is.na(bp_cfg)) ||
          any(bp_cfg < 1)    ||
          any(bp_cfg != as.integer(bp_cfg))) {
        
        offending <- names(bp_cfg)[
          which(bp_cfg < 1 | bp_cfg != as.integer(bp_cfg))]
        stop(
          "'bp_cfg' values must be positive whole numbers (>=1). ",
          "Offending element(s): ", paste(offending, collapse = ", ")
        )
      }
      
      invisible(NULL)
    },
    

    #' Check Top Tables
    #' 
    #' @noRd
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

      if (any(vapply(required_args, is.null, logical(1)))) {
        return(NULL)
      }

      # Helper function to check data frames in a list
      check_list_of_dataframes <- function(df_list) {
        if (!is.list(df_list) ||
          !all(vapply(df_list, is.data.frame, logical(1)))) {
          stop_call_false("Expected a list of dataframes")
        } else {
          if (length(df_list) != 2) {
            stop_call_false(paste(
              "top_tables must be a list of two lists if you want",
              "to cluster the hits of the limma results of",
              "avrg_diff_conditions or interaction_condition_time.",
              "The list must contain one of those two plus",
              "the time_effect limma result, in any order"
            ))
          }

          underscore_count <- vapply(
            names(df_list),
            function(name) sum(grepl("_", name)),
            integer(1)
          )
          if (sum(underscore_count == 1) != 1 ||
            sum(underscore_count == 4) != 1) {
            stop_call_false(paste(
              "top_tables must be a list of two lists if you want",
              "to cluster the hits of the limma results of",
              "avrg_diff_conditions or interaction_condition_time.",
              "The list must contain one of those two plus",
              "the time_effect limma result, in any order"
            ))
          }

          for (df in df_list) {
            self$check_dataframe(df)
          }
        }
      }

      # Check if top_tables is a list
      if (!is.list(top_tables)) {
        stop_call_false("top_tables must be a list")
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
          matches <- gregexpr("_", element_name)
          underscore_count <- sum(vapply(
            matches,
            function(x) if (x[1] == -1) 0 else length(x),
            integer(1)
          ))

          if (underscore_count != 1) {
            stop(
              paste(
                "Very likely you did not pass the time_effect result",
                "of limma but rather one of avrg_diff_conditions or",
                "interaction_condition_time, which cannot be passed alone",
                "(to cluster the hits of one of those, put only one of them",
                "at a time into a list with the time_effect result. Further",
                "note: Please do not edit those list element names by hand."
              ),
              call. = FALSE
            )
          }

          self$check_dataframe(element)
        } else {
          stop(
            paste(
              "top_tables must contain either data frames or lists of",
              "data frames"
            ),
            call. = FALSE
          )
        }
      }
    },


    #' Check Design Formula
    #' 
    #' @noRd
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
      meta <- self$args[["meta"]]
      meta_index <- self$args[["meta_index"]]

      # Not strictly required
      meta_batch_column <- self$args[["meta_batch_column"]]
      meta_batch2_column <- self$args[["meta_batch2_column"]]

      required_args <- list(
        formula,
        meta
      )

      if (any(vapply(required_args, is.null, logical(1)))) {
        return(NULL)
      }
      
      effects <- extract_effects(formula)

      # Check if the formula is a valid character string
      if (!is.character(formula) || length(formula) != 1) {
        stop("The design formula must be a valid character string.",
          call. = FALSE
        )
      }

      # Ensure the formula contains allowed characters only
      allowed_chars <- "^[~ 1A-Za-z0-9_+*:()-|]*$"
      if (!grepl(allowed_chars, formula)) {
        stop("The design formula contains invalid characters.",
          call. = FALSE
        )
      }

      # # Ensure that the formula begins with an intercept (~ 1)
      # # Ignore whitespace, check the start of the string
      # if (!grepl("^\\s*~\\s*1", formula)) {
      #   stop(
      #     paste(
      #       "The design formula must start with an intercept term '~ 1'.",
      #       "This is because spline curves are plotted onto the data",
      #       "which is not possible without an intercept"
      #     ),
      #     call. = FALSE
      #   )
      # }

      # Ensure the formula contains the intercept term 'X'
      if (!grepl("\\bTime\\b", formula)) {
        stop_call_false(
          "The design formula must include the term 'Time' as a fixed effect"
        )
      }
      
      formatted_formula <- gsub("1\\|", "", formula) 
      # Extract terms from the formula (removing interactions and functions)
      formula_terms <- unlist(strsplit(gsub(
        "[~+*:()]",
        " ",
        formatted_formula
        ), " "))
      formula_terms <- formula_terms[formula_terms != ""]

      # Remove '1' and 'X' from terms since they are not columns
      formula_terms <- setdiff(formula_terms, c("1"))

      # Check if the terms are present in the dataframe
      missing_columns <- setdiff(formula_terms, names(meta))
      if (length(missing_columns) > 0) {
        if (!is.null(meta_index)) {
          stop_call_false(sprintf(
            "%s (data/meta pair index: %s): %s",
            "The following design columns are missing in meta",
            meta_index,
            paste(missing_columns, collapse = ", ")
          ))
        } else {
          stop_call_false(paste(
            "The following design columns are missing in meta:",
            paste(missing_columns, collapse = ", ")
          ))
        }
      }

      # Convert formula to string for regex checking
      formula_str <- as.character(formula)

      # Check if batch column is provided and validate its presence in the formula
      if (!is.null(meta_batch_column)) {
        if (!grepl(meta_batch_column, formula_str)) {
          stop_call_false(
            paste(
              "The batch effect column", meta_batch_column,
              "is provided in the SplineOmics object as the field",
              "meta_batch_column but not present in the design formula. ",
              "Please ensure that if you specify a batch column, ",
              "it is included in the design formula to",
              "remove batch effects."
            )
          )
        }
      }

      # Check if the second batch column is provided and validate its presence
      if (!is.null(meta_batch2_column)) {
        if (!grepl(meta_batch2_column, formula_str)) {
          stop_call_false(
            paste(
              "The second batch effect column", meta_batch2_column,
              "is provided but not present in the design formula.",
              "Please ensure that if you specify a second batch column,",
              "it is included in the design formula",
              "to remove batch effects."
            )
          )
        }
      }

      return(TRUE)
    },
    
    
    #' Check validity of the `robust_fit` argument
    #' 
    #' @noRd
    #'
    #' @description
    #' This internal method validates the `robust_fit` argument provided to the
    #' function. It ensures that the value is either `TRUE`, `FALSE`, or `NULL`,
    #' which control whether robust modeling should be applied in downstream 
    #' analysis.
    #'
    #' If the value is not one of these types, the function halts execution 
    #' using `stop_call_false()` with an informative error message.
    #'
    #' @return No return value. Called for its side effect of input validation.
    #'
    check_robust_fit = function() {
      robust_fit <- self$args[["robust_fit"]]
      
      if (!is.null(robust_fit) && !is.logical(robust_fit)) {
        stop_call_false(
          "'robust_fit' must be either TRUE, FALSE, or NULL."
        )
      }
    },
    
      
    #' Validate the `dream_params` argument
    #' 
    #' @noRd
    #'
    #' @description
    #' This function checks the validity of the `dream_params` argument provided 
    #' in the class. If `dream_params` is present, it ensures that it contains  
    #' only the allowed optional elements `dof` and `KenwardRoger` in the correct 
    #' format. Unnamed elements or elements other than these two are not allowed.
    #'
    #' @return
    #' Returns `TRUE` if `dream_params` passes all checks. Otherwise, stops the 
    #' function and returns an error message using `stop_call_false`.
    #'
    check_dream_params = function() {
      dream_params <- self$args[["dream_params"]]
      
      required_args <- list(dream_params)
      
      if (any(vapply(required_args, is.null, logical(1)))) {
        return(NULL)
      }
      
      # Check that dream_params is a named list
      if (!is.list(dream_params) || is.null(names(dream_params))) {
        stop_call_false("dream_params must be a named list.")
      }
      
      # Define allowed elements and check for unexpected elements
      allowed_elements <- c("dof", "KenwardRoger")
      if (!all(names(dream_params) %in% allowed_elements)) {
        stop_call_false(
          "dream_params contains invalid elements. Only 'dof' and 'KenwardRoger'
      are allowed."
        )
      }
      
      # If 'dof' is provided, check that it is an integer greater than 1
      if ("dof" %in% names(dream_params)) {
        if (!is.numeric(dream_params[["dof"]]) 
            || dream_params[["dof"]] <= 1 
            || dream_params[["dof"]] != as.integer(dream_params[["dof"]])) {
          stop_call_false("'dof' must be an integer greater than 1.")
        }
      }
      
      # If 'KenwardRoger' is provided, check that it is a boolean
      if ("KenwardRoger" %in% names(dream_params)) {
        if (!is.logical(dream_params[["KenwardRoger"]]) 
            || length(dream_params[["KenwardRoger"]]) != 1) {
          stop_call_false("'KenwardRoger' must be a boolean.")
        }
      }
      
      # No unnamed elements should be present
      if (any(names(dream_params) == "")) {
        stop_call_false("Unnamed elements are not allowed in dream_params.")
      }
      
      return(TRUE)  
    },


    #' Validate and check all modes
    #' 
    #' @noRd
    #'
    #' @description
    #' This function iterates over the `modes` argument, sets each `mode` in
    #' `self$args`, and calls `check_mode()` to validate each mode. After each
    #' validation, the `mode` is removed from `self$args`.
    #'
    #' @return NULL if `modes` is missing; otherwise, checks all modes.
    #'
    check_modes = function() {
      modes <- self$args[["modes"]]

      required_args <- list(modes)

      if (any(vapply(required_args, is.null, logical(1)))) {
        return(NULL)
      }

      for (mode in modes) {
        self$args$mode <- mode

        self$check_mode()

        self$args$mode <- NULL
      }
    },


    #' Check the mode argument for validity
    #' 
    #' @noRd
    #'
    #' @description
    #' This function checks if the `mode` argument is provided and validates
    #' that it is either "isolated" or "integrated". If `mode` is missing or
    #' invalid, an error is thrown.
    #'
    #' @return NULL if `mode` is missing; otherwise, validates the mode.
    #'
    check_mode = function() {
      mode <- self$args[["mode"]]
      meta <- self$args[["meta"]]
      condition <- self$args[["condition"]]
      design <- self$args[["design"]]

      required_args <- list(
        mode,
        meta,
        condition,
        design
        )

      if (any(vapply(required_args, is.null, logical(1)))) {
        return(NULL)
      }

      if (mode != "isolated" && mode != "integrated") {
        stop_call_false(
          "mode must be either 'isolated' or 'integrated' and not '", mode, "'!"
        )
      }
      
      # Convert design to character, if it's a formula
      design_str <- if (inherits(design, "formula")) {
        deparse(design)
      } else {
        as.character(design)
      }
      design_str <- paste(design_str, collapse = " ")
      
      if (mode == "integrated" && !grepl("[*|:]", design_str)) {
        warning(
          "In 'integrated' mode, the design formula should contain an ",
          "interaction term (e.g., '*', ':') between condition and time. ",
          "Otherwise, model coefficients for non-reference levels may be ",
          "missing."
        )
      }
      
      if (mode == "isolated" && grepl("[*|:]", design_str)) {
        warning(
          "In 'isolated' mode, the design formula should not contain ",
          "interaction terms (e.g., '*', ':'). These are not required when ",
          "each condition is modeled separately."
        )
      }
      
      message(
        "\nMake sure that the design formula contains no interaction ",
        "between the condition and time for mode == isolated, and that ",
        "it contains an interaction for mode == integrated. Otherwise, ",
        "you will get an uncaught error of 'coefficients not estimable' or ",
        "'subscript out of bounds' (or any other hard to understand error)."
      )
      
    },


    #' Check Multiple Designs and Metas
    #' 
    #' @noRd
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
    #' @noRd
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
    #' @return Returns `NULL` if any required arguments are mising, otherwise,
    #'         called for side effects.
    #'
    check_spline_params = function() {
      spline_params <- self$args$spline_params
      mode <- self$args[["mode"]]
      meta <- self$args$meta
      condition <- self$args$condition

      required_args <- list(
        spline_params,
        mode,
        meta,
        condition
      )

      if (any(vapply(required_args, is.null, logical(1)))) {
        return(NULL)
      }

      self$check_spline_params_generally(spline_params)
      self$check_spline_params_mode_dependent(
        spline_params,
        mode,
        meta,
        condition
      )
    },


    #' Check Spline Test Configurations
    #' 
    #' @noRd
    #'
    #' @describeIn InputControl
    #' This method verifies the spline test configurations and associated
    #' metadata
    #' within the object's arguments. It performs a series of checks on the
    #' configurations, including column verification, spline type validation,
    #' and ensuring that the degrees of freedom (dof) are within acceptable
    #' ranges.
    #'
    #' @param spline_test_configs A configuration object for spline tests.
    #' @param metas A list of metadata corresponding to the data matrices.
    #'
    #' @return Returns `NULL` if any required arguments are mising, otherwise,
    #'         called for side effects.
    #'
    check_spline_test_configs = function() {
      spline_test_configs <- self$args$spline_test_configs
      metas <- self$args$metas

      required_args <- list(spline_test_configs, metas)

      if (any(vapply(required_args, is.null, logical(1)))) {
        return(NULL)
      }

      self$check_columns_spline_test_configs(spline_test_configs)

      self$check_spline_type_column(spline_test_configs)

      self$check_spline_type_params(spline_test_configs)

      self$check_max_and_min_dof(
        spline_test_configs,
        metas
      )
    },


    #' Check Limma Top Tables Structure
    #' 
    #' @noRd
    #'
    #' @description
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

      required_names <- c(
        "time_effect",
        "avrg_diff_conditions",
        "interaction_condition_time"
      )

      if (!all(names(top_tables) %in% required_names) ||
        length(top_tables) != 3) {
        stop("The list must contain exactly three named elements:
         'time_effect', 'avrg_diff_conditions',
         'interaction_condition_time'", call. = FALSE)
      }

      expected_cols1_3 <- c(
        "X1", "X2", "AveExpr", "F", "P.Value", "adj.P.Val",
        "feature_nr", "feature_names", "intercept"
      )

      expected_cols2 <- c(
        "logFC", "AveExpr", "t", "P.Value", "adj.P.Val",
        "B", "feature_nr", "feature_names", "intercept"
      )

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
    #' @noRd
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
      thresh_timeeffect <- self$args[["adj_pthresholds"]]
      thresh_avrdiff    <- self$args[["adj_pthresh_avrg_diff_conditions"]]
      thresh_interact   <- self$args[["adj_pthresh_interaction_condition_time"]]
      
      required_args <- list(
        adj_pthresholds                         = thresh_timeeffect,
        adj_pthresh_avrg_diff_conditions        = thresh_avrdiff,
        adj_pthresh_interaction_condition_time  = thresh_interact
      )
      
      # If any is missing, do nothing (preserve your current behavior)
      if (any(vapply(required_args, is.null, logical(1)))) {
        return(NULL)
      }
      
      validate_vec <- function(name, x) {
        if (!is.numeric(x) || !is.atomic(x) || !is.null(dim(x))) {
          stop(
            paste0("'", name, "' must be a numeric vector."),
            call. = FALSE
          )
        }
        
        if (length(x) == 0) {
          stop(
            paste0("'", name, "' must have length >= 1."),
            call. = FALSE
          )
        }
        
        if (any(!is.finite(x))) {
          idx <- which(!is.finite(x))
          stop(
            paste0(
              "'", name, "' must not contain NA/NaN/Inf. ",
              "Offending indices: ", paste(idx, collapse = ", "),
              ". Values: ", paste(as.character(x[idx]), collapse = ", "), "."
            ),
            call. = FALSE
          )
        }
        
        if (any(x <= 0)) {
          idx <- which(x <= 0)
          stop(
            paste0(
              "'", name, "' must have all elements > 0. ",
              "Offending indices: ", paste(idx, collapse = ", "),
              ". Values: ", paste(x[idx], collapse = ", "), "."
            ),
            call. = FALSE
          )
        }
        
        if (any(x >= 1)) {
          idx <- which(x >= 1)
          stop(
            paste0(
              "'", name, "' must have all elements < 1. ",
              "Offending indices: ", paste(idx, collapse = ", "),
              ". Values: ", paste(x[idx], collapse = ", "), "."
            ),
            call. = FALSE
          )
        }
        
        TRUE
      }
      
      # Validate all three
      invisible(lapply(names(required_args), function(nm) {
        validate_vec(nm, required_args[[nm]])
      }))
      
      return(TRUE)
    },


    #' Check adjusted p-value thresholds for limma category 2 and 3
    #' 
    #' @noRd
    #'
    #' @description
    #' This function checks that both adjusted p-value thresholds for
    #' average difference conditions and interaction condition time are
    #' non-null, floats, and in the range [0, 1].
    #'
    #' @return
    #' `NULL` if either argument is `NULL` or invalid.
    #' Otherwise, no return value (assumed valid inputs).
    #'
    check_adj_pthresh_limma_category_2_3 = function() {
      adj_pthresh_avrg_diff_conditions <-
        self$args[["adj_pthresh_avrg_diff_conditions"]]

      adj_pthresh_interaction_condition_time <-
        self$args[["adj_pthresh_interaction_condition_time"]]

      required_args <- list(
        adj_pthresh_avrg_diff_conditions,
        adj_pthresh_interaction_condition_time
      )

      if (any(vapply(required_args, is.null, logical(1)))) {
        return(NULL)
      }

      # Check that both arguments are numeric and in the range [0, 1]
      if (!all(vapply(required_args, function(arg) {
        is.numeric(arg) && arg >= 0 && arg <= 1
      }, logical(1)))) {
        stop_call_false(
          "Both adj_pthresh_avrg_diff_conditions and",
          "adj_pthresh_interaction_condition_time must",
          "be floats between 0 and 1."
        )
      }
    },


    #' Check Clusters
    #'
    #' @noRd
    #'
    #' @description
    #' Verifies the `nr_clusters` argument provided to the object.
    #'
    #' @details
    #' * `nr_clusters` must be a list with exactly two named elements whose
    #'   names match the two unique levels in `meta[[condition]]`.
    #' * Each element must be a positive integer (e.g. 3) or a vector of
    #'   positive integers (e.g. 2:6).
    #'   
    check_clusters = function() {
      nr_clusters <- self$args[["nr_clusters"]]
      meta        <- self$args[["meta"]]
      condition   <- self$args[["condition"]]
      
      # If any required field is missing, exit quietly
      if (any(vapply(list(nr_clusters, meta, condition), is.null, logical(1))))
        return(NULL)
      
      # --- Extract and validate condition levels (must be exactly 2) ---
      cond_vals <- meta[[condition]]
      if (is.null(cond_vals)) {
        stop_call_false("Column '", condition, "' not found in 'meta'.")
      }
      # Coerce to character to avoid factor level surprises
      cond_levels <- sort(unique(as.character(cond_vals)))
      n_levels <- length(cond_levels)
      if (n_levels != 2L) {
        stop_call_false(
          "Column '", condition, "' in 'meta' must have exactly 2 unique ",
          "values, but has ", n_levels, " (",
          paste(cond_levels, collapse = ", "), ")."
        )
      }
      
      # --- Top-level type, length, and names of nr_clusters ---
      if (!is.list(nr_clusters)) {
        stop_call_false(
          "'nr_clusters' must be a list with exactly two named elements."
        )
      }
      if (length(nr_clusters) != 2L) {
        stop_call_false(
          "'nr_clusters' must have length 2, but has length ",
          length(nr_clusters), "."
        )
      }
      nms <- names(nr_clusters)
      if (is.null(nms) || any(is.na(nms)) || any(nms == "")) {
        stop_call_false(
          "'nr_clusters' must be a named list; names must be non-empty."
        )
      }
      # Names must match the two condition levels exactly (order irrelevant)
      if (!setequal(nms, cond_levels)) {
        stop_call_false(
          "Names of 'nr_clusters' must match levels in '", condition, "'. ",
          "Expected {", paste(cond_levels, collapse = ", "), "}, got {",
          paste(nms, collapse = ", "), "}."
        )
      }
      
      # --- Element-wise validation ---
      for (i in seq_along(nr_clusters)) {
        vec <- nr_clusters[[i]]
        if (!is.numeric(vec) || any(vec != as.integer(vec))) {
          stop_call_false(
            "Element '", nms[i],
            "' of 'nr_clusters' must be an integer vector (e.g. 3 or 2:6)."
          )
        }
        if (any(vec <= 0)) {
          stop_call_false(
            "All integers inside 'nr_clusters' must be positive ",
            "(problem at element '", nms[i], "')."
          )
        }
        if (max(vec) > 100) {
          stop_call_false(
            "Element '", nms[i], "' of 'nr_clusters' has a value (",
            max(vec), ") greater than the allowed maximum of 100."
          )
        }
      }
      
      invisible(TRUE)
    },


    #' Check Plot Info
    #'
    #' @noRd
    #'
    #' @description
    #' Validate `plot_info` for treatment overlays. Require that
    #' `treatment_labels` and `treatment_timepoints` are both absent or both
    #' present. When present, each must be a named list. The names must be a
    #' subset of the values in `meta[[condition]]`. Each `treatment_labels`
    #' value must be a single character string. Each `treatment_timepoints`
    #' value must be a single numeric (integer or double).
    #'
    #' @details
    #' The method performs the following checks:
    #'
    #' * Ensure `plot_info`, `meta`, and `condition` exist; otherwise return
    #'   `NULL`.
    #' * If both treatment lists are absent, validation passes.
    #' * If exactly one of `treatment_labels` or `treatment_timepoints` is
    #'   provided, raise an error.
    #' * Both treatment objects must be lists with non-empty, unique names.
    #' * The two lists must have identical name sets and equal length.
    #' * All names must exist among values of `meta[[condition]]`.
    #' * Each element of `treatment_labels` must be length-1, non-NA char.
    #' * Each element of `treatment_timepoints` must be length-1, non-NA
    #'   numeric.
    #'
    #' If any of these checks fail, raise an informative error.
    #'
    #' @return
    #' `NULL` if `plot_info` is not provided. Otherwise, return
    #' `invisible(TRUE)` on success; errors are thrown on validation failures.
    #' 
    check_plot_info = function() {
      plot_info <- self$args$plot_info
      meta <- self$args$meta
      condition_column <- self$args[["condition"]]
      
      required <- list(plot_info, meta, condition_column)
      if (any(vapply(required, is.null, logical(1)))) {
        return(NULL)
      }
      
      tl <- plot_info$treatment_labels
      tt <- plot_info$treatment_timepoints
      
      if (is.null(tl) && is.null(tt)) {
        return(invisible(TRUE))
      }
      
      if (xor(is.null(tl), is.null(tt))) {
        stop(
          "Provide both treatment_labels and treatment_timepoints, or neither.",
          call. = FALSE
        )
      }
      
      if (!is.list(tl) || !is.list(tt)) {
        stop(
          "treatment_labels and treatment_timepoints must both be lists.",
          call. = FALSE
        )
      }
      
      nms_l <- names(tl)
      nms_t <- names(tt)
      
      if (is.null(nms_l) || is.null(nms_t) ||
          any(is.na(nms_l)) || any(is.na(nms_t)) ||
          any(nchar(nms_l) == 0) || any(nchar(nms_t) == 0)) {
        stop(
          paste(
            "Both treatment_labels and treatment_timepoints must have",
            "named elements with non-empty names."
          ),
          call. = FALSE
        )
      }
      
      if (any(duplicated(nms_l)) || any(duplicated(nms_t))) {
        stop(
          "Names in treatment lists must be unique.",
          call. = FALSE
        )
      }
      
      allowed <- unique(as.character(meta[[condition_column]]))
      allowed <- c(allowed, "double_spline_plots")
      
      invalid_l <- setdiff(nms_l, allowed)
      invalid_t <- setdiff(nms_t, allowed)
      
      if (length(invalid_l) > 0L || length(invalid_t) > 0L) {
        stop(
          paste0(
            "Invalid names in treatment lists. Allowed names are values of `",
            condition_column, "` and 'double_spline_plots'."
          ),
          call. = FALSE
        )
      }
      
      if (!setequal(nms_l, nms_t) || length(tl) != length(tt)) {
        stop(
          "Both lists must have identical name sets and the same length.",
          call. = FALSE
        )
      }
      
      ok_labels <- vapply(
        tl,
        function(x) is.character(x) && length(x) == 1L && !is.na(x),
        logical(1)
      )
      if (!all(ok_labels)) {
        stop(
          "Each treatment_labels element must be a single string.",
          call. = FALSE
        )
      }
      
      ok_tps <- vapply(
        tt,
        function(x) is.numeric(x) && length(x) == 1L && !is.na(x),
        logical(1)
      )
      if (!all(ok_tps)) {
        stop(
          "Each treatment_timepoints element must be a single numeric value.",
          call. = FALSE
        )
      }
      
      invisible(TRUE)
    },
  

    #' Check plot options
    #'
    #' @noRd
    #'
    #' @description
    #' This method checks if the `plot_options` list contains the required
    #' elements
    #' `meta_replicate_column` and `cluster_heatmap_columns`. It validates that
    #' `cluster_heatmap_columns` is either TRUE or FALSE, and that
    #' `meta_replicate_column` is a valid column name in the `meta` dataframe.
    #' If the checks fail, the script stops with an error message.
    #'
    check_plot_options = function() {
      plot_options <- self$args[["plot_options"]]
      meta <- self$args[["meta"]]

      required_args <- list(
        plot_options,
        meta
      )

      # Check if any required arguments are NULL
      if (any(vapply(required_args, is.null, logical(1)))) {
        return(NULL)
      }

      # Ensure at least one of the required elements is present in plot_options
      if (!any(c(
        "meta_replicate_column",
        "cluster_heatmap_columns"
        ) %in% names(plot_options))) {
        stop_call_false(
        "At least one of 'meta_replicate_column' or 'cluster_heatmap_columns' 
        must be present in plot_options"
        )
      }

      # Check if meta_replicate_column is present, and if so, validate it
      if ("meta_replicate_column" %in% names(plot_options)) {
        if (!is.character(plot_options[["meta_replicate_column"]]) ||
          length(plot_options[["meta_replicate_column"]]) != 1) {
          stop_call_false("'meta_replicate_column' must be a single string")
        }

        if (!plot_options[["meta_replicate_column"]] %in% colnames(meta)) {
          stop_call_false(
            "The value of 'meta_replicate_column' does not exist in",
            "the meta dataframe"
          )
        }
      }

      # Check if cluster_heatmap_columns is present, and if so, validate it
      if ("cluster_heatmap_columns" %in% names(plot_options)) {
        if (!is.logical(plot_options[["cluster_heatmap_columns"]]) ||
          length(plot_options[["cluster_heatmap_columns"]]) != 1) {
          stop_call_false("'cluster_heatmap_columns' must be TRUE or FALSE")
        }
      }
    },
    
    
    #' Check Raw Data Validity
    #' 
    #' @noRd
    #'
    #' @description
    #' This function verifies the `raw_data` argument to ensure that it is a 
    #' numeric matrix with the same dimensions as `data`. NA values are permitted 
    #' within `raw_data`. If `raw_data` is not a numeric matrix or does not match 
    #' the dimensions of `data`, the function stops execution with an error.
    #'
    #' @return Returns NULL if either `data` or `raw_data` is not provided or if 
    #' all checks pass. Stops execution and raises an error if `raw_data` does 
    #' not meet the criteria.
    #'
    check_raw_data = function() {
      data <- self$args[["data"]]
      raw_data <- self$args[["raw_data"]]
      
      required_args <- list(data, raw_data)
      
      if (any(vapply(required_args, is.null, logical(1)))) {
        return(NULL)
      }
      
      # Check if raw_data is a numeric matrix (NA values allowed)
      if (!is.matrix(raw_data) || !is.numeric(raw_data)) {
        stop_call_false(
          "'raw_data' must be a numeric matrix (NA values allowed)"
        )
      }
      
      # Check if raw_data has the same dimensions as data
      if (!all(dim(raw_data) == dim(data))) {
        stop_call_false(
          "'raw_data' must have the same dimensions as 'data'"
        )
      }
    },


    #' Check and Create Report Directory
    #' 
    #' @noRd
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
      report_dir <- self$args[["report_dir"]]

      if (is.null(report_dir)) {
        # some functions just have a different name.
        report_dir <- self$args$output_dir
      }

      required_args <- list(report_dir)

      if (any(vapply(required_args, is.null, logical(1)))) {
        return(NULL)
      }

      # Attempt to create the directory if it does not exist
      if (!file.exists(report_dir)) {
        tryCatch(
          {
            dir.create(report_dir, recursive = TRUE)
          },
          warning = function(w) {
            stop_call_false(
              sprintf(
                "Warning occurred while creating the directory: %s",
                w$message
              )
            )
          },
          error = function(e) {
            stop_call_false(
              sprintf(
                "Error occurred while creating the directory: %s",
                e$message
              )
            )
          }
        )
      }

      # Verify that the directory exists and is a directory
      if (!file.exists(report_dir) || !file.info(report_dir)$isdir) {
        stop_call_false(
          sprintf(
            "The specified path is not a valid directory: %s",
            report_dir
          )
        )
      }

      return(TRUE)
    },


    #' Check Genes Validity
    #' 
    #' @noRd
    #'
    #' @description
    #' This function checks the validity of the `data` and `genes` arguments
    #' within the `self$args` list. It ensures that `genes` is a character
    #' vector,
    #' that neither `data` nor `genes` is `NULL`, and that the length of `genes`
    #' matches the number of rows in `data`.
    #'
    #' @return Returns `TRUE` if all checks pass. Returns `NULL` if any required
    #' arguments are `NULL`. Throws an error if `genes` is not a character
    #' vector
    #' or if the length of `genes` does not match the number of rows in `data`.
    #'
    check_genes = function() {
      data <- self$args$data
      genes <- self$args$genes

      required_args <- list(data, genes)

      if (any(vapply(required_args, is.null, logical(1)))) {
        return(NULL)
      }

      if (!is.character(genes)) {
        stop(paste0("genes must be a character vector"), call. = FALSE)
      }

      if (length(genes) != nrow(data)) {
        stop(paste0("length(genes) must be equal to nrow(data)"), call. = FALSE)
      }
    },


    #' Check p-Adjustment Method
    #' 
    #' @noRd
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
      padjust_method <- self[["args"]][["padjust_method"]]

      required_args <- list(padjust_method)

      if (any(vapply(required_args, is.null, logical(1)))) {
        return(NULL)
      }

      supported_methods <- stats::p.adjust.methods
      if (!(is.character(padjust_method) &&
        padjust_method %in% supported_methods)) {
        stop(
          sprintf(paste(
            "padjust_method must be a character and one of the",
            "supported methods (%s).",
            paste(supported_methods, collapse = ", ")
          )),
          call. = FALSE
        )
      }
    },


    #' Check Report Information
    #' 
    #' @noRd
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
      report_info <- self$args[["report_info"]]

      required_args <- list(report_info)

      if (any(vapply(required_args, is.null, logical(1)))) {
        return(NULL)
      }

      mandatory_fields <- c(
        "omics_data_type",
        "data_description",
        "data_collection_date",
        "analyst_name",
        "contact_info",
        "project_name"
      )

      all_fields <- c(
        mandatory_fields,
        "method_description",
        "results_summary",
        "conclusions"
      )

      # Check if report_info is a named list
      if (!is.list(report_info) || is.null(names(report_info))) {
        stop_call_false("report_info must be a named list.")
      }

      # Check if all values in report_info are strings
      non_string_fields <- vapply(
        report_info,
        function(x) !is.character(x),
        logical(1)
      )

      if (any(non_string_fields)) {
        invalid_fields <- names(report_info)[non_string_fields]
        stop_call_false(paste(
          "The following fields in report_info must be strings:",
          paste(invalid_fields, collapse = ", ")
        ))
      }

      # Check if all mandatory fields are present
      missing_fields <- setdiff(mandatory_fields, names(report_info))
      if (length(missing_fields) > 0) {
        stop_call_false(paste(
          "Missing mandatory fields in report_info:",
          paste(missing_fields, collapse = ", ")
        ))
      }

      # Check if there are any extra fields not in all_fields
      extra_fields <- setdiff(names(report_info), all_fields)
      if (length(extra_fields) > 0) {
        stop_call_false(paste(
          "The following fields in report_info are not recognized:",
          paste(extra_fields, collapse = ", ")
        ))
      }

      # Check omics_data_type format
      if (!grepl("^[a-zA-Z_]+$", report_info[["omics_data_type"]])) {
        stop_call_false(paste(
          "The 'omics_data_type' field must contain only alphabetic",
          "letters and underscores."
        ))
      }

      excluded_fields <- c(
        "data_description",
        "method_description",
        "results summary",
        "conclusions"
      )
      excluded_limit <- 700

      check_long_fields <- function(data,
                                    excluded_fields,
                                    excluded_limit) {
        long_fields <- vapply(data, function(x) {
          if (any(names(data) %in% excluded_fields)) {
            any(nchar(x) > excluded_limit)
          } else {
            any(nchar(x) > 70)
          }
        }, logical(1))
        return(long_fields)
      }

      # Check if any field exceeds 70 characters
      long_fields <- check_long_fields(
        report_info,
        excluded_fields,
        excluded_limit
      )

      if (any(long_fields)) {
        too_long_fields <- names(report_info)[long_fields]
        stop(paste("The following fields have strings exceeding 70 characters:",
          paste(too_long_fields, collapse = ", "),
          sep = "\n"
        ), call. = FALSE)
      }

      return(TRUE)
    },


    #' Check Feature Name Columns
    #' 
    #' @noRd
    #'
    #' @description
    #' This function checks whether all elements of `feature_name_columns` are
    #' characters of length 1 and whether they are valid column names in the
    #' `annotation` data frame.
    #'
    #' @return Returns `NULL` if any required arguments are missing. Throws
    #' an error
    #' if any element of `feature_name_columns` is not a character of length 1
    #' or if
    #' any element is not a column name in `annotation`. Returns `TRUE` if all
    #'  checks
    #' pass.
    #'
    #' @details
    #' This function performs the following checks:
    #' 1. Ensures `feature_name_columns` and `annotation` are not `NULL`.
    #' 2. Verifies that each element in `feature_name_columns` is a character
    #'  with
    #'    a length of 1.
    #' 3. Checks that all elements of `feature_name_columns` are valid column
    #' names in the `annotation` data frame.
    #'
    check_feature_name_columns = function() {
      feature_name_columns <- self$args[["feature_name_columns"]]
      annotation <- self$args[["annotation"]]

      required_args <- list(
        feature_name_columns,
        annotation
      )

      if (any(vapply(required_args, is.null, logical(1)))) {
        return(NULL)
      }

      # Check if every element in feature_name_columns is a character
      # of length 1
      if (!all(vapply(
        feature_name_columns,
        function(x) is.character(x) && length(x) == 1,
        logical(1)
      ))
      ) {
        stop_call_false(
          paste(
            "All elements of feature_name_columns must be characters",
            "with length 1."
          )
        )
      }

      # Check if all elements of feature_name_columns are column names
      # in annotation
      if (!all(feature_name_columns %in% colnames(annotation))) {
        stop_call_false(
          paste(
            "All elements of feature_name_columns must be column names",
            "in annotation."
          )
        )
      }
    },


    #' Check Report
    #' 
    #' @noRd
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
      report <- self$args[["report"]]

      required_args <- list(report)

      if (any(vapply(required_args, is.null, logical(1)))) {
        return(NULL)
      }

      if (!rlang::is_bool(report)) {
        stop_call_false("report must be either Boolean TRUE or FALSE")
      }
    },
    
    
    #' Check Support Argument
    #'
    #' @noRd
    #'
    #' @description
    #' This method validates the `support` argument to ensure it is a properly
    #' formatted, non-negative integer. This argument is typically used to 
    #' control thresholds for data completeness or filtering criteria.
    #'
    #' @details
    #' The method performs the following checks:
    #'
    #' * Ensures that `support` is provided.
    #' * Confirms that `support` is a numeric scalar.
    #' * Verifies that `support` is an integer (not a float or character).
    #' * Checks that `support` is non-negative.
    #'
    #' If any of these checks fail, an informative error message is returned.
    #'
    #' @return NULL if `support` is not provided. Otherwise, performs
    #' checks and potentially raises errors if checks fail.
    #' 
    check_support = function() {
      support <- self$args[["support"]]
      
      if (is.null(support)) {
        return(NULL)
      }
      
      if (!is.numeric(support) 
          || length(support) != 1 
          || support != as.integer(support)) {
        stop(
          "support must be a single integer value!",
          call. = FALSE
        )
      }
      
      if (support < 0) {
        stop(
          "support must be non-negative!",
          call. = FALSE
        )
      }
    }
  )
)


# Level2Functions class --------------------------------------------------------


#' Level2Functions: A class providing level 2 functionalities
#' 
#' @noRd
#'
#' @description
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
    #'  data for error messages. Default is NA.
    #'
    #' @return Returns TRUE if all checks pass. Stops execution and returns an
    #' error
    #' message if any check fails.
    #'
    check_data = function(
        data,
        data_meta_index = NULL
        ) {
      
      all_zero <- function(x) {
        apply(x, 1, function(row) all(row == 0))
      }
      
      # Check for all-zero rows
      zero_rows <- which(all_zero(data))
      if (length(zero_rows) > 0) {
        row_names <- rownames(data)[zero_rows]
        stop_call_false(
          self$create_error_message(
            paste0(
              "Data must not contain rows where all values are zero! ",
              "Offending rows: ",
              paste(ifelse(
                is.null(row_names),
                zero_rows,
                row_names
                ),
                collapse = ", ")
            ),
            data_meta_index
          )
        )
      }
      
      # Check for all-NA rows
      na_rows <- which(apply(data, 1, function(row) all(is.na(row))))
      if (length(na_rows) > 0) {
        row_names <- rownames(data)[na_rows]
        stop_call_false(
          self$create_error_message(
            paste0(
              "Data must not contain rows where all values are NA! ",
              "Offending rows: ",
              paste(ifelse(
                is.null(row_names),
                na_rows,
                row_names
                ),
                collapse = ", ")
            ),
            data_meta_index
          )
        )
      }
      
      # Check for all-NA columns
      na_cols <- which(apply(data, 2, function(col) all(is.na(col))))
      if (length(na_cols) > 0) {
        col_names <- colnames(data)[na_cols]
        stop_call_false(
          self$create_error_message(
            paste0(
              "Data must not contain columns where all values are NA! ",
              "Offending columns: ",
              paste(ifelse(
                is.null(col_names),
                na_cols,
                col_names
                ),
                collapse = ", ")
            ),
            data_meta_index
          )
        )
      }

      # Ensure data is a numeric matrix
      if (!inherits(data, "matrix") || !is.numeric(data)) {
        stop_call_false(paste0(
          "Input 'data' must be a numeric matrix (class must include 'matrix'",
          "and typeof 'double'). ",
          "Actual class: ", paste(class(data), collapse = ", "),
          "; typeof: ", typeof(data), ". ",
          data_meta_index
        ))
      }
      
      # Check for missing values (only warn instead of stopping execution)
      if (any(is.na(data))) {
        message(
          self$create_error_message(
            paste(
              "Warning: The data contains missing values (NA).",
              "Ensure that imputation or handling of missing values",
              "aligns with your analysis requirements. Note that limma can",
              " handle missing values (it just ignores them), and therefore",
              " also SplineOmics can handle them."
            ),
            data_meta_index
          )
        )
      }
      
      # Check for non-negative values, allowing NAs
      if (any(data < 0, na.rm = TRUE)) {  # Ignore NA values in the check
        message(
          self$create_error_message(
            paste(
              "Hint: The data contains negative values. This may occur if the",
              "data has been transformed (e.g., log-transformed or normalized)",
              "and is valid in such cases. Ensure that the data preprocessing",
              "aligns with your analysis requirements."
            ),
            data_meta_index
          )
        )
      }
      
      # Check if row headers (rownames) are present and non-null
      row_headers <- rownames(data)
      if (is.null(row_headers)) {
        stop_call_false(
          self$create_error_message(
            "The data matrix must have row headers!",
            data_meta_index
          )
        )
      }
      
      # Check if all row headers are unique
      if (any(duplicated(row_headers))) {
        dupes <- unique(row_headers[duplicated(row_headers)])
        stop(
          self$create_error_message(
            paste0(
              "The data matrix has duplicated row headers: ",
              paste(dupes, collapse = ", "), ". ",
              "All row names must be unique."
            ),
            data_meta_index
          ),
          call. = FALSE
        )
      }

      return(TRUE)
    },


    #' Check Metadata
    #' 
    #' @noRd
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
    check_meta = function(
        meta,
        condition,
        meta_batch_column = NULL,
        meta_batch2_column = NULL,
        data_meta_index = NULL
        ) {
      
      if (!is.data.frame(meta) ||
        !"Time" %in% names(meta) ||
        !is.numeric(meta[["Time"]])) {
        stop_call_false(
          self$create_error_message(
            paste("meta must be a dataframe with the numeric column Time"),
            data_meta_index
          )
        )
      }

      if (any(is.na(meta))) {
        stop_call_false(
          self$create_error_message(
            "meta must not contain missing values.",
            data_meta_index
          )
        )
      }

      # Check if condition is a single character
      if (!is.character(condition) || length(condition) != 1) {
        stop_call_false("'condition' must be a single character")
      }

      # Check if condition is a column in the meta dataframe
      if (!condition %in% colnames(meta)) {
        stop_call_false(self$create_error_message(
          sprintf(
            "The condition '%s' is not a %s",
            condition,
            paste("column in meta")
          ),
          data_meta_index
        ))
      }
      
      # Ensure the condition column has at most two levels
      cond_vals <- meta[[condition]]
      # Count distinct observed values (ignore unused factor levels)
      n_levels <- length(unique(if (is.factor(cond_vals))
        droplevels(cond_vals) else cond_vals))
      if (n_levels > 2) {
        levs <- sort(unique(as.character(cond_vals)))
        msg <- sprintf(
          "The condition column '%s' must have at most two levels,
          but found %d: %s",
          condition,
          n_levels,
          paste(levs, collapse = ", ")
        )
        stop_call_false(self$create_error_message(msg, data_meta_index))
      }

      # Check if the factor column is of appropriate type
      if (!is.factor(meta[[condition]]) &&
        !is.character(meta[[condition]])) {
        stop(
          self$create_error_message(
            sprintf("The factor column '%s' must be of type
                 factor or character.", condition),
            data_meta_index
          ),
          call. = FALSE
        )
      }

      if (!is.character(meta_batch2_column)) {
        meta_batch2_column <- NULL
      }

      if (is.null(meta_batch_column) && !is.null(meta_batch2_column)) {
        stop(paste(
          "For removing the batch effect, batch2 can only be used when",
          "batch is used!"
        ), call. = FALSE)
      }

      if (!is.null(meta_batch_column)) {
        if (is.character(meta_batch_column)) {
          if (meta_batch_column == "Time" || meta_batch_column == condition) {
            stop(paste("meta_batch_column must not be == 'Time' or", condition),
              call. = FALSE
            )
          }

          self$check_batch_column(
            meta,
            meta_batch_column,
            data_meta_index
          )
        } else {
          stop(
            "meta_batch_column must be a character",
            call. = FALSE
          )
        }
      }

      if (!is.null(meta_batch2_column)) {
        if (is.character(meta_batch2_column)) {
          if (meta_batch2_column == "Time" || meta_batch2_column == condition) {
            stop(paste("meta_batch2_column must not be == 'Time' or", condition),
              call. = FALSE
            )
          }

          if (meta_batch_column == meta_batch2_column) {
            stop(
              paste(
                "meta_batch_column must not be equal to",
                "meta_batch2_column"
              ),
              call. = FALSE
            )
          }

          self$check_batch_column(
            meta,
            meta_batch2_column,
            data_meta_index
          )
        } else {
          stop(
            "meta_batch2_column must be a character",
            call. = FALSE
          )
        }
      }
      
      # Check if batch column is blocked (only one level)
      if (!is.null(meta_batch_column)) {
        colname <- as.character(meta_batch_column)
        n_levels <- length(unique(meta[[colname]]))
        if (n_levels <= 1) {
          stop(
            paste0(
              "meta_batch_column '", colname, "' is blocked (only ",
              n_levels, " level). Cannot be used."
            ),
            call. = FALSE
          )
        }
      }
      
      if (!is.null(meta_batch2_column)) {
        colname <- as.character(meta_batch2_column)
        n_levels <- length(unique(meta[[colname]]))
        if (n_levels <= 1) {
          stop(
            paste0(
              "meta_batch2_column '", colname, "' is blocked (only ",
              n_levels, " level). Cannot be used."
            ),
            call. = FALSE
          )
        }
      }
      
      return(TRUE)
    },


    #' Check Dataframe
    #' 
    #' @noRd
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
        feature_nr = "integer",
        feature_names = "character",
        intercept = "numeric"
      )

      # Check if all required columns are present
      missing_columns <- setdiff(names(required_columns), names(df))
      if (length(missing_columns) > 0) {
        stop(
          paste(
            "Missing columns in top_table:",
            paste(missing_columns, collapse = ", ")
          ),
          call. = FALSE
        )
      }

      # Check if columns have the correct type
      for (col in names(required_columns)) {
        if (!inherits(df[[col]], required_columns[[col]])) {
          stop(
            paste(
              "top_table column", col, "must be of type",
              required_columns[[col]]
            ),
            call. = FALSE
          )
        }
      }

      return(TRUE)
    },


    #' Check Spline Parameters Generally
    #' 
    #' @noRd
    #'
    #' @description
    #' Validates the general structure and contents of spline parameters.
    #'
    #' @param spline_params A list of spline parameters.
    #'
    #' @return No return value, called for side effects.
    #'
    check_spline_params_generally = function(spline_params) {
      allowed_fields <- c("spline_type", "degree", "dof")

      if (!all(names(spline_params) %in% allowed_fields)) {
        stop(
          paste(
            "spline_params contains invalid fields. Only 'spline_type',",
            "'degree', and 'dof' are allowed."
          ),
          call. = FALSE
        )
      }

      # Check if spline_type exists and contains valid values
      if ("spline_type" %in% names(spline_params)) {
        if (!all(spline_params$spline_type %in% c("b", "n"))) {
          stop(
            paste(
              "Elements of spline_type must be either 'b' for B-splines",
              "or 'n' for natural cubic splines."
            ),
            call. = FALSE
          )
        }
      } else {
        stop("spline_type is missing in spline_params.", call. = FALSE)
      }

      # Check if degree exists and is an integer vector for B-splines,
      # and NA for natural splines
      if ("degree" %in% names(spline_params)) {
        for (i in seq_along(spline_params$spline_type)) {
          if (spline_params$spline_type[i] == "b" &&
            (!is.integer(spline_params$degree[i]) ||
              is.na(spline_params$degree[i]))) {
            stop(
              paste(
                "Degree must be specified as an integer for",
                "B-splines in spline_params."
              ),
              call. = FALSE
            )
          }
          if (spline_params$spline_type[i] == "n" &&
            !is.na(spline_params$degree[i])) {
            stop(
              paste(
                "Degree must be NA for natural cubic",
                "splines in spline_params."
              ),
              call. = FALSE
            )
          }
        }
      } else if (all(spline_params$spline_type == "b")) {
        stop("degree is missing in spline_params.", call. = FALSE)
      }

      # Check if dof exists and is an integer vector
      if ("dof" %in% names(spline_params)) {
        if (!all(spline_params$dof == as.integer(spline_params$dof))) {
          stop_call_false("dof must be an integer vector in spline_params.")
        }
        # B-spline DoF must be > 2 - but only if not zero (i.e., not auto)
        for (i in seq_along(spline_params$spline_type)) {
          if (spline_params$spline_type[i] == "b" &&
              spline_params$dof[i] != 0 &&
              spline_params$dof[i] < 3) {
            stop_call_false(
              paste("B-splines require DoF > 2 at index", i,
                    "unless auto-selection (DoF = 0) is used.")
            )
          }
        }
      } else {
        stop("dof is missing in spline_params.", call. = FALSE)
      }
    },


    #' Check Spline Parameters Mode Dependent
    #' 
    #' @noRd
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
    check_spline_params_mode_dependent = function(
        spline_params,
        mode,
        meta,
        condition
        ) {

      if (mode == "integrated") {
        # every entry in spline_params must be a scalar (length-1)
        if (any(vapply(spline_params,
                       function(x) length(x) != 1,
                       logical(1)))) {
          stop_call_false(paste(
            "All parameters in spline_params must have exactly one element",
            "when mode is 'integrated'.",
            "Different spline parameters for the different levels is not",
            "supported for this mode."
          ))
        }

        # # Additional check for 'knots' and 'bknots' if they exist
        # if ("knots" %in% names(spline_params)) {
        #   # Check if 'knots' is atomic (i.e., treat vectors as a single unit) or NA
        #   if (!(is.atomic(spline_params$knots) || all(is.na(spline_params$knots)))) {
        #     stop(
        #       paste(
        #         "All elements in 'knots' in spline_params must be atomic",
        #         "or NA when mode is 'integrated'.",
        #         "Different spline parameters for different levels are not",
        #         "supported for this mode."
        #       ),
        #       call. = FALSE
        #     )
        #   }
        # }
        #
        # if ("bknots" %in% names(spline_params)) {
        #   # Check if 'bknots' is atomic (i.e., treat vectors as a single unit) or NA
        #   if (!(is.atomic(spline_params$bknots) ||
        #         all(is.na(spline_params$bknots)))) {
        #     stop(
        #       paste(
        #         "All elements in 'bknots' in spline_params must be atomic or",
        #         "NA when mode is 'integrated'.",
        #         "Different spline parameters for different levels are not",
        #         "supported for this mode."
        #       ),
        #       call. = FALSE
        #     )
        #   }
        # }
      } else if (mode == "isolated") {
        num_levels <- length(unique(meta[[condition]]))
        if (any(vapply(spline_params, length, integer(1)) != num_levels)) {
          stop_call_false(paste(
            "Each vector in spline_params must have as many",
            "elements as there are unique elements in the",
            condition,
            "column of meta when mode is 'isolated'."
          ))
        }
        # if ("knots" %in% names(spline_params)) {
        #   if (length(spline_params$knots) != num_levels) {
        #     stop(
        #       paste("'knots' in spline_params must have the same number of",
        #             "elements as there are unique elements in the ", condition,
        #             "column of meta when mode is 'isolated'."),
        #          call. = FALSE)
        #   }
        # }
        # if ("bknots" %in% names(spline_params)) {
        #   if (length(spline_params$bknots) != num_levels) {
        #     stop(
        #       paste("'bknots' in spline_params must have the same number of",
        #             "elements as there are unique elements in the ", condition,
        #             "column of meta when mode is 'isolated'."),
        #          call. = FALSE)
        #   }
        # }
      }
    },


    #' Check Columns in Spline Test Configurations
    #' 
    #' @noRd
    #'
    #' @description
    #' Validates that the spline test configurations contain the required 
    #' columns in the correct order.
    #'
    #' @param spline_test_configs A dataframe containing spline test
    #'  configurations.
    #'
    #' @return No return value, called for side effects.
    #'
    #' @keywords internal
    #'
    check_columns_spline_test_configs = function(spline_test_configs) {
      required_columns <- c(
        "spline_type",
        "degree",
        "dof"
        # "knots",
        # "bknots"
      )

      # Check for exact match of column names and order
      if (!identical(names(spline_test_configs), required_columns)) {
        # Find the missing or extra columns
        missing_columns <- setdiff(required_columns, names(spline_test_configs))
        extra_columns <- setdiff(names(spline_test_configs), required_columns)
        error_message <- "Error: Incorrect columns in dataframe. "

        # Append specific issues to the error message
        if (length(missing_columns) > 0) {
          error_message <- paste0(
            error_message, "Missing columns: ",
            paste(missing_columns, collapse = ", "), ". "
          )
        }
        if (length(extra_columns) > 0) {
          error_message <- paste0(
            error_message, "Extra columns: ",
            paste(extra_columns, collapse = ", "), ". "
          )
        }
        error_message <- paste0(
          error_message,
          "Expected columns in order: ",
          paste(required_columns, collapse = ", "), "."
        )

        stop(error_message)
      }
    },


    #' Check Spline Type Column
    #' 
    #' @noRd
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
    #' @keywords internal
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
          paste(unique(invalid_entries), collapse = ", ")
        )

        stop(error_message)
      }
    },


    #' Check Spline Type Parameters
    #' 
    #' @noRd
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
    #' @keywords internal
    #'
    check_spline_type_params = function(spline_test_configs) {
      for (i in seq_len(nrow(spline_test_configs))) {
        row <- spline_test_configs[i, ]
        switch(as.character(row$spline_type),
          "n" = {
            if (!is.na(row$degree)) {
              stop("degree must be NA for spline_type n")
            }

            # if (!((is.na(row$dof) && !is.na(row$knots)) ||
            #       (!is.na(row$dof) && is.na(row$knots)))) {
            #   stop("Either dof or knots must be NA, but not both, for
            #     spline_type n")
            # }

            if (!is.integer(row$dof)) {
              stop("dof must be an integer when it is not NA for
                    spline_type n")
            }

            # if (!is.na(row$knots) && !is.numeric(row$knots)) {
            #   stop("knots must be a numeric vector when it is not NA for
            #     spline_type n")
            # }
            # if (!is.na(row$bknots) && (!is.numeric(row$bknots) ||
            #                            length(row$bknots) != 2)) {
            #   stop("bknots must be a numeric vector of exactly two elements
            #     or NA for spline_type n")
            # }
          },
          "b" = {
            if (!is.integer(row$degree)) stop("degree must be an integer for
                                               spline_type b")
            if (!is.na(row$dof) && (!is.integer(row$dof) ||
              row$dof < row$degree)) {
              stop("dof must be an integer at least as big as degree for
                    spline_type b")
            }

            # if (!is.na(row$knots) && !is.numeric(row$knots)) {
            #   stop("knots must be a numeric vector when it is not NA for
            #     spline_type b")
            # }
            # if (!is.na(row$bknots) && (!is.numeric(row$bknots) ||
            #                            length(row$bknots) != 2)) {
            #   stop("bknots must be a numeric vector of exactly two elements
            #     or NA for spline_type b")
            # }
          },
          stop("spline_type must be either 'n' or 'b'")
        )
      }
      return(TRUE)
    },


    #' Check Maximum and Minimum Degrees of Freedom
    #' 
    #' @noRd
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
    #' @keywords internal
    #'
    check_max_and_min_dof = function(
        spline_test_configs,
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
    #' @noRd
    #'
    #' @description
    #' This function checks if the columns of a dataframe match the expected
    #' column names and their respective data types.
    #'
    #' @param df A dataframe to check.
    #' @param expected_cols A character vector of expected column names.
    #'
    #' @return This function does not return a value. It stops execution if the
    #' dataframe columns or their classes do not match the expected structure.
    #'
    check_columns = function(
        df,
        expected_cols) {
      actual_cols <- names(df)
      if (!all(expected_cols %in% actual_cols)) {
        stop_call_false("Dataframe columns do not match expected structure")
      }
      expected_classes <- c(
        "numeric",
        "numeric",
        "numeric",
        "numeric",
        "numeric",
        "numeric",
        "integer",
        "character",
        "numeric"
      )
      actual_classes <- vapply(df, class, character(1))
      if (!all(actual_classes == expected_classes)) {
        stop_call_false(
          "Dataframe column classes do not match expected classes"
          )
      }
    }
  )
)



# Level3Functions class --------------------------------------------------------


#' Level3Functions: A class for level 3 utility functions
#' 
#' @noRd
#'
#' @description
#' This class provides methods for creating error messages and checking
#' batch columns.
#'
#' @seealso \code{\link{Level2Functions}}
#'
Level3Functions <- R6::R6Class("Level3Functions",
  inherit = Level4Functions,
  public = list(

    #' Check the structure of a voom object
    #' 
    #' @noRd
    #'
    #' @description
    #' This function checks the structure of a `voom` object to ensure that it
    #' contains
    #' all the expected components and that these components have the correct
    #' types
    #' and dimensions. The function does not check the actual data within the
    #' matrices.
    #'
    #' @param voom_obj A list representing a `voom` object, typically created
    #' by the
    #'   `voom` function from the `limma` package.
    #'
    #' @details
    #' The function verifies that the `voom` object contains the following
    #' components:
    #' - `E`: A matrix of log2-counts per million (logCPM) values.
    #' - `weights`: A matrix of observation-specific weights that matches the
    #' dimensions of `E`.
    #' - `design`: A matrix representing the design matrix used in the linear
    #' modeling,
    #'   with the same number of rows as there are columns in `E`.
    #'
    #' The function also checks for optional components such as:
    #' - `genes`: A data frame of gene annotations.
    #' - `targets`: A data frame of target information.
    #' - `sample.weights`: A numeric vector of sample-specific weights.
    #'
    #' If any of these checks fail, the function stops and reports the issues.
    #' If the structure is valid, a message confirming the validity is printed.
    #'
    #' @return Boolean TRUE or FALSE. However, the function is mostly called for
    #' its side effects, which stop the script if the structure is not valid.
    #'
    check_voom_structure = function(voom_obj) {
      # Initialize a list to collect any issues found
      issues <- list()

      # Check if the input is a list
      if (!is.list(voom_obj)) {
        stop("The input is not a list. A voom object should be a list.")
      }

      # Check for the presence of the expected components
      expected_components <- c("E", "weights", "design")
      for (comp in expected_components) {
        if (!comp %in% names(voom_obj)) {
          issues <- c(issues, paste("Missing component:", comp))
        }
      }

      # Check that 'E' is a matrix
      if ("E" %in% names(voom_obj) && !is.matrix(voom_obj$E)) {
        issues <- c(issues, "'E' should be a matrix.")
      }

      # Check that 'weights' is a matrix and matches the dimensions of 'E'
      if ("weights" %in% names(voom_obj)) {
        if (!is.matrix(voom_obj$weights)) {
          issues <- c(issues, "'weights' should be a matrix.")
        } else if (!all(dim(voom_obj$weights) == dim(voom_obj$E))) {
          issues <- c(
            issues,
            "'weights' dimensions do not match 'E' dimensions."
          )
        }
      }

      # Check that 'design' is a matrix and has the correct number of rows
      if ("design" %in% names(voom_obj)) {
        if (!is.matrix(voom_obj$design)) {
          issues <- c(issues, "'design' should be a matrix.")
        } else if (nrow(voom_obj$design) != ncol(voom_obj$E)) {
          issues <- c(
            issues,
            "'design' matrix should have the same number of rows as the
           number of columns in 'E'."
          )
        }
      }

      # Optionally, check for the presence of other common components
      if ("genes" %in% names(voom_obj) && !is.data.frame(voom_obj$genes)) {
        issues <- c(issues, "'genes' should be a data frame.")
      }

      # Check for the presence of optional components like targets or
      # sample weights
      if ("targets" %in% names(voom_obj) && !is.data.frame(voom_obj$targets)) {
        issues <- c(issues, "'targets' should be a data frame.")
      }

      if ("sample.weights" %in% names(voom_obj) &&
        !is.numeric(voom_obj$sample.weights)) {
        issues <- c(issues, "'sample.weights' should be numeric.")
      }

      # Report results
      if (length(issues) > 0) {
        stop(
          "The voom object failed the structure check:\n",
          paste(
            issues,
            collapse = "\n"
          ),
          call. = FALSE
        )
      } else {
        return(TRUE)
      }
    },


    #' Check Batch Column
    #' 
    #' @noRd
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
    check_batch_column = function(
      meta,
      meta_batch_column,
      data_meta_index
    ) {
      if (!is.null(meta_batch_column) 
          && !(meta_batch_column %in% names(meta))) {
        stop(
          self$create_error_message(
            sprintf(
              "Batch effect column '%s' not found in meta",
              meta_batch_column
            ),
            data_meta_index
          ),
          call. = FALSE
        )
      }
      # else: valid column or NULL → silently continue
    }
  )
)



# Level4Functions class --------------------------------------------------------

#' Level4Functions: A class for level 3 utility functions
#' 
#' @noRd
#'
#' @description
#' This class provides methods for creating error messages and checking
#' batch columns.
#'
#' @seealso \code{\link{Level3Functions}}
#'
Level4Functions <- R6::R6Class("Level4Functions",
  public = list(

    #' Create Error Message
    #' 
    #' @noRd
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
    #' @return Returns a formatted error message string. If an index is 
    #' provided, the message includes the index; otherwise, it returns the 
    #' message as is.
    #'
    create_error_message = function(
        message,
        data_meta_index = NULL) {
      if (!is.null(data_meta_index)) {
        return(sprintf(
          "data/meta pair index %d: %s",
          data_meta_index,
          message
        ))
      } else {
        return(message)
      }
    }
  )
)



# Utility input control functions ----------------------------------------------


#' Check for Required Elements in the SplineOmics Object
#' 
#' @noRd
#'
#' @description
#' This function checks if the given object contains all the required named
#' elements for a specified function type. If any element is missing, it stops
#' the script and provides an informative error message.
#'
#' @param splineomics The object to be checked.
#' @param func_type A string specifying the function type. It can be one of
#' "cluster_hits", "create_limma_report", "run_limma_splines", or "explore_data"
#'
#' @return None. Stops execution if any required element is missing.
#'
check_splineomics_elements <- function(
    splineomics,
    func_type
    ) {
  
  required_elements <- switch(func_type,
    "explore_data" = c(
      "data",
      "meta",
      "condition",
      "report_info"
    ),
    "screen_limma_hyperparams" = c(
      "condition",
      "report_info",
      "padjust_method"
    ),
    "preprocess_rna_seq_data" = c(
      "data",
      "meta",
      "spline_params",
      "design",
      "condition"
    ),
    "run_limma_splines" = c(
      "data",
      "meta",
      "design",
      "mode",
      "condition",
      "spline_params",
      "padjust_method"
    ),
    "find_peaks_valleys" = c(
      "data",
      "meta",
      "meta_batch_column"
    ),
    "create_limma_report" = c(
      "limma_splines_result",
      "report_info"
    ),
    "cluster_hits" = c(
      "data",
      "meta",
      "design",
      "mode",
      "condition",
      "spline_params",
      "limma_splines_result"
    ),
    stop_call_false("Invalid function type provided.")
  )

  missing_elements <- required_elements[
    !vapply(
      splineomics[required_elements],
      function(x) !is.null(x),
      logical(1)
    )
  ]

  if (length(missing_elements) > 0) {
    stop_call_false(paste(
      "The following required elements for the function",
      func_type,
      "were not passed to the SplineOmics object:",
      paste(
        missing_elements,
        collapse = ", "
      ),
      "\nAll required elements for",
      func_type,
      "are:",
      paste(
        required_elements,
        collapse = ", "
      )
    ))
  }
}


#' Check for NULL Elements in Arguments
#' 
#' @noRd
#'
#' @description
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
  null_elements <- names(args)[vapply(args, is.null, logical(1))]

  if (length(null_elements) > 0) {
    stop(
      paste(
        "The following function arguments are NULL:",
        paste(null_elements, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}
