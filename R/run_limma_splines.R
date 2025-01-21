#' run_limma_splines()
#'
#' @description
#' This function performs a limma spline analysis to identify significant
#' time-dependent changes in features (e.g., proteins) within an omics
#' time-series dataset. It evaluates features within each condition level
#' and between levels by comparing average differences and interactions
#' between time and condition.
#'
#' @param splineomics An S3 object of class `SplineOmics` that contains the
#' following elements:
#' \itemize{
#'   \item \code{data}: The matrix of the omics dataset, with the feature
#'   names optionally as row headers.
#'   \item \code{rna_seq_data}: An object containing the preprocessed
#'   RNA-seq data,
#'   such as the output from `limma::voom` or a similar preprocessing pipeline.
#'   \item \code{meta}: A dataframe containing metadata corresponding to the
#'   \code{data}, must include a 'Time' column and the column specified by
#'   \code{condition}.
#'   \item \code{design}: A character string representing the limma design
#'   formula.
#'   \item \code{condition}: A character string specifying the column name
#'   in \code{meta} used to define groups for analysis.
#'   \item \code{spline_params}: A list of spline parameters used in the
#'   analysis, including:
#'     \itemize{
#'       \item \code{spline_type}: The type of spline (e.g., "n" for natural
#'       splines or "b" for B-splines).
#'       \item \code{dof}: Degrees of freedom for the spline.
#'       \item \code{knots}: Positions of the internal knots (for B-splines).
#'       \item \code{bknots}: Boundary knots (for B-splines).
#'       \item \code{degree}: Degree of the spline (for B-splines only).
#'     }
#' }
#'
#' @return The SplineOmics object, updated with a list with three elements:
#'         - `time_effect`: A list of top tables for each level with the time
#'                          effect.
#'         - `avrg_diff_conditions`: A list of top tables for each comparison
#'                                  between the levels. The comparison is the
#'                                  average difference of the values.
#'         - `interaction_condition_time`: A list of top tables for each
#'                                         comparison between levels. The
#'                                         comparison is the interaction between
#'                                         the condition and the time.
#'
#' @importFrom purrr partial map map_chr map2
#' @importFrom stats setNames
#' @importFrom utils combn
#'
#' @export
#'
run_limma_splines <- function(
    splineomics
    ) {
  
  check_splineomics_elements(
    splineomics = splineomics,
    func_type = "run_limma_splines"
  )

  args <- lapply(
    as.list(match.call()[-1]),
    eval,
    parent.frame()
  )

  check_null_elements(args)
  input_control <- InputControl$new(args)
  input_control$auto_validate()

  data <- splineomics[["data"]]
  rna_seq_data <- splineomics[["rna_seq_data"]]
  meta <- splineomics[["meta"]]
  spline_params <- splineomics[["spline_params"]]
  padjust_method <- splineomics[["padjust_method"]]
  design <- splineomics[["design"]]
  dream_params <- splineomics[["dream_params"]]
  mode <- splineomics[["mode"]]
  condition <- splineomics[["condition"]]

  feature_names <- rownames(data)

  rownames(data) <- NULL # To just have numbers describing the rows

  meta[[condition]] <- factor(meta[[condition]])
  levels <- levels(meta[[condition]])

  # Get hits for level (within level analysis)
  process_level_with_params <- purrr::partial(
    within_level,
    spline_params = spline_params,
    data = data,
    rna_seq_data = rna_seq_data,
    meta = meta,
    design = design,
    dream_params = dream_params,
    condition = condition,
    feature_names = feature_names,
    padjust_method = padjust_method,
    mode = mode
  )

  results_list <- purrr::imap(
    levels,
    process_level_with_params
  )

  within_level_top_table <-
    stats::setNames(
      purrr::map(results_list, "top_table"),
      purrr::map_chr(results_list, "name")
    )

  # Factor and Factor:Time comparisons between levels
  between_level_condition_only <- list()
  between_level_condition_time <- list() # Factor AND time

  if (mode == "integrated") {
    level_combinations <- utils::combn(levels, 2, simplify = FALSE)
    for (lev_combo in level_combinations) {
      result <- between_level(
        data = data,
        rna_seq_data = rna_seq_data,
        meta = meta,
        design = design,
        dream_params = dream_params,
        spline_params = spline_params,
        condition = condition,
        compared_levels = lev_combo,
        padjust_method = padjust_method,
        feature_names = feature_names
      )

      between_level_condition_only[[
        paste0(
          "avrg_diff_", lev_combo[1],
          "_vs_", lev_combo[2]
        )
      ]] <- result$condition_only

      between_level_condition_time[[
        paste0(
          "time_interaction_",
          lev_combo[1],
          "_vs_", lev_combo[2]
        )
      ]] <- result$condition_time
    }
  } else { # mode == "isolated"
    message(paste(
      "mode == 'integrated' necessary for between level",
      "comparisons. Returning emtpy lists the limma result categories 2 and 3
      (avrg diff conditions, and interaction condition time)."
    ))
  }

  message("\033[32mInfo\033[0m limma spline analysis completed successfully")

  limma_splines_result <- list(
    time_effect = within_level_top_table,
    avrg_diff_conditions = between_level_condition_only,
    interaction_condition_time = between_level_condition_time
  )

  splineomics <- update_splineomics(
    splineomics = splineomics,
    limma_splines_result = limma_splines_result
  )
}



# Level 1 internal functions ---------------------------------------------------


#' Between Level Analysis
#' 
#' @noRd
#'
#' @description
#' Performs a between-level analysis using LIMMA to compare specified levels
#' within a condition.
#'
#' @param data A matrix of data values.
#' @param rna_seq_data An object containing the preprocessed RNA-seq data,
#' such as the output from `limma::voom` or a similar preprocessing pipeline.
#' @param meta A dataframe containing metadata, including a 'Time' column.
#' @param design A design formula or matrix for the LIMMA analysis.
#' @param dream_params A named list or NULL. When not NULL, it must at least 
#' contain the named element 'random_effects', which must contain a string that
#' is a formula for the random effects of the mixed models by dream. 
#' Additionally, it can contain the named elements dof, which must be a int
#' bigger than 1, which is the degree of freedom for the dream topTable, and
#' the named element KenwardRoger, which must be a bool, specifying whether
#' to use that method or not.
#' @param spline_params A list of spline parameters for the analysis.
#' @param condition A character string specifying the condition.
#' @param compared_levels A vector of levels within the condition to compare.
#' @param padjust_method A character string specifying the p-adjustment method.
#' @param feature_names A non-empty character vector of feature names.
#'
#' @return A list containing top tables for the factor only and factor-time
#' contrast.
#'
#' @seealso
#' \code{\link[splines]{bs}}, \code{\link[splines]{ns}},
#' \code{\link[limma]{lmFit}}, \code{\link[limma]{eBayes}},
#' \code{\link[limma]{topTable}}, \code{\link{modify_limma_top_table}}
#'
#' @importFrom splines bs ns
#' @importFrom stats as.formula model.matrix
#' @importFrom limma lmFit eBayes topTable
#' @importFrom variancePartition dream eBayes topTable
#'
between_level <- function(
    data,
    rna_seq_data,
    meta,
    design,
    dream_params,
    spline_params,
    condition,
    compared_levels,
    padjust_method,
    feature_names
    ) {

  samples <- which(meta[[condition]] %in% compared_levels)
  data <- data[, samples]

  meta <- meta[meta[[condition]] %in% compared_levels, ]

  result <- design2design_matrix(
    meta = meta,
    spline_params = spline_params,
    level_index = 1,
    design = design
  )
  
  design_matrix <- result$design_matrix

  if (!is.null(rna_seq_data)) {
    data <- rna_seq_data
  }
  
  condition_only_contrast_coeff <- paste0(
    condition,
    compared_levels[2]
  )
  
  num_matching_columns <- sum(
    grepl(
      "^X\\d+$",
      colnames(design_matrix)
    )
  )
  
  interaction_condition_time_contrast_coeffs <- paste0(
    condition,
    compared_levels[2],
    ":X",
    seq_len(num_matching_columns)
  )

  if (!is.null(dream_params)) {
    colnames(data) <- rownames(meta)  # dream requires this format
    meta <- result$meta  # Only do it after, new meta has more columns
    
    # Apply the Kenward-Roger method if specified
    if (isTRUE(dream_params[["KenwardRoger"]])) {
      method <- "Kenward-Roger"
    } else {
      method <- NULL
    }

    fit <- variancePartition::dream(
      exprObj = data,
      formula = stats::as.formula(design),
      data = meta,
      random.formula = stats::as.formula(dream_params[["random_effects"]]),
      ddf = method
    )
    
    fit <- variancePartition::eBayes(fit)

    condition_only <- variancePartition::topTable(
      fit,
      coef = condition_only_contrast_coeff,
      adjust.method = padjust_method,
      number = Inf
    )
    
    condition_time <- variancePartition::topTable(
      fit,
      coef = interaction_condition_time_contrast_coeffs,
      adjust.method = padjust_method,
      number = Inf,
      sort.by = "F"
    )
    
  } else {
    fit <- limma::lmFit(
      data,
      design_matrix
    )
    
    fit <- limma::eBayes(fit)

    condition_only <- limma::topTable(
      fit,
      coef = condition_only_contrast_coeff,
      adjust.method = padjust_method,
      number = Inf
    )
  
    condition_time <- limma::topTable(
      fit,
      coef = interaction_condition_time_contrast_coeffs,
      adjust.method = padjust_method,
      number = Inf
    )
  }
  
  condition_only_resuls <- list(
    top_table = condition_only,
    fit = fit
  )
  
  top_table_condition_only <- process_top_table(
    condition_only_resuls,
    feature_names
  )
  
  condition_and_time_results <- list(
    top_table = condition_time,
    fit = fit
  )
  top_table_condition_and_time <- process_top_table(
    condition_and_time_results,
    feature_names
  )

  list(
    condition_only = top_table_condition_only,
    condition_time = top_table_condition_and_time
  )
}


#' Within level analysis
#' 
#' @noRd
#'
#' @description
#' Processes a single level within a condition, performing limma analysis
#' and generating the top table of results.
#'
#' @param level The level within the condition to process.
#' @param level_index The index of the level within the condition.
#' @param spline_params A list of spline parameters for the analysis.
#' @param data A matrix of data values.
#' @param rna_seq_data An object containing the preprocessed RNA-seq data,
#' such as the output from `limma::voom` or a similar preprocessing pipeline.
#' @param meta A dataframe containing the metadata for data.
#' @param design A design formula or matrix for the limma analysis.
#' @param dream_params A named list or NULL. When not NULL, it must at least 
#' contain the named element 'random_effects', which must contain a string that
#' is a formula for the random effects of the mixed models by dream. 
#' Additionally, it can contain the named elements dof, which must be a int
#' bigger than 1, which is the degree of freedom for the dream topTable, and
#' the named element KenwardRoger, which must be a bool, specifying whether
#' to use that method or not.
#' @param condition A character string specifying the condition.
#' @param feature_names A non-empty character vector of feature names.
#' @param padjust_method A character string specifying the p-adjustment method.
#' @param mode A character string specifying the mode
#'            ('isolated' or 'integrated').
#'
#' @return A list containing the name of the results and the top table of
#'          results.
#'
#' @seealso
#' \code{\link{within_level}}, \code{\link{process_top_table}}
#'
#' @importFrom stats relevel
#'
within_level <- function(
    level,
    level_index,
    spline_params,
    data,
    rna_seq_data,
    meta,
    design,
    dream_params,
    condition,
    feature_names,
    padjust_method,
    mode) {
  if (mode == "isolated") {
    samples <- which(meta[[condition]] == level)
    data_copy <- data[, samples]
    meta_copy <- meta[meta[[condition]] == level, , drop = FALSE]
  } else { # mode == "integrated"
    data_copy <- data
    meta_copy <- meta
    meta_copy[[condition]] <- stats::relevel(
      meta_copy[[condition]],
      ref = level
    )
    # spline_params must be uniform across all levels for integrated mode.
    level_index <- 1L
  }

  result <- process_within_level(
    data = data_copy,
    rna_seq_data = rna_seq_data,
    meta = meta_copy,
    design = design,
    dream_params = dream_params,
    spline_params = spline_params,
    level_index = level_index,
    padjust_method = padjust_method
  )

  top_table <- process_top_table(
    result,
    feature_names
  )

  results_name <- paste(
    condition,
    level,
    sep = "_"
  )

  list(
    name = results_name,
    top_table = top_table
  )
}


# Level 2 internal functions ---------------------------------------------------


#' Process Top Table
#' 
#' @noRd
#'
#' @description
#' Processes the top table from a LIMMA analysis, adding feature names and
#' intercepts.
#'
#' @param process_within_level_result List of lists containing the limma
#'                                    topTable, and fit. All of this is from
#'                                    one specific level.
#' @param feature_names A non-empty character vector of feature names.
#'
#' @return A dataframe containing the processed top table with added intercepts.
#'
#' @seealso
#' \link{modify_limma_top_table}, \link[limma]{lmFit}
#'
#' @importFrom stats coef
#'
process_top_table <- function(
    process_within_level_result,
    feature_names) {
  top_table <- process_within_level_result$top_table
  fit <- process_within_level_result$fit

  top_table <- modify_limma_top_table(
    top_table,
    feature_names
  )

  intercepts <- as.data.frame(stats::coef(fit)[, "(Intercept)", drop = FALSE])
  intercepts_ordered <- intercepts[top_table$feature_nr, , drop = FALSE]
  top_table$intercept <- intercepts_ordered[, 1]

  top_table
}


#' Process Within Level
#' 
#' @noRd
#'
#' @description
#' Performs a within-level analysis using limma to generate top tables and fit
#' objects based on the specified spline parameters. Performs the limma spline
#' analysis for a selected level of a factor
#'
#' @param data A matrix of data values.
#' @param rna_seq_data An object containing the preprocessed RNA-seq data,
#' such as the output from `limma::voom` or a similar preprocessing pipeline.
#' @param meta A dataframe containing metadata, including a 'Time' column.
#' @param design A design formula or matrix for the limma analysis.
#' @param dream_params A named list or NULL. When not NULL, it must at least 
#' contain the named element 'random_effects', which must contain a string that
#' is a formula for the random effects of the mixed models by dream. 
#' Additionally, it can contain the named elements dof, which must be a int
#' bigger than 1, which is the degree of freedom for the dream topTable, and
#' the named element KenwardRoger, which must be a bool, specifying whether
#' to use that method or not.
#' @param spline_params A list of spline parameters for the analysis.
#' @param level_index The index of the level within the factor.
#' @param padjust_method A character string specifying the p-adjustment method.
#'
#' @return A list containing the top table and the fit object from the limma
#' analysis.
#'
#' @seealso
#' \link[splines]{bs}, \link[splines]{ns}, \link[limma]{lmFit},
#' \link[limma]{eBayes}, \link[limma]{topTable}
#'
#' @importFrom splines bs ns
#' @importFrom stats as.formula model.matrix
#' @importFrom limma lmFit eBayes topTable
#' @importFrom variancePartition dream eBayes topTable
#'
process_within_level <- function(
    data,
    rna_seq_data,
    meta,
    design,
    dream_params,
    spline_params,
    level_index,
    padjust_method
    ) {
  
  result <- design2design_matrix(
    meta,
    spline_params,
    level_index,
    design
  )
  
  design_matrix <- result$design_matrix

  if (!is.null(rna_seq_data)) {
    data <- rna_seq_data
  }

  if (!is.null(dream_params)) {
    colnames(data) <- rownames(meta)  # dream wants it like this.
    meta <- result$meta  # Only do it after, new meta has more columns

    if (isTRUE(dream_params[["KenwardRoger"]])) {
      method <- "Kenward-Roger"
    } else {
      method = NULL
    }
    
    fit <- variancePartition::dream(
      exprObj = data,
      formula = stats::as.formula(design),
      data = meta,
      random.formula = stats::as.formula(dream_params[["random_effects"]]),
      ddf = method
    )
    
    fit <- variancePartition::eBayes(fit)
    
    num_matching_columns <- sum(grepl("^X\\d+$", colnames(design_matrix)))
    coeffs <- paste0("X", seq_len(num_matching_columns))
    
    if (!is.null(dream_params[["dof"]])) {
      dof <- dream_params[["dof"]]
    } else {
      dof = Inf
    }

    top_table <- variancePartition::topTable(
      fit,
      adjust.method = padjust_method,
      number = dof,
      coef = coeffs,
      sort.by = "F"
    )
  } else {
    fit <- limma::lmFit(
      data,
      design_matrix
    )
    fit <- limma::eBayes(fit)
    
    # Extract the top table based on coefficients
    num_matching_columns <- sum(grepl("^X\\d+$", colnames(design_matrix)))
    coeffs <- paste0("X", seq_len(num_matching_columns))
    
    top_table <- limma::topTable(
      fit,
      adjust.method = padjust_method,
      number = Inf,
      coef = coeffs
    )
  }

  attr(top_table, "adjust.method") <- padjust_method

  list(
    top_table = top_table,
    fit = fit
  )
}


#' Remove intercept from a formula
#' 
#' @noRd
#'
#' @description
#' This function modifies a given formula by replacing the first occurrence
#' of a standalone intercept (`1`) with `0`. It works even if the `1` is
#' preceded by a tilde (`~`), ensuring that the intercept is removed while
#' leaving other parts of the formula intact.
#'
#' @param formula A formula object. The formula can include an intercept (`1`)
#'   and other terms. If a `1` is found, it is replaced with `0`.
#'
#' @return A modified formula with the intercept removed. The first standalone
#'   occurrence of `1` will be replaced by `0`.
#'
remove_intercept <- function(formula) {
  # Deparse the formula into a single string
  formula_str <- paste(deparse(formula), collapse = " ")
  
  # Regular expression to match the first standalone 1
  # (with optional spaces around it)
  pattern <- "(~\\s*|\\s+)(1)(\\s|$)"
  
  # Replace the first occurrence of "1" with "0"
  formula_str <- sub(pattern, "\\10\\3", formula_str)
  
  new_formula <- as.formula(formula_str)
}


# Level 3 internal functions ---------------------------------------------------


#' Modify limma Top Table
#' 
#' @noRd
#'
#' @description
#' Modifies the limma top table to include feature indices and names.
#'
#' @param top_table A dataframe containing the top table results from limma
#' @param feature_names A character vector of feature names.
#'
#' @return A tibble with feature indices and names included.
#'
#' @importFrom tibble as_tibble
#' @importFrom dplyr relocate last_col mutate
#' @importFrom rlang .data
#'
modify_limma_top_table <- function(
    top_table,
    feature_names) {
  is_integer_string <- function(x) {
    return(grepl("^[0-9]+$", x))
  }

  # Because the row headers of a potential rna_seq_data object were not
  # converted to ints (written as strings) beforehand. This is run only when
  # the row headers are still "real" strings.
  if (!all(vapply(
    rownames(top_table),
    is_integer_string,
    logical(1)
  ))
  ) {
    rownames(top_table) <- vapply(
      rownames(top_table),
      function(id) {
        # Find the index of the current row name in feature_names
        index <- which(feature_names == id)
        # Return the index as a string
        return(as.character(index))
      },
      character(1)
    )
  }

  top_table <- tibble::as_tibble(
    top_table,
    rownames = "feature_nr"
  )

  # feature_nr <- NULL  # dummy declaration for the lintr and R CMD.

  # Convert feature_nr to integer
  top_table <- top_table |>
    dplyr::mutate(feature_nr = as.integer(.data$feature_nr)) |>
    dplyr::relocate(.data$feature_nr, .after = dplyr::last_col())

  # Sort and add feature names based on the feature_nr
  sorted_feature_names <- feature_names[top_table$feature_nr]
  top_table <- top_table |> dplyr::mutate(feature_names = sorted_feature_names)

  return(top_table)
}
