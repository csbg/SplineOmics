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
#'   \item\code{dream_params}: #' A named list or NULL. When not NULL, it can 
#'                                contain the following named 
#'      elements:
#'      - `dof`: An integer greater than 1, specifying the degrees of freedom 
#'      for  the dream topTable.
#'      - `KenwardRoger`: A boolean indicating whether to use the Kenward-Roger 
#'      approximation for mixed models.
#'      Note that random effects are now directly specified in the design 
#'      formula and not in `dream_params`.
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
  use_array_weights <- splineomics[["use_array_weights"]]
  
  # Because at first I enforced that X in the design formula stands for the time
  # and I heavily oriented my code towards that. But then I realised that it is
  # nonsense to encode the time as X, and now it is explicitly "Time" (because
  # meta must contain the exact name "Time" for this respective column).
  design <- gsub("Time", "X", design)  
  
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
  
  limma_splines_result <- list(
    time_effect = within_level_top_table
  )
  
  if (mode == "integrated") {
    # Step 1: Fit the global model once
    fit_obj <- fit_global_model(
      data = data,
      rna_seq_data = rna_seq_data,
      meta = meta,
      design = design,
      dream_params = dream_params,
      spline_params = spline_params,
      condition = condition,
      padjust_method = padjust_method,
      feature_names = feature_names,
      use_array_weights = use_array_weights
    )
    
    # Step 2: Extract pairwise contrasts for all level combinations
    contrast_results <- extract_between_level_contrasts(
      fit_obj = fit_obj,
      condition = condition
    )
    
    # Add contrast-based results
    limma_splines_result[["avrg_diff_conditions"]] <- 
      contrast_results[["condition_only"]]
    limma_splines_result[["interaction_condition_time"]] <- 
      contrast_results[["condition_time"]]
    
  } else { # mode == "isolated"
    message(paste(
      "mode == 'integrated' necessary for between level",
      "comparisons. Returning emtpy lists the limma result categories 2 and 3",
      "(avrg diff conditions, and interaction condition time)."
    ))
  }
  
  message("\033[32mInfo\033[0m limma spline analysis completed successfully")
  
  splineomics <- update_splineomics(
    splineomics = splineomics,
    limma_splines_result = limma_splines_result,
    homosc_violation_result = 
      if (exists("fit_obj")) fit_obj[["homosc_violation_result"]] else NULL
  )
}



# Level 1 internal functions ---------------------------------------------------


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
    mode
) {
  
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
  
  if (!is.null(rna_seq_data) && ncol(rna_seq_data$E) != nrow(meta_copy)) {
    stop_call_false(
      "Mismatch detected: rna_seq_data$E has ", ncol(rna_seq_data$E),
      " columns, but meta_copy has ", nrow(meta_copy), " rows. The most likely",
      "cause for this is that you selected mode == isolated, but passed the ",
      "full data in rna_seq_data. For RNA-seq data, you must pass the ",
      "respective datasets of the different conditions individually (for all ",
      "other omics datasets, this is handleded implicitly by SplineOmics: it ",
      "splits up the data and meta into the different conditions. However, ",
      "that is not possible with the RNA-seq data objects"
    )
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


#' Between Level Analysis
#' 
#' @noRd
#'
#' @description
#' Performs a between-level analysis using limma to compare specified levels
#' within a condition.
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
#' @param condition A character string of the column name of meta that contains
#'                  the levels of the experimental condition.
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
#' @importFrom variancePartition dream eBayes
#'
fit_global_model <- function(
    data,
    rna_seq_data,
    meta,
    design,
    dream_params,
    spline_params,
    condition,
    padjust_method,
    feature_names,
    use_array_weights = NULL
) {
  
  effects <- extract_effects(design)
  
  design2design_matrix_result <- design2design_matrix(
    meta = meta,
    spline_params = spline_params,
    level_index = 1,
    design = effects[["fixed_effects"]]
  )
  
  design_matrix <- design2design_matrix_result[["design_matrix"]]
  
  if (!is.null(rna_seq_data)) {
    data <- rna_seq_data   # Just having one variable makes the code easier
  } 
  
  # For RNA-seq data, this is handled when calling limma::voom (happens before)
  # This here is the implicit fall-back logic when the user has not explicitly 
  # decided whether to use the array_weights strategy or not (is NULL)
  if (is.null(rna_seq_data) && is.null(use_array_weights)) {
    homosc_violation_result <- check_homoscedasticity_violation(
      data = data,
      meta = meta,
      design = design,
      design2design_matrix_result = design2design_matrix_result,
      condition = condition,
      random_effects = effects[["random_effects"]] != ""  # Boolean flag
    )
    
    # If there is a considerable violation, select use_array_weights strategy
    use_array_weights <- homosc_violation_result[["violation"]]
  } else{
    homosc_violation_result <- NULL
  }
  
  if (effects[["random_effects"]] != "") {
    colnames(data) <- rownames(meta)  # dream requires this format
    
    # Apply the Kenward-Roger method if specified
    if (isTRUE(dream_params[["KenwardRoger"]])) {
      method <- "Kenward-Roger"
    } else {
      method <- NULL
    }
    
    if (use_array_weights) {
      aw <- limma::arrayWeights(      # vector of length = # samples
        object = data,
        design = design_matrix
      ) 
      weights_matrix <- matrix(
        rep(aw, each = nrow(data)),
        nrow = nrow(data),
        byrow = TRUE
      )
      fit <- variancePartition::dream(
        exprObj = data,
        formula = stats::as.formula(design),
        data = design2design_matrix_result[["meta"]], # Spline transformed meta.
        ddf = method,
        useWeights = TRUE,
        weightsMatrix = weights_matrix
      )
      fit <- variancePartition::eBayes(
        fit = fit,
        robust = TRUE
      ) 
    } else {
      fit <- variancePartition::dream(
        exprObj = data,
        formula = stats::as.formula(design),
        data = design2design_matrix_result[["meta"]], # Spline transformed meta.
        ddf = method
      )
      
      fit <- variancePartition::eBayes(fit = fit) 
    }
  } else {
    if (use_array_weights) {
      weights <- limma::arrayWeights(
        object = data,
        design = design_matrix
      )
      fit <- limma::lmFit(
        object = data,
        design = design_matrix,
        weights = weights
      )
      fit <- limma::eBayes(
        fit = fit,
        robust = TRUE
      )
    } else {
      fit <- limma::lmFit(
        object = data,
        design = design_matrix
      )
      fit <- limma::eBayes(fit = fit)
    }
  }
  
  list(
    fit = fit,
    design_matrix = design_matrix,
    meta = design2design_matrix_result[["meta"]],
    feature_names = feature_names,
    condition = condition,
    padjust_method = padjust_method,
    homosc_violation_result = homosc_violation_result
  )
}


#' Extract pairwise contrasts for a condition from a fitted limma model
#'
#' @noRd
#'
#' @description
#' Internal helper function that extracts condition-only and interaction
#' contrasts for all pairwise combinations of levels within a condition
#' factor. Works with fitted limma or dream models stored in a structured
#' list returned by a global model fitting function.
#'
#' @param fit_obj A list containing:
#'   - `fit`: fitted model object from `lmFit()` or `dream()`
#'   - `meta`: sample metadata used in the model
#'   - `design_matrix`: the design matrix used for fitting
#'   - `feature_names`: optional row annotations
#'   - `condition`: the condition factor used
#'   - `padjust_method`: method for p-value adjustment
#'
#' @param condition A string indicating the column in `meta` that contains
#'   the condition factor for which pairwise contrasts should be extracted.
#'
#' @return A named list with two elements:
#'   - `condition_only`: average condition-level differences
#'   - `condition_time`: condition × time interaction effects
#'
#' Each is itself a list of results for all level combinations.
#'
extract_between_level_contrasts <- function(
    fit_obj,
    condition
) {
  
  levels <- unique(fit_obj$meta[[condition]])
  level_combinations <- utils::combn(
    levels,
    2,
    simplify = FALSE
  )
  
  between_level_condition_only <- list()
  between_level_condition_time <- list()
  
  for (lev_combo in level_combinations) {
    contrast_result <- extract_contrast_for_pair(
      fit_obj = fit_obj,
      condition = condition,
      level_pair = lev_combo
    )
    
    between_level_condition_only[[
      paste0(
        "avrg_diff_",
        lev_combo[1],
        "_vs_",
        lev_combo[2]
      )
    ]] <- contrast_result$condition_only
    
    between_level_condition_time[[
      paste0(
        "time_interaction_",
        lev_combo[1],
        "_vs_",
        lev_combo[2]
      )
    ]] <- contrast_result$condition_time
  }
  
  list(
    condition_only = between_level_condition_only,
    condition_time = between_level_condition_time
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
#' @param dream_params A named list or NULL. When not NULL, it it can contain 
#' the named elements dof, which must be a int
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
#' @importFrom variancePartition dream eBayes
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
  
  effects <- extract_effects(design)
  
  result <- design2design_matrix(
    meta,
    spline_params,
    level_index,
    design = effects[["fixed_effects"]]
  )
  
  design_matrix <- result[["design_matrix"]]
  
  if (!is.null(rna_seq_data)) {
    data <- rna_seq_data
  }
  
  if (effects[["random_effects"]] != "") {
    colnames(data) <- rownames(meta)  # dream wants it like this.
    
    if (isTRUE(dream_params[["KenwardRoger"]])) {
      method <- "Kenward-Roger"
    } else {
      method <- "adaptive"  # Kenward-Roger for < 20 samples, else Satterthwaite
    }
    
    fit <- variancePartition::dream(
      exprObj = data,
      formula = stats::as.formula(design),
      data = result[["meta"]],
      ddf = method
    )
    
    fit <- variancePartition::eBayes(fit)
    
    num_matching_columns <- sum(grepl("^X\\d+$", colnames(design_matrix)))
    coeffs <- paste0("X", seq_len(num_matching_columns))
    
    if (!is.null(dream_params[["dof"]])) {
      dof <- dream_params[["dof"]]
    } else {
      dof <- Inf
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


#' Extract contrasts for a single pair of condition levels
#'
#' @noRd
#'
#' @description
#' Extracts both the condition-only and condition-time interaction
#' contrasts for a given pair of condition levels from a fitted limma
#' model. Designed to be called from a wrapper function that loops over
#' all pairwise combinations.
#'
#' @param fit_obj A list containing:
#'   - `fit`: fitted model object from `lmFit()` or `dream()`
#'   - `design_matrix`: the design matrix used in fitting
#'   - `feature_names`: optional row annotations
#'   - `padjust_method`: method for p-value adjustment
#'
#' @param condition A string giving the name of the condition factor used
#'   in the design.
#'
#' @param level_pair A character vector of length 2 giving the levels to
#'   compare.
#'
#' @return A named list with two elements:
#'   - `condition_only`: the top table for the condition effect
#'   - `condition_time`: the top table for the interaction effect
#'   
extract_contrast_for_pair <- function(
    fit_obj,
    condition,
    level_pair
) {
  
  fit <- fit_obj[["fit"]]
  design_matrix <- fit_obj[["design_matrix"]]
  feature_names <- fit_obj[["feature_names"]]
  padjust_method <- fit_obj[["padjust_method"]]
  
  # Get properly constructed contrasts
  contrasts <- build_contrasts_for_pair(
    condition = condition,
    level_pair = level_pair,
    design_matrix = design_matrix
  )
  
  # Extract condition-only top table
  condition_only <- limma::topTable(
    fit,
    coef = contrasts$condition,
    adjust.method = padjust_method,
    number = Inf
  )
  
  # Extract condition-time interaction top table
  condition_time <- limma::topTable(
    fit,
    coef = contrasts$interaction,
    adjust.method = padjust_method,
    number = Inf,
    sort.by = "F"
  )
  
  # Wrap results in lists to pass through process_top_table
  condition_only_result <- list(
    top_table = condition_only,
    fit = fit
  )
  condition_time_result <- list(
    top_table = condition_time,
    fit = fit
  )
  
  list(
    condition_only = process_top_table(
      condition_only_result,
      feature_names
    ),
    condition_time = process_top_table(
      condition_time_result,
      feature_names
    )
  )
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
    feature_names
    ) {
  
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

  # Convert feature_nr to integer
  top_table <- top_table |>
    dplyr::mutate(feature_nr = as.integer(.data$feature_nr)) |>
    dplyr::relocate(feature_nr, .after = dplyr::last_col())
  
  # Sort and add feature names based on the feature_nr
  sorted_feature_names <- feature_names[top_table$feature_nr]
  top_table <- top_table |> dplyr::mutate(feature_names = sorted_feature_names)
  
  return(top_table)
}


#' Construct contrast names for a pair of condition levels
#'
#' @noRd
#'
#' @description
#' Builds the coefficient names used for extracting the condition-only and
#' condition-time interaction effects between two condition levels. Assumes
#' a design matrix where one condition level is used as the reference.
#'
#' @param condition A string giving the name of the condition factor.
#'
#' @param level_pair A character vector of length 2 specifying the two levels
#'   to be compared.
#'
#' @param design_matrix The design matrix used for model fitting, whose column
#'   names determine which condition level is modeled explicitly.
#'
#' @return A named list with two elements:
#'   - `condition`: name of the coefficient for the condition-only effect
#'   - `interaction`: vector of coefficient names for the interaction terms
#'   
build_contrasts_for_pair <- function(
    condition,
    level_pair,
    design_matrix
) {
  
  cols <- colnames(design_matrix)
  
  # Try each level
  level1 <- make.names(paste0(
    condition,
    level_pair[1]
  ))
  level2 <- make.names(paste0(
    condition, level_pair[2]
  ))
  
  if (any(grepl(level1, cols))) {
    modeled_level <- level_pair[1]
  } else if (any(grepl(level2, cols))) {
    modeled_level <- level_pair[2]
  } else {
    stop("Neither level appears in design_matrix columns.")
  }
  
  # Build contrasts
  condition_contrast <- make.names(paste0(
    condition,
    modeled_level
  ))
  
  # Time interaction columns like X1, X2
  spline_cols <- grep(
    "^X\\d+$",
    cols,
    value = TRUE
  )
  interaction_contrasts <- paste0(
    condition_contrast,
    ":",
    spline_cols
  )
  
  list(
    condition = condition_contrast,
    interaction = interaction_contrasts
  )
}