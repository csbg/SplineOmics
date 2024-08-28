#' run_limma_splines.R contains the exported package function run_limma_splines  
#' and all the functions that make up the functionality of run_limma_splines.
#' run_limma_splines performs a limma analysis, using splines, to assign a 
#' p-value to every feature of a time series omics dataset, to find out which
#' features are significantly changed over the time course.


# Exported function: run_limma_splines() ---------------------------------------

 
#' Run Limma Analysis with Spline Interpolation for Hyperparameter Screening
#'
#' This function conducts differential expression analysis using the Limma 
#' package, 
#' incorporating spline interpolation to model the effect of various 
#' experimental 
#' factors across different levels. It supports both isolated and integrated 
#' modes 
#' for within-level analysis and between-level comparison, adjusting for 
#' multiple 
#' degrees of freedom corresponding to the factors under investigation.
#'
#' @param splineomics A SplineOmics object, containing data, meta, design, 
#'                    condition, and spline_params.
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
  
  design <- splineomics[["design"]]
  condition <- splineomics[["condition"]]
  mode <- determine_analysis_mode(
    design,
    condition
    )
  
  args <- lapply(
    as.list(match.call()[-1]),
    eval,
    parent.frame()
    )
  
  args$mode <- mode
  check_null_elements(args)
  input_control <- InputControl$new(args)
  input_control$auto_validate()
  
  data <- splineomics[["data"]]
  preprocess_rna_seq <- splineomics[["preprocess_rna_seq"]]
  normalization_fun <- splineomics[["normalization_fun"]]
  meta <- splineomics[["meta"]]
  spline_params <- splineomics[["spline_params"]]
  padjust_method <- splineomics[["padjust_method"]]
  
  feature_names <- rownames(data)
  data_copy <- data
  rownames(data_copy) <- NULL   # To just have numbers describing the rows

  meta[[condition]] <- factor(meta[[condition]])
  levels <- levels(meta[[condition]])
  
  # Get hits for level (within level analysis) 
  process_level_with_params <- purrr::partial(
    within_level, 
    spline_params = spline_params,
    data = data_copy, 
    preprocess_rna_seq = preprocess_rna_seq,
    normalization_fun = normalization_fun,
    meta = meta, 
    design = design, 
    condition = condition, 
    feature_names = feature_names, 
    padjust_method = padjust_method, 
    mode = mode
    )

  results_list <- purrr::map2(
    levels, 
    seq_along(levels), 
    process_level_with_params
    )

  within_level_top_table <- 
    stats::setNames(
      purrr::map(results_list, "top_table"), 
      purrr::map_chr(results_list, "name")
      )

  # For RNA-seq data, voom$E data matrices must be passed to cluster_hits()
  voom_matrices <- lapply(
    results_list,
    function(x) x$voom_data_matrix_level
    )

  if (!any(sapply(voom_matrices, is.null))) {
    if (args$mode == "isolated") {
      data <- do.call(rbind, voom_matrices)   # Combine from all levels
    } else {    # mode == "integrated"
      # All levels contain the full data. Can just take the first one.
      data <- voom_matrices[[1]]  
    }
    rownames(data) <- feature_names  # Readd the original row headers.
  }
  
  # Factor and Factor:Time comparisons between levels 
  between_level_condition_only <- list()        
  between_level_condition_time <- list()        # Factor AND time
  
  if (mode == "integrated") {
    level_combinations <- utils::combn(levels, 2, simplify = FALSE)
    for (lev_combo in level_combinations) {
      result <- between_level(
        data = data_copy, 
        preprocess_rna_seq = preprocess_rna_seq,
        normalization_fun = normalization_fun,
        meta = meta, 
        design = design, 
        spline_params = spline_params,
        condition = condition, 
        compared_levels = lev_combo, 
        padjust_method = padjust_method, 
        feature_names = feature_names
        )
      
      between_level_condition_only[[
        paste0(
          "avrg_diff_" ,lev_combo[1],
          "_vs_", lev_combo[2]
          )
        ]] <- result$condition_only
      
      between_level_condition_time[[
        paste0(
          "time_interaction_" ,
          lev_combo[1],
          "_vs_", lev_combo[2]
          )
        ]] <- result$condition_time
    }
  } else { # mode == "isolated"
    message(paste(
      "mode == 'integrated' necessary for between level",
      "comparisons. Returning emtpy lists for ttslc_factor_only",
      "and ttslc_factor_time (ttslc means 'top tables level comparison')."
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
    data = data,   # In case voom_data_matrix has been generated.
    limma_splines_result = limma_splines_result
 )
}



# Level 1 internal functions ---------------------------------------------------


#' Between Level Analysis
#'
#' @description
#' Performs a between-level analysis using LIMMA to compare specified levels 
#' within a condition.
#'
#' @param data A matrix of data values.
#' @param preprocess_rna_seq Boolean specifying whether to preprocess RNA seq
#' @param normalization_fun Function for normalizing RNA-seq raw-counts.
#' @param meta A dataframe containing metadata, including a 'Time' column.
#' @param design A design formula or matrix for the LIMMA analysis.
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
#' @importFrom splines bs
#' @importFrom splines ns
#' @importFrom stats as.formula
#' @importFrom stats model.matrix
#' @importFrom limma lmFit
#' @importFrom limma eBayes
#' @importFrom limma topTable
#' 
between_level <- function(
    data, 
    preprocess_rna_seq,
    normalization_fun,
    meta, 
    design, 
    spline_params, 
    condition, 
    compared_levels,
    padjust_method, 
    feature_names
    ) {
  
  samples <- which(meta[[condition]] %in% compared_levels)
  data <- data[, samples]
  
  meta <- subset(meta, meta[[condition]] %in% compared_levels)
  
  design_matrix <- design2design_matrix(
    meta = meta,
    spline_params = spline_params,
    level_index = 1,
    design = design
    )
  
  if (preprocess_rna_seq) {
    data <- preprocess_rna_seq_data(
      raw_counts = data,
      design_matrix = design_matrix,
      normalization_fun
    )
  }
  
  fit <- limma::lmFit(
    data,
    design_matrix
    )
  fit <- limma::eBayes(fit)

  factor_only_contrast_coeff <- paste0(
    condition,
    compared_levels[2]
    )
  condition_only <- limma::topTable(
    fit,
    coef = factor_only_contrast_coeff,
    adjust.method = padjust_method,
    number = Inf
    )
  
  condition_only_resuls <- list(
    top_table = condition_only,
    fit = fit
    )
  top_table_condition_only <- process_top_table(
    condition_only_resuls, 
    feature_names
    )
  
  
  num_matching_columns <- sum(
    grepl(
      "^X\\d+$",
      colnames(design_matrix)
      )
    )
  
  factor_time_contrast_coeffs <- paste0(
    condition,
    compared_levels[2],
    ":X", 
    seq_len(num_matching_columns)
    )
  
  condition_time <- limma::topTable(
    fit,
    coef = factor_time_contrast_coeffs,
    adjust.method = padjust_method,
    number = Inf
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
#' @description
#' Processes a single level within a condition, performing LIMMA analysis 
#' and generating the top table of results.
#'
#' @param level The level within the condition to process.
#' @param level_index The index of the level within the condition.
#' @param spline_params A list of spline parameters for the analysis.
#' @param data A matrix of data values.
#' @param preprocess_rna_seq Boolean specifying whether to preprocess RNA seq
#' @param normalization_fun Function to normalize RNA-seq raw counts.
#' @param meta A dataframe containing metadata.
#' @param design A design formula or matrix for the LIMMA analysis.
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
    preprocess_rna_seq,
    normalization_fun,
    meta, 
    design, 
    condition, 
    feature_names, 
    padjust_method, 
    mode
    ) {
  
  if (mode == "isolated") {
    samples <- which(meta[[condition]] == level)
    data_copy <- data[, samples]
    meta_copy <- subset(meta, meta[[condition]] == level)
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
    preprocess_rna_seq = preprocess_rna_seq,
    normalization_fun = normalization_fun,
    meta = meta_copy, 
    design = design, 
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
    top_table = top_table,
    voom_data_matrix_level = result$voom_data_matrix
    )
}


# Level 2 internal functions ---------------------------------------------------


#' Perform default preprocessing of raw RNA-seq counts
#'
#' @description
#' This function is called when `preprocess_rna_seq` is `TRUE`. It performs the 
#' default preprocessing steps for raw RNA-seq counts, including creating a 
#' `DGEList` object, normalizing the counts, and applying the `voom` 
#' transformation.
#'
#' @param raw_counts A matrix of raw RNA-seq counts (genes as rows, samples as
#'  columns).
#' @param design_matrix A design matrix used in the linear modeling, typically
#'  specifying the experimental conditions.
#' @param normalize_func An optional normalization function. If provided, this 
#' function will be used to normalize the `DGEList` object. If not provided,
#'  TMM normalization (via `edgeR::calcNormFactors`) will be used by default.
#'
#' @return A `voom` object, which includes the log2-counts per million (logCPM)
#'  matrix and observation-specific weights.
#' 
#' @importFrom limma voom
#'   
preprocess_rna_seq_data <- function(
    raw_counts,
    design_matrix,
    normalize_func = NULL
) {
  
  message("Preprocessing RNA-seq data (normalization + voom)...")
  
  # Check if edgeR is installed; if not, prompt the user
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    message("The 'edgeR' package is not installed.")
    
    # Prompt user for action
    repeat {
      user_input <- readline(prompt = 
                               "What would you like to do?\n
                               1: Automatically install edgeR\n
                               2: Manually install edgeR\n
                               3: Cancel\n
                               Please enter 1, 2, or 3: "
                             )
      
      if (user_input == "1") {
        # Try to install edgeR automatically from Bioconductor
        message("Attempting to install 'edgeR' automatically 
                from Bioconductor...")
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager")
        }
        tryCatch(
          {
            BiocManager::install("edgeR", update = FALSE)
          },
          error = function(e) {
            stop(
            "Automatic installation of 'edgeR' failed.
            Please install it manually and try again.",
            call. = FALSE
            )
          }
        )
        break  # Exit the loop if installation is successful
      } else if (user_input == "2") {
        stop(
        "Please install 'edgeR' manually using 
         BiocManager::install('edgeR') and then re-run the function.",
         call. = FALSE
        )
      } else if (user_input == "3") {
        stop("Operation canceled by the user.", call. = FALSE)
      } else {
        message("Invalid input. Please enter 1, 2, or 3.")
      }
    }
  }
  
  # Load edgeR functions after installation check
  require(edgeR)
  
  # Step 1: Create DGEList object from raw counts
  y <- edgeR::DGEList(counts = raw_counts)
  
  # Step 2: Apply the normalization function (either user-provided or default)
  if (!is.null(normalize_func) && is.function(normalize_func)) {
    y <- normalize_func(y)   # user provided normalisation function
  } else {
    # Default: Normalize the counts using TMM normalization
    y <- edgeR::calcNormFactors(y)
  }
  
  # Step 3: Apply voom transformation to get logCPM values and weights
  voom_obj <- limma::voom(
    y,
    design_matrix
  )
  
  return(voom_obj)
}


#' Process Top Table
#'
#' @description
#' Processes the top table from a LIMMA analysis, adding feature names and 
#' intercepts.
#'
#' @param top_table_and_fit A list containing the top table and the fit object 
#' from LIMMA.
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
    feature_names
    ) {

  top_table <- process_within_level_result$top_table
  fit <- process_within_level_result$fit

  top_table <- modify_limma_top_table(
    top_table,
    feature_names
    )
  
  intercepts <- as.data.frame(stats::coef(fit)[, "(Intercept)", drop = FALSE])
  intercepts_ordered <- intercepts[match(top_table$feature_nr, 
                                         rownames(intercepts)), , 
                                   drop = FALSE]
  top_table$intercept <- intercepts_ordered[, 1]
  
  top_table
}


#' Process Within Level
#'
#' @description
#' Performs a within-level analysis using limma to generate top tables and fit 
#' objects based on the specified spline parameters. Performs the limma spline 
#' analysis for a selected level of a factor
#'
#' @param data A matrix of data values.
#' @param preprocess_rna_seq Boolean specifying whether to preprocess RNA seq
#' @param normalization_fun Function for normalizing RNA-seq raw counts.
#' @param meta A dataframe containing metadata, including a 'Time' column.
#' @param design A design formula or matrix for the limma analysis.
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
#' @importFrom splines bs
#' @importFrom splines ns
#' @importFrom stats as.formula
#' @importFrom stats model.matrix
#' @importFrom limma lmFit
#' @importFrom limma eBayes
#' @importFrom limma topTable
#' 
process_within_level <- function(
    data,
    preprocess_rna_seq,
    normalization_fun,
    meta,
    design,
    spline_params,
    level_index,
    padjust_method
    ) {

  design_matrix <- design2design_matrix(
    meta,
    spline_params,
    level_index,
    design
    )

  if (preprocess_rna_seq) {
    data <- preprocess_rna_seq_data(
      raw_counts = data,
      design_matrix = design_matrix,
      normalize_func = normalization_fun 
    )
    voom_data_matrix <- data$E
  } else {
    voom_data_matrix = NULL
  }
  
  fit <- limma::lmFit(
    data,
    design_matrix
    )
  fit <- limma::eBayes(fit)

  num_matching_columns <- sum(grepl(
    "^X\\d+$",
    colnames(design_matrix)
    ))
  coeffs <- paste0("X", seq_len(num_matching_columns))
  
  top_table <- limma::topTable(
    fit,
    adjust.method = padjust_method, 
    number = Inf,
    coef = coeffs
    )

  attr(top_table, "adjust.method") <- padjust_method

  list(
    top_table = top_table,
    fit = fit,
    voom_data_matrix = voom_data_matrix
    )
}



# Level 3 internal functions ---------------------------------------------------


#' Modify limma Top Table
#'
#' @description
#' Modifies the limma top table to include feature indices and names.
#'
#' @param top_table A dataframe containing the top table results from limma
#' @param feature_names A character vector of feature names.
#'
#' @return A tibble with feature indices and names included.
#' 
#' @importFrom tidyr as_tibble
#' @importFrom dplyr relocate last_col mutate
#' 
modify_limma_top_table <- function(
    top_table, 
    feature_names
    ) {

  top_table <- tidyr::as_tibble(
    top_table, 
    rownames = "feature_nr"
    )
  
  feature_nr <- NULL  # dummy declaration for the lintr and R CMD.

  # Convert feature_nr to integer
  top_table <- top_table |> 
    dplyr::mutate(feature_nr = as.integer(feature_nr)) |>
    dplyr::relocate(feature_nr, .after = dplyr::last_col())
  
  # Sort and add feature names based on the feature_nr
  sorted_feature_names <- feature_names[top_table$feature_nr]
  top_table <- top_table |> dplyr::mutate(feature_names = sorted_feature_names)
  
  return(top_table)
}
