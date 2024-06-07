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
#' @param data A dataframe containing the expression data, 
#'             where rows represent features (e.g., genes) and columns represent 
#'             samples. 
#' @param meta A dataframe containing the metadata for the samples, 
#'             must include a column for each factor in `factors`.
#' @param design A character string specifying the design formula for the Limma 
#' model.
#' @param condition The name of the factor column of meta. The factor divides
#' the experiment into different levels (unique elements of the factor column). 
#' @param spline_params A named list, specifying the parameters for the splines,
#' such as type (ns or B-splines, dof, knots, etc).
#' @param padjust_method A character string specifying the method for adjusting 
#'                       p-values for multiple testing. Defaults to "BH" 
#'                       (Benjamini-Hochberg).
#'
#' @return A list containing three elements: 
#'         - `top_tables`: A list of top tables generated from within-level 
#'         analysis 
#'                         for each factor and level.
#'         - `ttslc_factor_only`: A list of top tables from between-level 
#'         comparisons 
#'                                for each factor, excluding time effects.
#'         - `ttslc_factor_time`: A list of top tables from between-level 
#'         comparisons 
#'                                for each factor, including time effects.
#'
#' @examples
#' \dontrun{
#'   results <- run_limma_splines(data, meta, "~ 1 + Factor*X + Time", c(2L), 
#'                                c("Factor"), c("Gene1", "Gene2"), "isolated")
#' }
#' 
#' @importFrom purrr partial map map_chr map2
#' @importFrom stats setNames
#' @importFrom utils combn
#' @importFrom stringr str_split
#' 
#' @export
#' 
run_limma_splines <- function(data,
                              meta,
                              design,
                              condition, 
                              spline_params = list(spline_type = c("n"),
                                                   dof = c(2L)),
                              padjust_method = "BH") {
  
  mode <- determine_analysis_mode(design,
                                  condition)
  
  args <- lapply(as.list(match.call()[-1]), eval, parent.frame())
  args$mode <- mode
  check_null_elements(args)
  input_control <- InputControl$new(args)
  input_control$auto_validate()

  matrix_and_feature_names <- process_data(data)
  data <- matrix_and_feature_names$data
  feature_names <- matrix_and_feature_names$feature_names
  
  meta[[condition]] <- factor(meta[[condition]])
  levels <- levels(meta[[condition]])
  
  # Get hits for level (within level analysis) ---------------------------------
  process_level_with_params <- purrr::partial(within_level, 
                                              spline_params = spline_params,
                                              data = data, 
                                              meta = meta, 
                                              design = design, 
                                              condition = condition, 
                                              feature_names = feature_names, 
                                              padjust_method = padjust_method, 
                                              mode = mode)
  
  results_list <- purrr::map2(levels, 
                              seq_along(levels), 
                              process_level_with_params)
  
  within_level_top_table <- 
    stats::setNames(purrr::map(results_list, "top_table"), 
                    purrr::map_chr(results_list, "name"))
  
  
  # Factor and Factor:Time comparisons between levels --------------------------
  between_level_condition_only <- list()        
  between_level_condition_time <- list()        # Factor AND time
  
  if (mode == "integrated") {
    level_combinations <- utils::combn(levels, 2, simplify = FALSE)
    for (lev_combo in level_combinations) {
      result <- between_level(data = data, 
                              meta = meta, 
                              design = design, 
                              spline_params = spline_params,
                              condition = condition, 
                              compared_levels = lev_combo, 
                              padjust_method = padjust_method, 
                              feature_names = feature_names)
      
      between_level_condition_only[[paste0("avrg_diff_" ,lev_combo[1],
                                           "_vs_", lev_combo[2])]] <- 
        result$condition_only
      between_level_condition_time[[paste0("time_interaction_" ,
                                           lev_combo[1],
                                           "_vs_", lev_combo[2])]] <- 
        result$condition_time
    }
  } else { # mode == "isolated"
    message(paste0("mode == 'integrated' necessary for between level ",
                 "comparisons. Returning emtpy lists for ttslc_factor_only ",
                 "and ttslc_factor_time (ttslc means 'top tables level ",
                 "comparison')."))
  }
  
 list(time_effect = within_level_top_table, 
      avrg_diff_conditions = between_level_condition_only,
      interaction_condition_time = between_level_condition_time)
}



# Level 1 internal functions ---------------------------------------------------


#' Between Level Analysis
#'
#' @description
#' Performs a between-level analysis using LIMMA to compare specified levels 
#' within a condition.
#'
#' @param data A matrix of data values.
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
between_level <- function(data, 
                          meta, 
                          design, 
                          spline_params, 
                          condition, 
                          compared_levels,
                          padjust_method, 
                          feature_names) {
  
  samples <- which(meta[[condition]] %in% compared_levels)
  data <- data[, samples]
  meta <- subset(meta, meta[[condition]] %in% compared_levels)
  
  design_matrix <- design2design_matrix(meta = meta,
                                        spline_params = spline_params,
                                        level_index = 1,
                                        design = design)
  
  fit <- limma::lmFit(data, design_matrix)
  fit <- limma::eBayes(fit)

  factor_only_contrast_coeff <- paste0(condition, compared_levels[2])
  condition_only <- limma::topTable(fit, coef = factor_only_contrast_coeff,
                               adjust.method = padjust_method, number = Inf)
  
  condition_only_resuls <- list(top_table = condition_only,
                                fit = fit)
  top_table_condition_only <- process_top_table(condition_only_resuls, 
                                                feature_names)
  
  
  num_matching_columns <- sum(grepl("^X\\d+$", colnames(design_matrix)))
  factor_time_contrast_coeffs <- paste0(condition, compared_levels[2], ":X", 
                                        seq_len(num_matching_columns))
  condition_time <- limma::topTable(fit, coef = factor_time_contrast_coeffs,
                                    adjust.method = padjust_method,
                                    number = Inf)
  
  condition_and_time_results <- list(top_table = condition_time,
                                     fit = fit)
  top_table_condition_and_time <- process_top_table(condition_and_time_results, 
                                                    feature_names)
  
  list(condition_only = top_table_condition_only,
       condition_time = top_table_condition_and_time)
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
within_level <- function(level, 
                         level_index,
                         spline_params,
                         data, 
                         meta, 
                         design, 
                         condition, 
                         feature_names, 
                         padjust_method, 
                         mode) {
  
  if (mode == "isolated") {
    samples <- which(meta[[condition]] == level)
    data_copy <- data[, samples]
    meta_copy <- subset(meta, meta[[condition]] == level)
  } else { # mode == "integrated"
    data_copy <- data
    meta_copy <- meta
    meta_copy[[condition]] <- stats::relevel(meta_copy[[condition]], 
                                             ref = level)
    level_index <- 1L   # spline_params must be uniform across all levels for
    # integrated mode.
  }
  
  result <- process_within_level(data_copy, 
                                 meta_copy, 
                                 design, 
                                 condition, 
                                 level,
                                 spline_params,
                                 level_index, 
                                 padjust_method)

  top_table <- process_top_table(result, 
                                 feature_names)
  
  results_name <- paste(condition, level, sep = "_")
  list(name = results_name, top_table = top_table)
}


# Level 2 internal functions ---------------------------------------------------


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
process_top_table <- function(top_table_and_fit, 
                              feature_names) {
  
  top_table <- top_table_and_fit$top_table
  fit <- top_table_and_fit$fit
  
  top_table <- modify_limma_top_table(top_table, feature_names)
  
  intercepts <- as.data.frame(stats::coef(fit)[, "(Intercept)", drop = FALSE])
  intercepts_ordered <- intercepts[match(top_table$feature_index, 
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
#' @param meta A dataframe containing metadata, including a 'Time' column.
#' @param design A design formula or matrix for the limma analysis.
#' @param factor A character string specifying the factor.
#' @param level The level within the factor to process.
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
process_within_level <- function(data,
                                 meta,
                                 design,
                                 factor,
                                 level,
                                 spline_params,
                                 level_index,
                                 padjust_method) {

  design_matrix <- design2design_matrix(meta,
                                        spline_params,
                                        level_index,
                                        design)
  
  fit <- limma::lmFit(data, design_matrix)
  fit <- limma::eBayes(fit)

  num_matching_columns <- sum(grepl("^X\\d+$", colnames(design_matrix)))
  coeffs <- paste0("X", seq_len(num_matching_columns))
  
  top_table <- limma::topTable(fit, adjust.method = padjust_method, 
                               number = Inf, coef = coeffs)

  attr(top_table, "adjust.method") <- padjust_method
  
  list(top_table = top_table, fit = fit)
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
#' @seealso
#' \link[tidyr]{as_tibble}, \link[dplyr]{relocate}, \link[dplyr]{mutate}
#' 
#' @importFrom tidyr as_tibble
#' @importFrom dplyr relocate last_col mutate
#' 
modify_limma_top_table <- function(top_table, 
                                   feature_names) {
  
  feature_index <- dplyr::sym("feature_index")
  
  top_table <- tidyr::as_tibble(top_table, rownames = "feature_index")
  
  top_table <- top_table %>% 
    dplyr::relocate(feature_index, .after = dplyr::last_col()) %>%
    dplyr::mutate(feature_index = as.integer(feature_index))
  
  sorted_feature_names <- feature_names[top_table$feature_index]
  top_table <- top_table %>% dplyr::mutate(feature_names = sorted_feature_names)
  
  top_table
}
