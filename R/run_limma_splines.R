# Import libraries -------------------------------------------------------------
library(limma)
library(splines)
library(purrr)



# Internal functions level 2 ---------------------------------------------------


#' Converts top_table to a tibble and adds new column feature_name
modify_limma_top_table <- function(top_table, 
                                   feature_names) {
  top_table <- as_tibble(top_table, rownames = "feature_index") %>% 
    relocate(feature_index, .after = last_col()) %>%
    mutate(feature_index = as.integer(feature_index))
  
  sorted_feature_names <- feature_names[top_table$feature_index]
  top_table <- top_table %>% mutate(feature_names = sorted_feature_names)
  
  top_table
}


#' Performs the limma spline analysis for a selected level of a factor
within_level <- function(data, 
                         meta, 
                         design, 
                         factor, 
                         level, 
                         DoF, 
                         padjust_method) {
  meta$X <- ns(meta$Time, df=DoF, intercept = FALSE)
  design_matrix <- model.matrix(as.formula(design), data = meta)
  
  fit <- lmFit(data, design_matrix)
  fit <- eBayes(fit)
  
  coeffs <- paste0("X", seq_len(DoF))
  
  top_table <- topTable(fit, adjust=padjust_method, number=Inf, coef=coeffs)
  
  list(top_table = top_table, fit = fit)
}



# Internal functions level 1 ---------------------------------------------------


#' Controls the function input and raises an error if not in expected format
control_inputs_run_limma <- function(data, 
                                     meta, 
                                     design, 
                                     DoFs, 
                                     condition, 
                                     feature_names, 
                                     mode, 
                                     padjust_method) {
  if (!is.matrix(data)) {
    stop("data must be a matrix")
  }
  
  # Ensure meta is a dataframe with a 'Time' column
  if (!is.data.frame(meta) || 
      !"Time" %in% names(meta) || 
      !is.numeric(meta$Time)) {
    stop("meta must be a dataframe and contain the column 'Time' 
         with numeric values.")
  }
  
  # Check that design is a single string
  if (!is.character(design) || length(design) != 1) {
    stop("design must be a single string.")
  }
  
  # Validate DoFs as an integer vector
  if (!is.integer(DoFs) || length(DoFs) != length(unique(meta[[condition]]))) {
    stop("DoFs must be an integer vector and its length must match the number of 
         levels in meta$condition")
  }
  
  # Ensure factors is a non-empty character vector
  if (!is.character(condition) && !(length(condition) == 1)) {
    stop("factors must be a non-empty character vector.")
  }
  
  # Check that feature_names is a non-empty character vector
  if (!is.character(feature_names) || length(feature_names) == 0) {
    stop("feature_names must be a non-empty character vector.")
  }
  
  # Ensure mode is either 'integrated' or 'isolated'
  if (!(mode == "integrated" || mode == "isolated")) {
    stop("mode must be either 'integrated' or 'isolated'.")
  }
  
  supported_methods <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", 
                         "fdr", "none")
  if (!(is.character(padjust_method) && padjust_method %in% supported_methods)) 
    {
    stop("padjust_method must be a character and one of the supported methods (
         holm, hochberg, hommel, bonferroni, BH, BY, fdr, none).")
  }
}


#' Performs limma spline analysis and gets the top_table for the level 
#' comparison: Both for factor and for factor:time
between_level <- function(data, 
                          meta, 
                          design, 
                          DoFs, 
                          condition, 
                          compared_levels,
                          padjust_method, 
                          feature_names) {
  samples <- which(meta[[condition]] %in% compared_levels)
  data <- data[, samples]
  meta <- subset(meta, meta[[condition]] %in% compared_levels)
  
  meta$X <- ns(meta$Time, df = DoFs[1], intercept = FALSE)
  design_matrix <- model.matrix(as.formula(design), data = meta)

  fit <- lmFit(data, design_matrix)
  fit <- eBayes(fit)

  factor_only_contrast_coeff <- paste0(condition, compared_levels[2])
  ttlc_factor_only <- topTable(fit, coef = factor_only_contrast_coeff,
                               adjust.method = padjust_method, number = Inf)
  ttlc_factor_only <- modify_limma_top_table(ttlc_factor_only, feature_names)
  
  factor_time_contrast_coeffs <- paste0(condition, compared_levels[2], ":X", 
                                        seq_len(DoFs[2]))
  ttlc_factor_time <- topTable(fit, coef = factor_time_contrast_coeffs,
                               adjust.method = padjust_method, number = Inf)
  ttlc_factor_time <- modify_limma_top_table(ttlc_factor_time, feature_names)
  
  list(ttlc_factor_only = ttlc_factor_only, ttlc_factor_time = ttlc_factor_time)
}


#' Converts top_table to a tibble and adds the columns feature_names and 
#' (Intercept) and appends it to the top_table list.
process_top_table <- function(top_table_and_fit, 
                              feature_names) {
  top_table <- top_table_and_fit$top_table
  fit <- top_table_and_fit$fit
  
  top_table <- modify_limma_top_table(top_table, feature_names)
  
  intercepts <- as.data.frame(coef(fit)[, "(Intercept)", drop = FALSE])
  intercepts_ordered <- intercepts[match(top_table$feature_index, 
                                         rownames(intercepts)), , 
                                   drop = FALSE]
  top_table$intercept <- intercepts_ordered[, 1]
  
  top_table
}


process_level <- function(level, 
                          DoF, 
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
    meta_copy[[condition]] <- relevel(meta_copy[[condition]], ref = level)
  }
  
  result <- within_level(data_copy, 
                         meta_copy, 
                         design, 
                         condition, 
                         level, 
                         DoF, 
                         padjust_method)
  top_table <- process_top_table(result, 
                                 feature_names)
  
  results_name <- paste(condition, level, sep = "_")
  list(name = results_name, top_table = top_table)
}



# Export functions -------------------------------------------------------------


#' Run Limma Analysis with Spline Interpolation for Hyperparameter Screening
#'
#' This function conducts differential expression analysis using the Limma package, 
#' incorporating spline interpolation to model the effect of various experimental 
#' factors across different levels. It supports both isolated and integrated modes 
#' for within-level analysis and between-level comparison, adjusting for multiple 
#' degrees of freedom corresponding to the factors under investigation.
#'
#' @param data A matrix or dataframe containing the expression data, 
#'             where rows represent features (e.g., genes) and columns represent samples.
#' @param meta A dataframe containing the metadata for the samples, 
#'             must include a column for each factor in `factors`.
#' @param design A character string specifying the design formula for the Limma model.
#' @param DoFs An integer vector specifying the degrees of freedom for spline interpolation 
#'             for each factor; length should match that of `factors`.
#' @param factors A character vector specifying the names of experimental factors 
#'                to be analyzed, as represented in `meta`.
#' @param feature_names A character vector of feature names (e.g., gene names) corresponding 
#'                      to the rows in `data`.
#' @param mode A character string specifying the analysis mode: 'isolated' for analyzing 
#'             each level of a factor separately, or 'integrated' for comparing levels 
#'             within a factor.
#' @param padjust_method A character string specifying the method for adjusting 
#'                       p-values for multiple testing. Defaults to "BH" (Benjamini-Hochberg).
#'
#' @return A list containing three elements: 
#'         - `top_tables`: A list of top tables generated from within-level analysis 
#'                         for each factor and level.
#'         - `ttslc_factor_only`: A list of top tables from between-level comparisons 
#'                                for each factor, excluding time effects.
#'         - `ttslc_factor_time`: A list of top tables from between-level comparisons 
#'                                for each factor, including time effects.
#'
#' @examples
#' \dontrun{
#'   results <- run_limma_splines(data, meta, "~ 1 + Factor*X + Time", c(2L), 
#'                                c("Factor"), c("Gene1", "Gene2"), "isolated")
#' }
#'
#' @export
#' @importFrom limma lmFit eBayes topTable
#' people should use R 4.1, so that they can use the native pipe ( %>% )
#' @importFrom stats setNames
#' 
run_limma_splines <- function(data,
                              meta,
                              design,
                              DoFs,
                              condition, 
                              feature_names,
                              mode = c("isolated", "integrated"),
                              padjust_method = "BH") {
  mode <- match.arg(mode)
  control_inputs_run_limma(data, 
                           meta, 
                           design, 
                           DoFs, 
                           condition, 
                           feature_names, 
                           mode, 
                           padjust_method)
  
  meta[[condition]] <- factor(meta[[condition]])
  levels <- levels(meta[[condition]])
  
  # Get hits for level (within level analysis) ---------------------------------
  process_with_data_meta <- purrr::partial(process_level, 
                                           data = data, 
                                           meta = meta, 
                                           design = design, 
                                           condition = condition, 
                                           feature_names = feature_names, 
                                           padjust_method = padjust_method, 
                                           mode = mode)
  
  results_list <- map2(levels, DoFs, process_with_data_meta)
  top_tables <- setNames(map(results_list, "top_table"), 
                         map_chr(results_list, "name"))
  
  
  # Factor and Factor:Time comparisons between levels --------------------------
  ttslc_factor_only <- list()        # ttslc = top_tables_level_comparison
  ttslc_factor_time <- list()        # Factor AND time
  
  if (mode == "integrated") {
    level_combinations <- combn(levels, 2, simplify = FALSE)
    for (lev in level_combinations) {
      lev_DoFs <- DoFs[match(lev, unique(meta[[condition]]))]
      
      result <- between_level(data, 
                              meta, 
                              design, 
                              lev_DoFs, 
                              condition, 
                              lev, 
                              padjust_method, 
                              feature_names)
      
      ttslc_factor_only[[paste0(lev[1], "_vs_", lev[2])]] <- 
        result$ttlc_factor_only
      ttslc_factor_time[[paste0(lev[1], "_vs_", lev[2])]] <- 
        result$ttlc_factor_time
    }
  } else { # mode == "isolated"
    print("mode == 'integrated' necessary for between level comparison. 
          Returning emtpy lists for ttslc_factor_only and ttslc_factor_time
          (ttslc means top tables level comparison).")
  }
  
  list(top_tables = top_tables, 
       ttslc_factor_only = ttslc_factor_only,
       ttslc_factor_time = ttslc_factor_time)
}
