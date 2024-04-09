# Import libraries -------------------------------------------------------------
library(limma)
library(splines)



# Internal functions level 2 ----------------------------------------------


modify_limma_top_table <- function(top_table, feature_names) {
  top_table <- as_tibble(top_table, rownames = "feature_index") %>% 
    select(-feature_index, everything(), feature_index) %>%
    mutate(feature_index = as.numeric(feature_index))
  
  sorted_feature_names <- feature_names[top_table$feature_index]
  top_table <- top_table %>%
    mutate(feature_names = sorted_feature_names)
  
  return(top_table)
}


# Internal functions level 1 ---------------------------------------------------


control_inputs_run_limma <- function(data, meta, design, DoFs, factors, 
                                     feature_names, mode) {
  # Check if data is a dataframe and not a matrix, then convert to matrix
  if (is.data.frame(data) && !is.matrix(data)) {
    data <- as.matrix(data)
  }
  
  # Ensure meta is a dataframe with a 'Time' column
  if (!is.data.frame(meta) || !"Time" %in% names(meta)) {
    stop("meta must be a dataframe and contain the column 'Time'.")
  }
  
  # Check that design is a single string
  if (!is.character(design) || length(design) != 1) {
    stop("design must be a single string.")
  }
  
  # Validate DoFs as an integer vector
  if (!is.integer(DoFs)) {
    stop("DoFs must be an integer vector.")
  }
  
  # Ensure factors is a non-empty character vector
  if (!is.character(factors) || length(factors) == 0) {
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
}


between_level <- function(data, meta, design, DoF, factor, level2,
                          padjust_method, feature_names) {
  meta$X <- ns(meta$Time, df=DoF, intercept = FALSE)
  design_matrix <- model.matrix(as.formula(design), data = meta)

  fit <- lmFit(data, design_matrix)
  fit <- eBayes(fit)

  factor_only_contrast_coeff <- paste0(factor, level2)
  ttlc_factor_only <- topTable(fit, coef=factor_only_contrast_coeff,
                                   adjust.method=padjust_method, number=Inf)
  ttlc_factor_only <- modify_limma_top_table(ttlc_factor_only, feature_names)
  
  factor_time_contrast_coeffs <- paste0(factor, level2, ":X", seq_len(DoF))
  ttlc_factor_time <- topTable(fit, coef=factor_time_contrast_coeffs,
                               adjust.method=padjust_method, number=Inf)
  ttlc_factor_time <- modify_limma_top_table(ttlc_factor_time, feature_names)
  
  return(list(ttlc_factor_only = ttlc_factor_only, 
              ttlc_factor_time = ttlc_factor_time))
}


within_level <- function(data, meta, design, factor, level, DoF, 
                         padjust_method) {
  meta$X <- ns(meta$Time, df=DoF, intercept = FALSE)
  design_matrix <- model.matrix(as.formula(design), data = meta)
  
  fit <- lmFit(data, design_matrix)
  fit <- eBayes(fit)
  
  coef_names <- paste0("X", seq_len(DoF))

  top_table <- topTable(fit, adjust=padjust_method, number=Inf, coef=coef_names)
  
  return(list(top_table = top_table, fit = fit))
}


process_top_table <- function(top_table, top_tables, feature_names, fit, factor, 
                              level) {
  top_table <- modify_limma_top_table(top_table, feature_names)
  
  intercepts <- as.data.frame(coef(fit)[, "(Intercept)", drop = FALSE])
  intercepts_ordered <- intercepts[match(top_table$feature_index, 
                                         rownames(intercepts)), , 
                                   drop = FALSE]
  top_table$intercept <- intercepts_ordered[, 1]
  
  results_name <- paste(factor, level, sep = "_")
  top_tables[[results_name]] <- top_table
  
  return(top_tables)
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
#' @importFrom dplyr %>%
#' @importFrom stats setNames
#' 
run_limma_splines <- function(data, meta, design, DoFs, factors, 
                              feature_names, mode, padjust_method = "BH") {
  
  control_inputs_run_limma(data, meta, design, DoFs, factors, feature_names, 
                           mode)
  
  top_tables <- list()
  ttslc_factor_only <- list()        # ttslc = top_tables_level_comparison
  ttslc_factor_time <- list()        # Factor AND time
  
  for (factor in factors) {             # for example fermentation phase
    meta[[factor]] <- factor(meta[[factor]])
    levels <- levels(meta[[factor]])
    
    # Within level analysis ----------------------------------------------------
    for (i in 1:(length(levels))) { 
      level <- levels[i]  # level is for example exponential or stationary phase
      
      DoF <- DoFs[i]
      if (is.na(DoF)) {
        DoF <- DoFs[1]
      }
      
      if (mode == "isolated") {
        samples <- which(meta[[factor]] == level)
        data_copy <- data[, samples]
        meta_copy <- subset(meta, meta[[factor]] == level)
      } else if (mode == "integrated") {
        data_copy <- data
        meta_copy <- meta
        meta_copy[[factor]] <- relevel(meta_copy[[factor]], ref = level)
        
        level2 <- levels[i + 1]
        if (!is.na(level2)) {
          result <- between_level(data, meta, design, DoF, factor, level2, 
                          padjust_method, feature_names)
          
          ttslc_factor_only[[paste0(level, "_vs_", level2)]] <- 
            result$ttlc_factor_only
          ttslc_factor_time[[paste0(level, "_vs_", level2)]] <- 
            result$ttlc_factor_time
        }
      }
      
      result <- within_level(data_copy, meta_copy, design, factor, level, DoF, 
                             padjust_method)
      top_table <- result$top_table
      fit <- result$fit
      
      top_tables <- process_top_table(top_table, top_tables, feature_names,
                                      fit, factor, level)
    }
  }
  
  if (mode == "isolated") {
    print("mode == 'integrated' necessary for between level comparison. 
          Returning emtpy list for ttslc_factor_only and ttslc_factor_time
          (ttslc = top_tables_level_comparison).")
  }
  
  return(list(top_tables = top_tables, 
              ttslc_factor_only = ttslc_factor_only,
              ttslc_factor_time = ttslc_factor_time))
}
