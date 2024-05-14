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
#' @param data A matrix or dataframe containing the expression data, 
#'             where rows represent features (e.g., genes) and columns represent 
#'             samples.
#' @param meta A dataframe containing the metadata for the samples, 
#'             must include a column for each factor in `factors`.
#' @param design A character string specifying the design formula for the Limma 
#' model.
#' @param DoFs An integer vector specifying the degrees of freedom for spline 
#' interpolation 
#'             for each factor; length should match that of `factors`.
#' @param factors A character vector specifying the names of experimental 
#' factors 
#'                to be analyzed, as represented in `meta`.
#' @param feature_names A character vector of feature names (e.g., gene names) 
#' corresponding 
#'                      to the rows in `data`.
#' @param mode A character string specifying the analysis mode: 'isolated' for 
#' analyzing 
#'             each level of a factor separately, or 'integrated' for comparing 
#'             levels 
#'             within a factor.
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
#' @importFrom purrr partial
#' @importFrom purrr map
#' @importFrom purrr map_chr
#' @importFrom purrr map2
#' @importFrom stats setNames
#' @importFrom utils combn
#' 
#' @export
#' 
run_limma_splines <- function(data,
                              meta,
                              design,
                              spline_params,
                              condition, 
                              feature_names,
                              mode = c("isolated", "integrated"),
                              padjust_method = "BH") {
  mode <- match.arg(mode)
  control_inputs_run_limma(data = data, 
                           meta = meta, 
                           design = design, 
                           spline_params = spline_params, 
                           condition = condition, 
                           feature_names = feature_names, 
                           mode = mode, 
                           padjust_method = padjust_method)
  
  meta[[condition]] <- factor(meta[[condition]])
  levels <- levels(meta[[condition]])
  
  # Get hits for level (within level analysis) ---------------------------------
  process_level_with_params <- purrr::partial(process_level, 
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
  
  top_tables <- stats::setNames(purrr::map(results_list, "top_table"), 
                                purrr::map_chr(results_list, "name"))
  
  
  # Factor and Factor:Time comparisons between levels --------------------------
  ttslc_factor_only <- list()        # ttslc = top_tables_level_comparison
  ttslc_factor_time <- list()        # Factor AND time
  
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
      
      ttslc_factor_only[[paste0(lev_combo[1], "_vs_", lev_combo[2])]] <- 
        result$ttlc_factor_only
      ttslc_factor_time[[paste0(lev_combo[1], "_vs_", lev_combo[2])]] <- 
        result$ttlc_factor_time
    }
  } else { # mode == "isolated"
    print("mode == 'integrated' necessary for between level comparisons. 
          Returning emtpy lists for ttslc_factor_only and ttslc_factor_time
          (ttslc means 'top tables level comparison').")
  }
  
  list(top_tables = top_tables, 
       ttslc_factor_only = ttslc_factor_only,
       ttslc_factor_time = ttslc_factor_time)
}



# Level 1 internal functions ---------------------------------------------------


#' Controls the function input and raises an error if not in expected format
control_inputs_run_limma <- function(data, 
                                     meta, 
                                     design, 
                                     spline_params, 
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
  
  if (!(nrow(meta) == ncol(data))) {
    stop("The number of rows in meta does not match the number of columns in 
        data. This is required.\n")
  }
  
  # Check that design is a single string
  if (!is.character(design) || length(design) != 1) {
    stop("design must be a single string.")
  }
  
  check_spline_params(spline_params, mode)
  
  # Ensure factors is a non-empty character vector
  if (!is.character(condition) && !(length(condition) == 1)) {
    stop("factors must be a non-empty character vector.")
  }
  
  # Check that feature_names is a non-empty character vector
  if (!is.character(feature_names) || length(feature_names) == 0) {
    stop("feature_names must be a non-empty character vector.")
  }

  supported_methods <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", 
                         "fdr", "none")
  if (!(is.character(padjust_method) && padjust_method %in% supported_methods)) 
  {
    stop("padjust_method must be a character and one of the supported methods (
         holm, hochberg, hommel, bonferroni, BH, BY, fdr, none).")
  }
}


#' Calculate Differential Expression Between Levels
#'
#' Calculates differential expression between specified levels of a condition
#' using linear models and limma analysis. Returns top table results for 
#' contrasts between factors only and factors with time interactions.
#'
#' @param data A matrix of expression values.
#' @param meta Metadata containing experimental design information.
#' @param design A formula specifying the linear model design.
#' @param spline_params Parameters for spline generation.
#' @param condition The name of the condition column in the metadata.
#' @param compared_levels The levels of the condition to compare.
#' @param padjust_method The method for adjusting p-values.
#' @param feature_names Names of the features in the data matrix.
#' @return A list containing top tables for contrasts between factors only 
#'         and factors with time interactions.
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
  
  args <- list(x = meta$Time, intercept = FALSE)
  
  if (!is.null(spline_params$DoFs)) {
    args$df <- spline_params$DoFs[1]
  } else {
    args$knots <- spline_params$knots[[1]]
  }
  
  if (!is.null(spline_params$bknots)) {
    args$Boundary.knots <- spline_params$bknots[[1]]
  }
  
  
  if (spline_params$spline_type[1] == "b") {
    args$degree <- spline_params$degrees[1]
    meta$X <- do.call(splines::bs, args)
  } else {                                          # natural cubic splines
    meta$X <- do.call(splines::ns, args)
  }
  
  design_matrix <- stats::model.matrix(stats::as.formula(design), data = meta)
  
  fit <- limma::lmFit(data, design_matrix)
  fit <- limma::eBayes(fit)
  
  factor_only_contrast_coeff <- paste0(condition, compared_levels[2])
  ttlc_factor_only <- limma::topTable(fit, coef = factor_only_contrast_coeff,
                               adjust.method = padjust_method, number = Inf)
  ttlc_factor_only <- modify_limma_top_table(ttlc_factor_only, feature_names)
  
  num_matching_columns <- sum(grepl("^X\\d+$", colnames(design_matrix)))
  factor_time_contrast_coeffs <- paste0(condition, compared_levels[2], ":X", 
                                        seq_len(num_matching_columns))
  ttlc_factor_time <- limma::topTable(fit, coef = factor_time_contrast_coeffs,
                               adjust.method = padjust_method, number = Inf)
  ttlc_factor_time <- modify_limma_top_table(ttlc_factor_time, feature_names)
  
  list(ttlc_factor_only = ttlc_factor_only, ttlc_factor_time = ttlc_factor_time)
}


#' Converts top_table to a tibble and adds the columns feature_names and 
#' (Intercept) and appends it to the top_table list.
#' 
#' @importFrom stats coef
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

#' @importFrom stats relevel
process_level <- function(level, 
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
  
  result <- within_level(data_copy, 
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


#' Converts top_table to a tibble and adds new column feature_name
#' @importFrom tidyr as_tibble
#' @importFrom dplyr relocate
#' @importFrom dplyr mutate
modify_limma_top_table <- function(top_table, 
                                   feature_names) {
  top_table <- tidyr::as_tibble(top_table, rownames = "feature_index") %>% 
    dplyr::relocate(feature_index, .after = last_col()) %>%
    dplyr::mutate(feature_index = as.integer(feature_index))
  
  sorted_feature_names <- feature_names[top_table$feature_index]
  top_table <- top_table %>% dplyr::mutate(feature_names = sorted_feature_names)
  
  top_table
}


#' Performs the limma spline analysis for a selected level of a factor
#' @importFrom splines bs
#' @importFrom splines ns
#' @importFrom stats as.formula
#' @importFrom stats model.matrix
#' @importFrom limma lmFit
#' @importFrom limma eBayes
#' @importFrom limma topTable
#' 
within_level <- function(data,
                         meta,
                         design,
                         factor,
                         level,
                         spline_params,
                         level_index,
                         padjust_method) {
  
  args <- list(x = meta$Time, intercept = FALSE)
  
  if (!is.null(spline_params$DoFs)) {
    args$df <- spline_params$DoFs[level_index]
  } else {
    args$knots <- spline_params$knots[[level_index]]
  }
  
  if (!is.null(spline_params$bknots)) {
    args$Boundary.knots <- spline_params$bknots[[level_index]]
  }
  
  
  if (spline_params$spline_type[level_index] == "b") {
    args$degree <- spline_params$degrees[level_index]
    meta$X <- do.call(splines::bs, args)
  } else {                                          # natural cubic splines
    meta$X <- do.call(splines::ns, args)
  }
  
  design_matrix <- stats::model.matrix(stats::as.formula(design), data = meta)
  
  fit <- limma::lmFit(data, design_matrix)
  fit <- limma::eBayes(fit)

  column_names <- colnames(design_matrix)
  num_matching_columns <- sum(grepl("^X\\d+$", colnames(design_matrix)))
  coeffs <- paste0("X", seq_len(num_matching_columns))
  
  top_table <- limma::topTable(fit, adjust=padjust_method, number=Inf, 
                               coef=coeffs)
  
  list(top_table = top_table, fit = fit)
}
