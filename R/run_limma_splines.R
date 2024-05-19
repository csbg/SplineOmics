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
                              condition, 
                              feature_names,
                              mode = "integrated",
                              spline_params = list(spline_type = c("n"),
                                                   dof = c(2L)),
                              padjust_method = "BH") {
  
  control_inputs_run_limma(data = data, 
                           meta = meta, 
                           design = design, 
                           condition = condition, 
                           feature_names = feature_names, 
                           spline_params = spline_params, 
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
    message(paste0("mode == 'integrated' necessary for between level ",
                 "comparisons. Returning emtpy lists for ttslc_factor_only ",
                 "and ttslc_factor_time (ttslc means 'top tables level ",
                 "comparison')."))
  }
  
  list(top_tables = top_tables, 
       ttslc_factor_only = ttslc_factor_only,
       ttslc_factor_time = ttslc_factor_time)
}



# Level 1 internal functions ---------------------------------------------------


#' Control Inputs for Running LIMMA
#'
#' @description
#' Validates the inputs for running LIMMA analysis, ensuring all required 
#' structures and parameters are correctly formatted.
#'
#' @param data A matrix of data values.
#' @param meta A dataframe containing metadata, including a 'Time' column with 
#' numeric values.
#' @param design A design formula or matrix for the LIMMA analysis.
#' @param spline_params A list of spline parameters for the analysis.
#' @param condition A character string specifying the condition.
#' @param feature_names A non-empty character vector of feature names.
#' @param mode A character string specifying the mode 
#'            ('isolated' or 'integrated').
#' @param padjust_method A character string specifying the p-adjustment method.
#'
#' @return No return value, called for side effects.
#'
#' @examples
#' data <- matrix(runif(100), nrow = 10)
#' meta <- data.frame(Time = seq(1, 10))
#' design <- ~ 1
#' spline_params <- list(df = 3)
#' condition <- "example_condition"
#' feature_names <- c("feature1", "feature2")
#' mode <- "isolated"
#' padjust_method <- "BH"
#' control_inputs_run_limma(data, meta, design, spline_params, condition, 
#'                          feature_names, mode, padjust_method)
#'
#' @seealso
#' \code{\link{check_design_formula}}, \code{\link{check_mode}}, 
#' \code{\link{check_spline_params}}
#' 
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
  
  check_design_formula(design, meta)
  
  check_condition(condition, meta)
  
  # Check that feature_names is a non-empty character vector
  if (!is.character(feature_names) || length(feature_names) == 0) {
    stop("feature_names must be a non-empty character vector.")
  }
  
  check_mode(mode)
  
  check_spline_params(spline_params, mode)
  
  check_padjust_method(padjust_method)
}


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
#' @examples
#' data <- matrix(runif(100), nrow = 10)
#' meta <- data.frame(Time = seq(1, 10), condition = rep(c("A", "B"), each = 5))
#' design <- "~ 1"
#' spline_params <- list(spline_type = c("n"), dof = list(3))
#' condition <- "condition"
#' compared_levels <- c("A", "B")
#' padjust_method <- "BH"
#' feature_names <- c("feature1", "feature2")
#' between_level(data, meta, design, spline_params, condition, compared_levels, 
#'               padjust_method, feature_names)
#'
#' @seealso
#' \code{\link{splines::bs}}, \code{\link{splines::ns}}, 
#' \code{\link{limma::lmFit}}, \code{\link{limma::eBayes}}, 
#' \code{\link{limma::topTable}}, \code{\link{modify_limma_top_table}}
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
  
  if (!is.null(spline_params$dof)) {
    args$df <- spline_params$dof[1]
  } else {
    args$knots <- spline_params$knots[[1]]
  }
  
  if (!is.null(spline_params$bknots)) {
    args$Boundary.knots <- spline_params$bknots[[1]]
  }
  
  
  if (spline_params$spline_type[1] == "b") {
    args$degree <- spline_params$degree[1]
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
#' @examples
#' top_table_and_fit <- list(
#'   top_table = data.frame(feature_index = 1:10, adj.P.Val = runif(10)), 
#'   fit = limma::lmFit(matrix(runif(100), nrow = 10), 
#'                      model.matrix(~1, data = data.frame(Time = seq(1, 10))))
#' )
#' feature_names <- c("feature1", "feature2")
#' process_top_table(top_table_and_fit, feature_names)
#'
#' @seealso
#' \code{\link{modify_limma_top_table}}, \code{\link{limma::lmFit}}
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

 
#' Process Level
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
#' @examples
#' data <- matrix(runif(100), nrow = 10)
#' meta <- data.frame(Time = seq(1, 10), condition = rep(c("A", "B"), each = 5))
#' design <- "~ 1"
#' spline_params <- list(spline_type = c("n"), dof = list(3))
#' level <- "A"
#' level_index <- 1
#' condition <- "condition"
#' feature_names <- c("feature1", "feature2")
#' padjust_method <- "BH"
#' mode <- "isolated"
#' process_level(level, level_index, spline_params, data, meta, design, 
#'               condition, feature_names, padjust_method, mode)
#'
#' @seealso
#' \code{\link{within_level}}, \code{\link{process_top_table}}
#' 
#' @importFrom stats relevel
#' 
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


#' Check Design Formula
#'
#' @description
#' Validates the design formula ensuring it is a valid character string, 
#' contains allowed characters, includes the intercept term 'X', and references 
#' columns present in the metadata.
#'
#' @param formula A character string representing the design formula.
#' @param meta A dataframe containing metadata.
#'
#' @return TRUE if the design formula is valid, otherwise an error is thrown.
#'
#' @examples
#' meta <- data.frame(Time = seq(1, 10), condition = rep(c("A", "B"), each = 5))
#' check_design_formula("~ Time + condition * X", meta)
#'
#' @seealso
#' \code{\link{stats::model.matrix}}
#' 
check_design_formula <- function(formula, 
                                 meta) {
  
  # Check if the formula is a valid character string
  if (!is.character(formula) || length(formula) != 1) {
    stop("The design formula must be a valid character string.")
  }
  
  # Ensure the formula contains allowed characters only
  allowed_chars <- "^[~ 1A-Za-z0-9_+*:()-]*$"
  if (!grepl(allowed_chars, formula)) {
    stop("The design formula contains invalid characters.")
  }
  
  # Ensure the formula contains the intercept term 'X'
  if (!grepl("\\bX\\b", formula)) {
    stop("The design formula must include the term 'X'.")
  }
  
  # Extract terms from the formula (removing interactions and functions)
  formula_terms <- unlist(strsplit(gsub("[~+*:()]", " ", formula), " "))
  formula_terms <- formula_terms[formula_terms != ""]
  
  # Remove '1' and 'X' from terms since they are not columns
  formula_terms <- setdiff(formula_terms, c("1", "X"))
  
  # Check if the terms are present in the dataframe
  missing_columns <- setdiff(formula_terms, names(meta))
  if (length(missing_columns) > 0) {
    stop(paste("The following columns are missing in the dataframe:", 
               paste(missing_columns, collapse = ", ")))
  }
  
  return(TRUE)
}

 
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
#' @examples
#' top_table <- data.frame(
#'   logFC = rnorm(10), 
#'   AveExpr = rnorm(10), 
#'   t = rnorm(10), 
#'   P.Value = runif(10), 
#'   adj.P.Val = runif(10)
#' )
#' feature_names <- paste0("Feature_", 1:10)
#' modify_limma_top_table(top_table, feature_names)
#'
#' @seealso
#' \code{\link{tidyr::as_tibble}}, \code{\link{dplyr::relocate}}, 
#' \code{\link{dplyr::mutate}}
#' 
#' @importFrom tidyr as_tibble
#' @importFrom dplyr relocate
#' @importFrom dplyr mutate
#' 
modify_limma_top_table <- function(top_table, 
                                   feature_names) {
  
  top_table <- tidyr::as_tibble(top_table, rownames = "feature_index") %>% 
    dplyr::relocate(feature_index, .after = last_col()) %>%
    dplyr::mutate(feature_index = as.integer(feature_index))
  
  sorted_feature_names <- feature_names[top_table$feature_index]
  top_table <- top_table %>% dplyr::mutate(feature_names = sorted_feature_names)
  
  top_table
}


#' Within Level Analysis
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
#' @examples
#' data <- matrix(runif(100), nrow = 10)
#' meta <- data.frame(Time = seq(1, 10), factor = rep(c("A", "B"), each = 5))
#' design <- "~ 1"
#' spline_params <- list(spline_type = c("n"), dof = list(3))
#' level <- "A"
#' level_index <- 1
#' padjust_method <- "BH"
#' within_level(data, meta, design, "factor", level, spline_params, 
#'              level_index, padjust_method)
#'
#' @seealso
#' \code{\link{splines::bs}}, \code{\link{splines::ns}}, 
#' \code{\link{limma::lmFit}}, \code{\link{limma::eBayes}}, 
#' \code{\link{limma::topTable}}
#' 
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
  
  if (!is.null(spline_params$dof)) {
    args$df <- spline_params$dof[level_index]
  } else {
    args$knots <- spline_params$knots[[level_index]]
  }
  
  if (!is.null(spline_params$bknots)) {
    args$Boundary.knots <- spline_params$bknots[[level_index]]
  }
  
  
  if (spline_params$spline_type[level_index] == "b") {
    args$degree <- spline_params$degree[level_index]
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
