#' limma_hyperparams_screen.R contains the exported package function 
#' limma_hyperparams_screen  and all the functions that make up the 
#' functionality of limma_hyperparams_screen. limma_hyperparams_screen runs the 
#' other package function, run_limma_splines, for a time series omics dataset,
#' for different hyperparameters. Such are for example degree of freedom of the
#' spline, type of spline, limma design formula, and different versions of the 
#' data (full data vs. outliers removed). This can result in several 
#' combinations, and it is tedious analyzing the combination in an organised 
#' and structured manner. Therefore, this function streamlines that process.


# Exported function: limma_hyperparams_screen() --------------------------------


#' Limma Hyperparameters Screening
#'
#' This function screens through various combinations of hyperparameters for 
#' limma analysis,
#' including designs, modes, and degrees of freedom. It validates inputs, 
#' generates results for all
#' combinations, and plots the outcomes. Finally, it may also be involved in 
#' generating an HTML report
#' as part of a larger analysis workflow.
#'
#' @param datas A list of data frames containing the datasets to be analyzed.
#' @param datas_descr A description object for the data.
#' @param metas A list of data frames containing metadata for each dataset in 
#' `datas`.
#' @param designs A character vector of design formulas for the limma analysis.
#' @param modes A character vector indicating the mode of analysis for each 
#' design.
#' @param condition A single character string specifying the condition.
#' @param spline_test_configs A configuration object for spline tests.
#' @param feature_names A character vector of feature names to be analyzed.
#' @param report_info An object containing report information.
#' @param report_dir A non-empty string specifying the report directory.
#' @param adj_pthresholds A numeric vector of p-value thresholds for 
#' significance determination.
#' @param meta_batch_column A character string specifying the meta batch column.
#' @param time_unit A character string specifying the time unit label for plots.
#' @param padjust_method A character string specifying the method for p-value 
#' adjustment.
#'
#' @return Returns a list of plots generated from the limma analysis results.
#'         Each element in the list corresponds to a different combination of 
#'         hyperparameters.
#' 
#' @importFrom here here
#'
#' @export
#'
limma_hyperparams_screen <- function(datas, 
                                     datas_descr,
                                     metas, 
                                     designs, 
                                     condition, 
                                     spline_test_configs,
                                     report_info,
                                     report_dir = here::here(),
                                     adj_pthresholds = c(0.05),
                                     meta_batch_column = NA,  # batch-effect
                                     time_unit = "m",    # For the plot labels
                                     padjust_method = "BH") {

  data_matrices <- list()
  for (data in datas) {
    matrix_and_feature_names <- process_data(data)
    data <- matrix_and_feature_names$data
    data_matrices <- c(data_matrices, list(data))
    feature_names <- matrix_and_feature_names$feature_names
  }
  datas <- data_matrices

  modes <- c()
  for (design in designs) {
    mode <- determine_analysis_mode(design,
                                    condition)
    modes <- c(modes, mode)
  }
  
  control_inputs_hyperpara_screen(datas = datas,
                                  datas_descr = datas_descr,
                                  metas = metas,
                                  designs = designs,
                                  modes = modes,
                                  condition = condition,
                                  spline_test_configs = spline_test_configs,
                                  feature_names = feature_names,
                                  report_info = report_info,
                                  report_dir = report_dir,
                                  adj_pthresholds = adj_pthresholds,
                                  meta_batch_column = meta_batch_column,
                                  time_unit = time_unit,
                                  padjust_method = padjust_method)
  
  top_tables_combos <- 
    get_limma_combos_results(datas = datas, 
                             metas = metas, 
                             designs = designs, 
                             modes = modes, 
                             condition = condition, 
                             spline_test_configs = spline_test_configs,
                             feature_names = feature_names, 
                             adj_pthresholds = adj_pthresholds, 
                             padjust_method = padjust_method)
  
  combo_pair_plots <- 
    plot_limma_combos_results(top_tables_combos = top_tables_combos, 
                              datas = datas, 
                              metas = metas,
                              spline_test_configs = spline_test_configs,
                              meta_batch_column = meta_batch_column,
                              designs,
                              time_unit = time_unit)
  
  timestamp <- format(Sys.time(), "%d_%m_%Y-%H_%M_%S")
  
  generate_reports(combo_pair_plots = combo_pair_plots, 
                   report_info = report_info,
                   report_dir = report_dir,
                   timestamp = timestamp)
  
  generate_reports_meta(datas_descr = datas_descr, 
                        designs = designs, 
                        modes = modes, 
                        spline_test_configs = spline_test_configs,
                        report_dir = report_dir,
                        timestamp = timestamp)
}


# Level 1 internal functions ---------------------------------------------------


#' Control Inputs for Hyperparameter Screening
#'
#' @description
#' Validates and checks the inputs for hyperparameter screening. Ensures all 
#' provided data structures meet the required format and criteria.
#'
#' @param datas A list of matrices.
#' @param datas_descr A description object for the data.
#' @param metas A list of metadata corresponding to the data matrices.
#' @param designs A list of design matrices.
#' @param modes A character vector containing 'isolated' or 'integrated'.
#' @param condition A single character string specifying the condition.
#' @param spline_test_configs A configuration object for spline tests.
#' @param feature_names A character vector of feature names.
#' @param report_info An object containing report information.
#' @param report_dir A non-empty string specifying the report directory.
#' @param adj_pthresholds A numeric vector with elements > 0 and < 1.
#' @param meta_batch_column A character string specifying the meta batch column.
#' @param time_unit A character string specifying the time unit label for plots.
#' @param padjust_method A single character string specifying the p-adjustment 
#' method.
#'
#' @return No return value, called for side effects.
#'
control_inputs_hyperpara_screen <- function(datas, 
                                            datas_descr,
                                            metas, 
                                            designs, 
                                            modes, 
                                            condition, 
                                            spline_test_configs,
                                            feature_names, 
                                            report_info,
                                            report_dir,
                                            adj_pthresholds, 
                                            meta_batch_column,
                                            time_unit,
                                            padjust_method) {
  
  if(length(datas) != length(metas) || length(designs) != length(modes)) {
    stop(paste0("datas and metas, as well as designs and modes must have the ", 
         "same length."),
         call. = FALSE)
  }
  
  for (i in seq_along(datas)) {
    check_data_and_meta(data = datas[[i]], 
                        meta = metas[[i]], 
                        condition = condition, 
                        meta_batch_column = meta_batch_column,
                        data_meta_index = i)
  }

  check_datas_descr(datas_descr)
  
  for (design in designs) {
    for (j in seq_along(metas)) {
      meta <- metas[[j]]
      check_design_formula(formula = design, 
                           meta = meta,
                           meta_index = j)
    }
  }
    
  if (!is.character(modes) || !all(modes %in% c("isolated", "integrated"))) {
    stop("'modes' must be a character vector containing only 'isolated' or 
         'integrated'.",
         call. = FALSE)
  }
  
  check_spline_test_configs(spline_test_configs, metas)
  
  if (!is.character(feature_names)) {
    stop("'feature_names' must be a character vector.",
         call. = FALSE)
  }
  
  check_report_info(report_info)
  
  check_and_create_report_dir(report_dir)
  
  check_adj_pthresholds(adj_pthresholds)
  
  check_time_unit(time_unit)

  check_padjust_method(padjust_method)
}


#' Generate LIMMA Combination Results
#'
#' @description
#' Computes results for various combinations of data, design matrices, and 
#' spline configurations using the LIMMA method.
#'
#' @param datas A list of matrices.
#' @param metas A list of metadata corresponding to the data matrices.
#' @param designs A list of design matrices.
#' @param modes A character vector containing 'isolated' or 'integrated'.
#' @param condition A single character string specifying the condition.
#' @param spline_test_configs A configuration object for spline tests.
#' @param feature_names A character vector of feature names.
#' @param adj_pthresholds A numeric vector with elements > 0 and < 1.
#' @param padjust_method A single character string specifying the p-adjustment 
#' method.
#'
#' @return A list of results for each combination of data, design, and spline 
#' configuration.
#'                          
#' @importFrom magrittr %>%
#' @importFrom tidyr expand_grid
#' @importFrom dplyr mutate
#' @importFrom purrr pmap
#' @importFrom purrr set_names
#' @importFrom rlang sym
#'                          
get_limma_combos_results <- function(datas, 
                                     metas, 
                                     designs, 
                                     modes, 
                                     condition, 
                                     spline_test_configs,
                                     feature_names, 
                                     adj_pthresholds, 
                                     padjust_method) {
  
  combos <- tidyr::expand_grid(
    data_index = seq_along(datas),
    design_index = seq_along(designs),
    spline_config_index = seq_along(spline_test_configs$spline_type),
    pthreshold = adj_pthresholds
  ) %>% 
    dplyr::mutate(id = paste0("Data_", !!rlang::sym("data_index"), 
                              "_Design_", !!rlang::sym("design_index"), 
                              "_SConfig_", !!rlang::sym("spline_config_index"), 
                              "_PThresh_", !!rlang::sym("pthreshold")))
  
    purrr::pmap(combos, 
                process_combo,
                datas = datas,
                metas = metas,
                designs = designs,
                modes = modes,
                condition = condition, 
                spline_test_configs = spline_test_configs,
                feature_names = feature_names, 
                padjust_method = padjust_method) %>% 
    purrr::set_names(combos$id)
}


#' Plot limma Combination Results
#'
#' @description
#' Generates plots for pairwise comparisons of hyperparameter combinations 
#' using limma results.
#'
#' @param top_tables_combos A list of top tables for each combination.
#' @param datas A list of matrices.
#' @param metas A list of metadata corresponding to the data matrices.
#' @param spline_test_configs A configuration object for spline tests.
#' @param meta_batch_column A character string specifying the meta batch column.
#' @param time_unit A single character, such as s, m, h, or d, specifying the
#' time_unit that should be used for the plots (s = seconds, m = minutes,
#' h = hours, d = days). This single character will be converted to a string 
#' that is a little bit more verbose, such as sec in square brackets for s.
#'
#' @return A list of results including hit comparison plots and composite 
#' spline plots for each pair of combinations.
#'                           
#' @importFrom stringr str_extract
#' @importFrom progress progress_bar
#' @importFrom purrr set_names
#' @importFrom purrr map
#'                           
plot_limma_combos_results <- function(top_tables_combos,
                                      datas,
                                      metas,
                                      spline_test_configs,
                                      meta_batch_column,
                                      designs,
                                      time_unit = time_unit) {
  
  names_extracted <- stringr::str_extract(names(top_tables_combos),
                                          "Data_\\d+_Design_\\d+")
  
  combos_separated <- lapply(unique(names_extracted), function(id) {
    top_tables_combos[names_extracted == id]
  })
  
  names(combos_separated) <- unique(names_extracted)
  
  combos <- names(combos_separated)
  combo_pairs <- combn(combos, 2, simplify = FALSE)
  
  print("Generating the plots for all pairwise hyperparams-combo comparisons")
  progress_ticks <- length(combo_pairs)
  pb <- progress::progress_bar$new(total = progress_ticks, 
                                   format = "[:bar] :percent")
  pb$tick(0)
  
  time_unit_labels <- c("s" = "[sec]", "m" = "[min]", "h" = "[hours]", 
                        "d" = "[days]")
  time_unit_label <- time_unit_labels[time_unit]
  
  if (!is.na(meta_batch_column)) {
    
    # Takes the shortcut approach without the specific design_matrix
    datas <- remove_batch_effect(datas = datas, 
                                 metas = metas,
                                 condition = condition,
                                 meta_batch_column = meta_batch_column)
  }
  
  
  combo_pair_results <- purrr::set_names(
    purrr::map(combo_pairs, function(pair) {
      
      combo_pair <- combos_separated[pair]
      
      hitcomp <- gen_hitcomp_plots(combo_pair)
      
      composites <- purrr::map(combo_pair, function(combo) {
        composite <- gen_composite_spline_plots(combo,
                                                datas,
                                                metas,
                                                spline_test_configs,
                                                time_unit_label)
      })
      pb$tick()
      list(hitcomp = hitcomp, composites = composites)
    }
    ), purrr::map(combo_pairs, function(pair) paste(pair[1], "vs", pair[2],
                                             sep = "_"))
  )
}


#' Generate Reports
#'
#' @description
#' Builds HTML reports for all pairwise hyperparameter combination comparisons.
#'
#' @param combo_pair_plots A list of plots for each pair of combinations.
#' @param report_info An object containing report information.
#' @param report_dir A non-empty string specifying the report directory.
#' @param timestamp A timestamp to include in the reports.
#'
#' @return No return value, called for side effects.
#'                  
#' @importFrom progress progress_bar
#' @importFrom purrr imap
#'                  
generate_reports <- function(combo_pair_plots, 
                             report_info,
                             report_dir,
                             timestamp) {
  
  print("Building .html reports for all pairwise hyperparams-combo comparisons")
  progress_ticks <- length(combo_pair_plots)
  pb <- progress::progress_bar$new(total = progress_ticks, 
                                   format = "[:bar] :percent")
  
  result <- purrr::imap(combo_pair_plots, ~{
    process_combo_pair(.x, .y, report_info, report_dir, timestamp)
    pb$tick()  
  })
}


#' Generate Reports Metadata
#'
#' @description
#' Generates a metadata table for the LIMMA hyperparameter screen reports and 
#' saves it as an HTML file with custom styling.
#'
#' @param datas_descr A description object for the data.
#' @param designs A list of design matrices.
#' @param modes A character vector containing 'isolated' or 'integrated'.
#' @param spline_test_configs A configuration object for spline tests.
#' @param report_dir A non-empty string specifying the report directory.
#' @param timestamp A timestamp to include in the report filename.
#'
#' @return No return value, called for side effects.
#'                       
#' @importFrom here here
#' @importFrom knitr kable
#' @importFrom kableExtra kable_styling
#'                       
generate_reports_meta <- function(datas_descr, 
                                  designs, 
                                  modes, 
                                  spline_test_configs,
                                  report_dir,
                                  timestamp) {
  
  formatted_spline_configs <- flatten_spline_configs(spline_test_configs)
  
  # Combine the hyperparameters and their descriptions into two vectors
  hyperparameters <- c(paste0("Data_", seq_along(datas_descr)), 
                       paste0("Design_", seq_along(designs)),
                       paste0("SConfig_", seq_along(formatted_spline_configs)))
  descriptions <- c(datas_descr, paste(designs, "(mode:", modes, ")"), 
                    unlist(formatted_spline_configs))
  
  table_df <- data.frame(hyperparameter = hyperparameters, 
                         description = descriptions, 
                         stringsAsFactors = FALSE)
  
  filename <- sprintf("hyperparams_screen_meta_table_%s.html", timestamp)
  file_path <- here::here(report_dir, filename)
  
  custom_css <- "
  <style>
  table {
    font-size: 32px; /* Significantly larger font size */
    margin-left: auto; /* Center table horizontally */
    margin-right: auto;
  }
  th, td {
    border: 1px solid #cccccc;
    padding: 12px; /* Increased padding for more space between table cells */
    text-align: left;
  }
  th {
    background-color: #f2f2f2;
  }
  </style>
  "
  
  # Generate the HTML table with custom CSS for larger font size and save it
  html_table <- 
    paste0(custom_css, 
           knitr::kable(table_df, format = "html", escape = FALSE) %>% 
             kableExtra::kable_styling(
               bootstrap_options = c("striped", "hover")))
  
  writeLines(html_table, con = file_path)
  
  cat("Meta table for the limma hyperparameter screen reports saved to:",
      file_path, "\n")
}



# Level 2 internal functions ---------------------------------------------------


#' Check Data Descriptions
#'
#' @description
#' Validates that the data descriptions are character vectors with each element 
#' not exceeding 80 characters in length.
#'
#' @param datas_descr A character vector of data descriptions.
#'
#' @return No return value, called for side effects.
#'
#' @seealso
#' \code{\link{stop}} for error handling.
#' 
check_datas_descr <- function(datas_descr) {
  
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
}


#' Check Spline Test Configurations
#'
#' @description
#' Validates the spline test configurations by checking required columns, 
#' spline type column, spline type parameters, and maximum and minimum degrees 
#' of freedom against the metadata.
#'
#' @param spline_test_configs A configuration object for spline tests.
#' @param metas A list of metadata corresponding to the data matrices.
#'
#' @return No return value, called for side effects.
#'
#' @seealso
#' \code{\link{check_colums_spline_test_configs}}, 
#' \code{\link{check_spline_type_column}}, 
#' \code{\link{check_spline_type_params}}, 
#' \code{\link{check_max_and_min_dof}}
#' 
check_spline_test_configs <- function(spline_test_configs, 
                                      metas) {
  
  check_colums_spline_test_configs(spline_test_configs)
  
  check_spline_type_column(spline_test_configs)
  
  check_spline_type_params(spline_test_configs)
  
  check_max_and_min_dof(spline_test_configs, metas)
}


#' Process Combination
#'
#' @description
#' Processes a single combination of data, design, spline configuration, and 
#' p-threshold to generate LIMMA spline results.
#'
#' @param data_index Index of the data in the datas list.
#' @param design_index Index of the design in the designs list.
#' @param spline_config_index Index of the spline configuration in the 
#' spline_test_configs list.
#' @param pthreshold The p-value threshold for significance.
#' @param datas A list of data matrices
#' @param metas A list of metadata corresponding to the data matrices.
#' @param designs A list of design matrices.
#' @param modes A character vector containing 'isolated' or 'integrated'.
#' @param condition A single character string specifying the condition.
#' @param spline_test_configs A configuration object for spline tests.
#' @param feature_names A character vector of feature names.
#' @param padjust_method A single character string specifying the p-adjustment 
#' method.
#' @param ... Additional arguments.
#'
#' @return A list of top tables from the LIMMA spline analysis.
#'
#' @seealso
#' \code{\link{create_spline_params}}, 
#' \code{\link{run_limma_splines}}
#' 
process_combo <- function(data_index, 
                          design_index, 
                          spline_config_index, 
                          pthreshold, 
                          datas,
                          metas,
                          designs,
                          modes,
                          condition,
                          spline_test_configs,
                          feature_names,
                          padjust_method,
                          ...) {
  
  data <- datas[[data_index]]
  meta <- metas[[data_index]]
  design <- designs[[design_index]]
  mode <- modes[[design_index]]
  
  spline_params <- create_spline_params(spline_test_configs = 
                                          spline_test_configs, 
                                        index = spline_config_index, 
                                        meta = meta, 
                                        condition = condition, 
                                        mode = mode)
  
  # Because either DoF or knots are specified, and only optionally bknots
  # If they are not specified, their value is NA.
  spline_params <- Filter(is_not_na, spline_params)
  
  data <- as.data.frame(data)
  rownames(data) <- feature_names
  
  # suppressMessages will not affect warnings and error messages!
  result <- suppressMessages(run_limma_splines(data = data, 
                                               meta = meta, 
                                               design = design, 
                                               spline_params = spline_params, 
                                               condition = condition,
                                               padjust_method = padjust_method))
  
  
  result$top_tables
}


#' Remove Batch Effect
#'
#' @description
#' Removes batch effects from the data matrices using the specified batch 
#' column in the metadata.
#'
#' @param datas A list of matrices.
#' @param metas A list of metadata corresponding to the data matrices.
#' @param meta_batch_column A character string specifying the meta batch column.
#'
#' @return A list of matrices with batch effects removed where applicable.
#'
#' @seealso
#' \link[limma]{removeBatchEffect}
#' 
#' @importFrom limma removeBatchEffect
#'
remove_batch_effect <- function(datas, 
                                metas,
                                meta_batch_column,
                                condition) {
  
  results <- list()
  for (i in seq_along(datas)) {
    data <- datas[[i]]
    meta <- metas[[i]]
    
    # This is the shortcut approach. It would be better to remove the batch
    # effect using the design_matrix, but it is challenging to program this here
    data_corrected <- removeBatchEffect(x = data, 
                                        batch = meta[[meta_batch_column]],
                                        group = meta[[condition]])
    results[[i]] <- data_corrected
  }
  return(results)
}


#' Create New Hit Comparison Object
#'
#' @description
#' Creates a new hit comparison object with specified condition names.
#'
#' @param cond1name A character string for the first condition name 
#' (max length 25).
#' @param cond2name A character string for the second condition name 
#' (max length 25).
#'
#' @return An object of class "hitcomp" containing empty data lists 
#' and condition names.
#' 
hc_new <- function(cond1name = "Condition 1", 
                   cond2name = "Condition 2") {
  
  if (!is.character(cond1name) || nchar(cond1name) > 25) 
    stop("cond1name max len = 25")
  if (!is.character(cond2name) || nchar(cond2name) > 25) 
    stop("cond2name max len = 25")
  
  res <- list(
    data = list(
      list(),
      list()
    ),
    condition_names = c(cond1name, cond2name)
  )
  class(res) <- "hitcomp"
  return(res)
}


#' Add Data to Hit Comparison Object
#'
#' @description
#' Adds a new entry to the hit comparison object for a specified condition.
#'
#' @param hc_obj An object of class "hitcomp".
#' @param top_table A dataframe containing the top table data.
#' @param params_id A character string identifying the parameters 
#' (max length 70).
#' @param condition An integer (1 or 2) specifying the condition to which the 
#' data belongs.
#' @param threshold A numeric value specifying the adjusted p-value threshold.
#'
#' @return The updated hit comparison object.
#'              
hc_add <- function(hc_obj, 
                   top_table, 
                   params_id, condition = 1, 
                   threshold = 0.05){
  
  # Validate input
  if (!is.data.frame(top_table)) stop("top_table must be a dataframe.")
  if (!is.character(params_id) || nchar(params_id) > 70) {
    stop(paste("max len = 70. params_id '", params_id, "' is too long.", 
               sep = ""))
  }
  if (!(condition %in% c(1, 2))) stop("condition must be either 1 or 2.")
  if (!is.numeric(threshold)) stop("threshold must be numeric.")
  
  # Create the list to append
  new_entry <- list(
    DataFrame = top_table,
    Parameters = list(
      params_id = params_id,
      adj_p_value_threshold = threshold
    )
  )
  
  # Append to the appropriate global list in the package environment
  hc_obj$data[[condition]] <- append(hc_obj$data[[condition]], list(new_entry))
  return(hc_obj)
}


#' Generate Venn Heatmap
#'
#' @description
#' Creates a Venn heatmap to visualize the overlap of hits between two 
#' conditions stored in a hit comparison object.
#'
#' @param hc_obj An object of class "hitcomp" containing hit data for two 
#' conditions.
#'
#' @return A list containing the Venn heatmap plot and the number of hits.
#'
#' @seealso
#' \code{\link{store_hits}}, \link[pheatmap]{pheatmap}
#' 
#' @importFrom tidyr expand_grid unnest_longer replace_na pivot_wider
#' @importFrom tibble enframe column_to_rownames
#' @importFrom dplyr mutate select left_join
#' @importFrom pheatmap pheatmap
#' @importFrom purrr flatten_chr
#' 
hc_vennheatmap <- function(hc_obj) {
  
  hits_1 <- store_hits(hc_obj$data[[1]])
  hits_2 <- store_hits(hc_obj$data[[2]])
  
  df <- tidyr::expand_grid(
    features = union(
      flatten_chr(hits_1),
      flatten_chr(hits_2)
    ),
    params = union(
      names(hits_1),
      names(hits_2)
    )
  )
  
  df_1 <-
    hits_1 %>% 
    tibble::enframe("params", "features") %>% 
    tidyr::unnest_longer(!!rlang::sym("features")) %>%
    dplyr::mutate(x1 = 1)
  
  df_2 <-
    hits_2 %>% 
    tibble::enframe("params", "features") %>% 
    tidyr::unnest_longer(!!rlang::sym("features")) %>%
    dplyr::mutate(x2 = 2)
  
  venn_matrix <- 
    df %>% 
    dplyr::left_join(df_1, by = c("features", "params")) %>% 
    dplyr::left_join(df_2, by = c("features", "params")) %>% 
    tidyr::replace_na(list(x1 = 0, x2 = 0)) %>% 
    dplyr::mutate(x = !!rlang::sym("x1") + !!rlang::sym("x2")) %>% 
    dplyr::select(!c(!!rlang::sym("x1"), !!rlang::sym("x2"))) %>%
    tidyr::pivot_wider(names_from = !!rlang::sym("params"), 
                       values_from = !!rlang::sym("x")) %>% 
    tibble::column_to_rownames("features") %>% 
    as.matrix()
  
  venn_matrix <- venn_matrix[, order(colnames(venn_matrix))]
  
  color_palette <- c("white", "blue", "yellow", "green")
  breaks <- c(-0.5, 0.5, 1.5, 2.5, 3.5)
  
  plot_title <- sprintf("0 -> none, 1 -> %s, 2 -> %s, 3 -> both", 
                        hc_obj$condition_names[[1]],
                        hc_obj$condition_names[[2]])
  
  vennheatmap_plot <- pheatmap::pheatmap(venn_matrix, color = color_palette,
                                         breaks = breaks,
                                         cluster_cols = FALSE,
                                         cluster_rows = TRUE,
                                         show_rownames = TRUE,
                                         show_colnames = TRUE,
                                         border_color = NA,
                                         main = plot_title,
                                         silent = TRUE)
  
  return(list(vennheatmap = vennheatmap_plot, 
              nrhits = nrow(venn_matrix)))
}


#' Generate Barplot for Hit Comparison Object
#'
#' @description
#' Creates a barplot to visualize the number of significant features for each 
#' parameter set in the hit comparison object.
#'
#' @param hc_obj An object of class "hitcomp" containing hit data for two 
#' conditions.
#'
#' @return A ggplot2 object representing the barplot.
#'
#' @seealso
#' \code{\link{store_hits}}, \code{\link{ggplot2}}
#' 
#' @importFrom ggplot2 geom_col geom_text facet_wrap scale_y_continuous aes
#' @importFrom ggplot2 theme_minimal element_text element_blank expansion xlab
#' @importFrom ggplot2 theme
#' @importFrom purrr map_int set_names 
#' @importFrom tibble enframe
#' @importFrom dplyr bind_rows vars
#' 
hc_barplot <- function(hc_obj) {
  
  plot_data <- 
    list(
      store_hits(hc_obj$data[[1]]) %>% 
        purrr::map_int(length) %>% 
        tibble::enframe("params", "n_hits"),
      store_hits(hc_obj$data[[2]]) %>% 
        purrr::map_int(length) %>% 
        tibble::enframe("params", "n_hits")
    ) %>% 
    purrr::set_names(hc_obj$condition_names) %>% 
    dplyr::bind_rows(.id = "condition")
  
  ggplot2::ggplot(plot_data, aes(x = !!rlang::sym("params"), 
                                 y = !!rlang::sym("n_hits"))) +
    geom_col() +
    geom_text(
      aes(label = !!rlang::sym("n_hits")),
      vjust = -0.5,
      
    ) +
    xlab("Parameters") +
    scale_y_continuous(
      "number of significant features",
      expand = expansion(mult = c(0, .2))
    ) +
    facet_wrap(vars(!!rlang::sym("condition")), ncol = 1) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      panel.grid.major.x = element_blank()
    ) +
    NULL
}


#' Generate Hit Comparison Plots
#'
#' @description
#' Generates Venn heatmap and barplot for a given combination pair of top 
#' tables.
#'
#' @param combo_pair A list containing two combinations of top tables.
#'
#' @return A list containing the Venn heatmap plot, the number of hits divided 
#' by 16, the barplot, and a length indicator for the barplot.
#'
#' @seealso
#' \code{\link{hc_new}}, \code{\link{hc_add}}, \code{\link{hc_vennheatmap}}, 
#' \code{\link{hc_barplot}}
#' 
gen_hitcomp_plots <- function(combo_pair) {
  
  combo_pair_combined <- c(combo_pair[[1]], combo_pair[[2]])
  
  hitcomp <- hc_new(names(combo_pair)[1], names(combo_pair)[2])
  
  combo_names <- sapply(names(combo_pair_combined), function(name) {
    # Using sub to extract everything from 'DoF' onwards
    match <- sub(".*(SConfig.*)$", "\\1", name)
    return(match)
  })
  
  condition_nr <- 1L
  for (i in 1:length(combo_pair_combined)) {
    combo_top_tables <- combo_pair_combined[[i]]
    combo_name <- combo_names[i]
    # Extract the part where DoF and p-value threshold are written.
    pthreshold <- as.numeric(sub(".*_([^_]+)$", "\\1", combo_name))
    
    if (i == length(combo_pair_combined)/2 + 1) {
      condition_nr <- 2L
    }
    
    for (top_table_name in names(combo_top_tables)) {
      top_table <- combo_top_tables[[top_table_name]]
      id <- paste(top_table_name, combo_name, sep = "_")
      hitcomp <- hc_add(hitcomp, top_table, id, condition_nr, pthreshold)
    }
  }
  
  result <- hc_vennheatmap(hitcomp)
  barplot <- hc_barplot(hitcomp)
  
  list(vennheatmap = result$vennheatmap, 
       vennheatmap_len = c(as.integer(result$nrhits/16)),
       barplot = barplot,
       barplot_len = c(2L)
  )
}

 
#' Generate Composite Spline Plots
#'
#' @description
#' Creates composite spline plots for significant and non-significant features 
#' across multiple levels within a condition.
#' One half of one condition comparison HTML 
#' (composite spline plots for one 'condition' inside one condition comparison)
#'
#' @param internal_combos A list containing combinations of top tables.
#' @param datas A list of matrices.
#' @param metas A list of metadata corresponding to the data matrices.
#' @param spline_test_configs A configuration object for spline tests.
#' @param time_unit_label A character string specifying the time unit label 
#' for plots.
#'
#' @return A list containing the composite spline plots and their lengths.
#'
#' @seealso
#' \code{\link{plot_composite_splines}}
#' 
#' @importFrom utils tail
#' 
gen_composite_spline_plots <- function(internal_combos, 
                                       datas, 
                                       metas,
                                       spline_test_configs,
                                       time_unit_label) {
  
  plots <- list()
  plots_len <- integer(0)
  
  # all the combos of DoF and adj. p-value threshold for one condition
  for (combo_name in names(internal_combos)) {
    top_tables_levels <- internal_combos[[combo_name]]
    
    data <- datas[[as.integer(strsplit(combo_name, "_")[[1]][2])]]
    meta <- metas[[as.integer(strsplit(combo_name, "_")[[1]][2])]]
    
    pthresh <- as.numeric(utils::tail(strsplit(combo_name, "_")[[1]], 1)) 
    
    # for one given combo of DoF and adj. p-value threshold, within one 
    # condition, there are multiple levels (for example exp and stat)
    for (top_table_name in names(top_tables_levels)) {
      top_table <- top_tables_levels[[top_table_name]]
      
      parts <- strsplit(top_table_name, "_")[[1]]
      condition <- parts[1]
      level <- parts[2]
      meta_level <- meta[meta[[condition]] == level, ]
      data_level <- data[, which(meta[[condition]] == level)]
      
      # Show 6 significant and 6 non significant splines, each within a 
      # composite plot (the 6 individual plots combined with patchwork)
      for(type in c('significant', 'not_significant')) {
        
        if (type == "significant") {
          filtered_rows <- top_table[top_table$adj.P.Val < pthresh, ]
          selected_rows <- if(nrow(filtered_rows) > 6) {
            filtered_rows[sample(nrow(filtered_rows), 6), ]
          } else {
            filtered_rows
          }
          indices <- as.integer(selected_rows$feature_index)
        } else if (type == "not_significant") {
          filtered_rows <- top_table[top_table$adj.P.Val >= pthresh, ]
          selected_rows <- if(nrow(filtered_rows) > 6) {
            filtered_rows[sample(nrow(filtered_rows), 6), ]
          } else {
            filtered_rows
          }
          indices <- as.integer(selected_rows$feature_index)
        }
        
        # One composite spline plot for each unique combo between DoF and 
        # adj. p-value threshold (for one level within one condition)
        # This fun just generates a single composite plot
        result <- plot_composite_splines(data_level, 
                                         meta_level,
                                         spline_test_configs,
                                         top_table, 
                                         combo_name, 
                                         indices,
                                         type,
                                         time_unit_label)
        
        if (is.list(result)) {
          plot_name <- paste(combo_name, top_table_name, type, sep = "_")
          plots[[plot_name]] <- result$composite_plot
          plots_len <- c(plots_len, result$composite_plot_len)
        }
      }
    }
  }
  
  list(composite_plots = plots, 
       composite_plots_len = plots_len)
}


#' Process Combination Pair
#'
#' @description
#' Processes a combination pair to generate plots and compile them into an 
#' HTML report.
#'
#' @param combo_pair A list containing hit comparison and composite spline 
#' plots.
#' @param combo_pair_name A character string for naming the combination pair.
#' @param report_info An object containing report information.
#' @param report_dir A non-empty string specifying the report directory.
#' @param timestamp A timestamp to include in the report filename.
#'
#' @return No return value, called for side effects.
#'
#' @seealso
#' \code{\link{generate_report_html}}
#' 
process_combo_pair <- function(combo_pair, 
                               combo_pair_name, 
                               report_info,
                               report_dir,
                               timestamp) {
  
  plots <- list()
  plots_len <- integer(0)
  
  hitcomp <- combo_pair$hitcomp
  
  plots[[1]] <- hitcomp$vennheatmap
  plots[[2]] <- hitcomp$vennheatmap
  plots[[3]] <- hitcomp$barplot
  
  plots_len <- c(plots_len, 
                 2, 
                 hitcomp$vennheatmap_len, 
                 hitcomp$barplot_len)
  
  composites <- combo_pair$composites
  
  for (composite in composites) {
    for (plot in composite$composite_plots) {
      plots[[length(plots) + 1]] <- plot
    }
    
    for (len in composite$composite_plots_len) {
      plots_len <- c(plots_len, len)
    }
  }
  
  # Function is in splinetime_general_fun.R
  generate_report_html(plots = plots, 
                       plots_sizes = plots_len, 
                       report_info = report_info,
                       report_type = "limma_hyperparams_screen",
                       filename = combo_pair_name,
                       timestamp = timestamp,
                       report_dir = report_dir, )
}


#' Create Spline Parameters
#'
#' @description
#' Generates spline parameters based on the configuration, metadata, condition, 
#' and mode.
#'
#' @param spline_test_configs A configuration object for spline tests.
#' @param index Index of the spline configuration to process.
#' @param meta A dataframe containing metadata.
#' @param condition A character string specifying the condition.
#' @param mode A character string specifying the mode.
#'
#' @return A list of processed spline parameters.
#'
#' @seealso
#' \code{\link{process_config_column}}
#' 
create_spline_params <- function(spline_test_configs, 
                                 index, 
                                 meta, 
                                 condition, 
                                 mode) {
  
  num_levels <- length(unique(meta[[condition]]))  
  
  result <- lapply(spline_test_configs, 
                   process_config_column, 
                   index, 
                   num_levels, 
                   mode)
}


#' Flatten Spline Configurations
#'
#' @description
#' Flattens and formats spline configurations into a list of formatted strings.
#'
#' @param spline_configs A list of spline configuration objects.
#'
#' @return A list of formatted strings representing each spline configuration.
#' 
flatten_spline_configs <- function(spline_configs) {
  
  formatted_layers <- list()
  
  names_of_spline_configs <- names(spline_configs)
  for (i in 1:length(spline_configs$spline_type)) {
    ith_elements <- lapply(spline_configs, function(x) x[[i]])
    
    formatted_strings <- list()
    
    Map(function(name, element) {
      formatted_strings <<- 
        c(formatted_strings, paste0(name, " = ", 
                                    paste(element, collapse = ", ")))
      
    }, names_of_spline_configs, ith_elements)
    
    final_string <- paste(formatted_strings, collapse = ", ")
    formatted_layers[[length(formatted_layers) + 1]] <- final_string 
    
  }
  formatted_layers
}


# Level 3 internal functions  --------------------------------------------------


#' Check Columns in Spline Test Configurations
#'
#' @description
#' Validates that the spline test configurations contain the required columns 
#' in the correct order.
#'
#' @param spline_test_configs A dataframe containing spline test configurations.
#'
#' @return No return value, called for side effects.
#' 
check_colums_spline_test_configs <- function(spline_test_configs) {
  
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
}


#' Check Spline Type Column
#'
#' @description
#' Validates that the 'spline_type' column in the spline test configurations 
#' contains only 'n' or 'b'.
#'
#' @param spline_test_configs A dataframe containing spline test configurations.
#'
#' @return No return value, called for side effects.
#' 
check_spline_type_column <- function(spline_test_configs) {
  
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
}


#' Check Spline Type Parameters
#'
#' @description
#' Validates the parameters for each row in the spline test configurations 
#' based on the spline type.
#'
#' @param spline_test_configs A dataframe containing spline test configurations.
#'
#' @return TRUE if all checks pass, otherwise an error is thrown.
#' 
check_spline_type_params <- function(spline_test_configs) {
  
  for (i in seq_len(nrow(spline_test_configs))) {
    row <- spline_test_configs[i,]
    switch(as.character(row$spline_type),
           "n" = {
             if (!is.na(row$degree)) stop("degree must be NA for spline_type n")
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
}


#' Check Maximum and Minimum Degrees of Freedom
#'
#' @description
#' Validates the degrees of freedom (DoF) for each row in the spline test 
#' configurations based on the metadata.
#'
#' @param spline_test_configs A dataframe containing spline test configurations.
#' @param metas A list of metadata corresponding to the data matrices.
#'
#' @return No return value, called for side effects.
#' 
#' @importFrom stats setNames
#' 
check_max_and_min_dof <- function(spline_test_configs, 
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
}


#' Check if Not All Values are NA
#'
#' @description
#' Determines if a given atomic vector contains at least one non-NA value.
#'
#' @param x An atomic vector or any other object.
#'
#' @return TRUE if the vector contains at least one non-NA value or if the 
#' object is not atomic; FALSE otherwise.
#' 
is_not_na <- function(x) {
  
  if (is.atomic(x)) {
    return(!all(is.na(x)))
  } else {
    return(TRUE)
  }
}


#' Process Configuration Column
#'
#' @description
#' Processes a configuration column based on the given mode and number of 
#' levels.
#'
#' @param config_column A configuration column from the spline test 
#' configurations.
#' @param index Index of the configuration to process.
#' @param num_levels Number of unique levels in the metadata condition.
#' @param mode A character string specifying the mode 
#' ('integrated' or 'isolated').
#'
#' @return A vector or list with the processed configuration values.
#' 
process_config_column <- function(config_column, 
                                  index, 
                                  num_levels, 
                                  mode) {
  
  if (mode == "integrated") {
    if (is.list(config_column)) {
      config_column[[index]]
    } else {
      config_column[index]
    }
  } else {     # mode = 'isolated'
    if (is.list(config_column)) {
      rep(config_column[[index]], num_levels)
    } else {
      rep(config_column[index], num_levels)
    }
  } 
}


#' Store Hits
#'
#' @description
#' Stores the feature indices for significant hits based on the adjusted p-value 
#' threshold for each condition.
#'
#' @param condition A list containing dataframes and parameters for each 
#' condition.
#'
#' @return A list where each element is a vector of feature indices that meet 
#' the significance threshold.
#' 
#' @importFrom dplyr pull
#' 
store_hits <- function(condition) {
  
  hits_cond <- list()
  
  for (item in condition) {
    df <- item$DataFrame
    adj_p_value_treshold <- item$Parameters$adj_p_value_threshold
    key_name <- item$Parameters$params_id

    hits_cond[[key_name]] <- df %>%
      filter(df[["adj.P.Val"]] < adj_p_value_treshold) %>%
      pull(!!sym("feature_index")) %>%
      as.character()
  }
  
  return(hits_cond)
}


#' Plot Composite Splines
#'
#' @description
#' Generates composite spline plots for significant and non-significant 
#' features based on the specified indices.
#'
#' @param data A matrix of data values.
#' @param meta A dataframe containing metadata.
#' @param spline_test_configs A configuration object for spline tests.
#' @param top_table A dataframe containing the top table results.
#' @param top_table_name A character string specifying the name of the 
#' top table.
#' @param indices A vector of indices specifying which features to plot.
#' @param type A character string specifying the type of features ('significant' 
#' or 'not_significant').
#' @param time_unit_label A string shown in the plots as the unit for the time,
#' such as min or hours.
#'
#' @return A list containing the composite plot and its length if plots are 
#' generated, FALSE otherwise.
#'
#' @seealso
#' \link[splines]{bs}, \link[splines]{ns}, \link[ggplot2]{ggplot2}, 
#' \link[patchwork]{wrap_plots}
#' 
#' @importFrom splines ns
#' @importFrom ggplot2 ggplot geom_point geom_line theme_minimal 
#' @importFrom ggplot2 scale_x_continuous labs annotate theme
#' @importFrom patchwork wrap_plots plot_annotation
#' 
plot_composite_splines <- function(data, 
                                   meta, 
                                   spline_test_configs,
                                   top_table, 
                                   top_table_name, 
                                   indices,
                                   type,
                                   time_unit_label) {
  
  plot_list <- list()
  config_index <- 
    as.integer(sub(".*SConfig_([0-9]+)_.*", "\\1", top_table_name)) 

  smooth_timepoints <- seq(meta$Time[1], 
                           meta$Time[length(meta$Time)], 
                           length.out = 100)
  
  args <- list(x = smooth_timepoints, intercept = FALSE)
  
  if (!is.na(spline_test_configs$dof[[config_index]])) {
    args$df <- spline_test_configs$dof[[config_index]]
  } else {
    args$knots <- spline_test_configs$knots[[config_index]]
  }
  
  if (!is.na(spline_test_configs$bknots[[config_index]])) {
    args$Boundary.knots <- spline_test_configs$bknots[[config_index]]
  }
  
  
  if (spline_test_configs$spline_type[config_index] == "b") {
    args$degree <- spline_test_configs$degree[[config_index]]
    X <- do.call(splines::bs, args)
  } else {                                          # natural cubic splines
    X <- do.call(splines::ns, args)
  }
  
  # Generate all the individual plots
  for(index in indices) {
    # X <- splines::ns(smooth_timepoints, df = DoF, intercept = FALSE)
    
    DoF <- which(names(top_table) == "AveExpr") - 1
    spline_coeffs <- 
      as.numeric(top_table[top_table$feature_index == 
                             index, paste0("X", 1:DoF)])
    
    intercept <- 
      as.numeric(top_table$intercept[top_table$feature_index == index])
    
    fitted_values <- X %*% spline_coeffs + intercept
    
    plot_data <- data.frame(Time = meta$Time, 
                            Intensity = as.vector(t(data[index, ])))
    
    plot_spline <- data.frame(Time = smooth_timepoints, Fitted = fitted_values)
    
    # Calculate the extension for the x-axis
    x_max <- max(meta$Time)
    x_extension <- x_max * 0.05  # Extend the x-axis by 5% of its maximum value
    
    p <- ggplot2::ggplot() +
      geom_point(data = plot_data, aes(x = !!rlang::sym("Time"),
                                       y = !!rlang::sym("Intensity")), 
                 color = 'blue') +
      geom_line(data = plot_spline, aes(x = !!rlang::sym("Time"),
                                        y = !!rlang::sym("Fitted")), 
                color = 'red') +
      theme_minimal() +
      scale_x_continuous(limits = c(min(meta$Time), x_max + x_extension))
    labs(x = paste0("Time ", time_unit_label), y = "Intensity")
    
    title <- top_table$feature_names[index]  
    if (is.na(title)) {
      title <- paste("Feature:", index)
    }
    
    p <- p + labs(title = title, 
                  x = paste0("Time ", time_unit_label), y = "Intensity") +
      theme(plot.title = element_text(size = 4),
            axis.title.x = element_text(size = 8), 
            axis.title.y = element_text(size = 8)) +
      annotate("text", x = x_max + (x_extension / 2), y = 
                 max(fitted_values, na.rm = TRUE),
               label = "",
               hjust = 0.5, vjust = 1, size = 3.5, angle = 0, color = "black")

    plot_list[[length(plot_list) + 1]] <- p
  }
  
  if(length(plot_list) > 0) {
    # Generate the combined plot
    num_plots <- length(plot_list)
    ncol <- 3
    composite_plot_len <- as.integer(ceiling(num_plots / ncol))
    
    composite_plot <- patchwork::wrap_plots(plot_list, ncol = 3) + 
      patchwork::plot_annotation(title = paste(top_table_name, type, 
                                               sep = " | "), 
                      theme = theme(plot.title = element_text(hjust = 0.5, 
                                                              size = 14)))
    
    return(list(composite_plot = composite_plot, 
                composite_plot_len = composite_plot_len))
  } else {
    return(FALSE)
  }
}


#' Build Hyperparameters Screening Report
#'
#' @description
#' Constructs an HTML report for hyperparameter screening by embedding plots 
#' and their respective sizes into the provided header section.
#'
#' @param header_section A character string containing the HTML header section.
#' @param plots A list of ggplot2 plot objects.
#' @param plots_sizes A list of integers specifying the number of rows for each 
#' plot.
#' @param output_file_path A character string specifying the path to save the 
#' HTML report.
#'
#' @return No return value, called for side effects.
#'
#' @seealso
#' \code{\link{plot2base64}}
#' 
build_hyperparams_screen_report <- function(header_section, 
                                            plots, 
                                            plots_sizes, 
                                            output_file_path) {
  
  index <- 1
  for (plot in plots) {
    plot_nrows <- plots_sizes[[index]]
    index <- index + 1
    img_tag <- plot2base64(plot, plot_nrows)
    header_section <- paste(header_section, img_tag, sep="\n")
  }
  
  # Close the HTML document
  html_content <- paste(header_section, "</body></html>", sep="\n")
  
  dir_path <- dirname(output_file_path)
  
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  writeLines(html_content, output_file_path)
}
