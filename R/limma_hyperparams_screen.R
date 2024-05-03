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
#' @param metas A list of data frames containing metadata for each dataset in 
#' `datas`.
#' @param designs A character vector of design formulas for the limma analysis.
#' @param modes A character vector indicating the mode of analysis for each 
#' design.
#' @param factors A character vector of factors to be considered in the 
#' analysis.
#' @param DoFs An integer vector of degrees of freedom to be used for each 
#' factor in `factors`.
#' @param feature_names A character vector of feature names to be analyzed.
#' @param pthresholds A numeric vector of p-value thresholds for significance 
#' determination.
#' @param padjust_method A character string specifying the method for p-value 
#' adjustment.
#'
#' @return Returns a list of plots generated from the limma analysis results.
#'         Each element in the list corresponds to a different combination of 
#'         hyperparameters.
#'
#' @examples
#' \dontrun{
#'   results <- limma_hyperparams_screen(datas, metas, designs, modes, factors,
#'                                       DoFs = c(2L), feature_names, 
#'                                       pthresholds = c(0.05),
#'                                       padjust_method = "BH")
#' }
#' 
#' @importFrom here here
#'
#' @export
#'
limma_hyperparams_screen <- function(datas, 
                                     datas_descr,
                                     metas, 
                                     designs, 
                                     modes, 
                                     condition, 
                                     spline_configs,
                                     feature_names, 
                                     report_dir = here::here(),
                                     adj_pthresh = c(0.05),
                                     padjust_method = "BH") {
  
  control_inputs_hyperpara_screen(datas,
                                  datas_descr,
                                  metas,
                                  designs,
                                  modes,
                                  condition,
                                  spline_configs,
                                  feature_names,
                                  adj_pthresh,
                                  padjust_method)
  
  top_tables_combos <- get_limma_combos_results(datas, 
                                                metas, 
                                                designs, 
                                                modes, 
                                                condition, 
                                                spline_configs,
                                                feature_names, 
                                                adj_pthresh, 
                                                padjust_method)
  
  combo_pair_plots <- plot_limma_combos_results(top_tables_combos, 
                                                datas, 
                                                metas)
  
  timestamp <- format(Sys.time(), "%d_%m_%Y-%H_%M_%S")
  
  generate_reports(combo_pair_plots, 
                   report_dir,
                   timestamp)
  
  generate_reports_meta(datas_descr, 
                        designs, 
                        modes, 
                        spline_configs,
                        report_dir,
                        timestamp)
}


# Level 1 internal functions ---------------------------------------------------


control_inputs_hyperpara_screen <- function(datas, 
                                            datas_descr,
                                            metas, 
                                            designs, 
                                            modes, 
                                            condition, 
                                            spline_params,
                                            feature_names, 
                                            pthresholds, 
                                            padjust_method) {
  
  if (!is.list(datas) || any(!sapply(datas, is.matrix))) {
    stop("'datas' must be a list of matrices.")
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
  
  if (!is.list(metas) || any(!sapply(metas, is.data.frame))) {
    stop("'metas' must be a list of dataframes.")
  }
  
  if (length(metas) != length(datas)) {
    stop("The lists 'metas' and 'datas' must have the same number of elements.")
  }
  
  for (i in seq_along(metas)) {
    # Extract the current meta dataframe and data matrix
    meta <- metas[[i]]
    data <- datas[[i]]
    
    # Check if the number of rows in meta equals the number of columns in data
    if (!(nrow(meta) == ncol(data))) {
      stop(paste("Mismatch found for pair", i, ": The number of rows in meta 
                 does not match the number of columns in data. This is 
                 required."))
    }
  }
  
  if (!is.character(designs) || 
      any(nchar(designs) > 75) || 
      all(grepl("X", designs) == FALSE)) {
    
    long_or_missing_x_indices <- 
      which(nchar(designs) > 75 | grepl("X", designs) == FALSE)
    long_elements <- 
      designs[long_or_missing_x_indices[nchar(
        designs[long_or_missing_x_indices]) > 75]]
    missing_x_elements <- 
      designs[long_or_missing_x_indices[grepl(
        "X", designs[long_or_missing_x_indices]) == FALSE]]
    
    error_message_parts <- character()
    if (length(long_elements) > 0) {
      error_message_parts <- c(error_message_parts, sprintf(
        "Some 'designs' elements exceed 75 characters. Offending element(s) at 
        indices %s: '%s'.",
        paste(which(nchar(designs) > 75), collapse=", "),
        paste(long_elements, collapse="', '")
      ))
    }
    if (length(missing_x_elements) > 0) {
      error_message_parts <- c(error_message_parts, sprintf(
        "Some 'designs' elements do not contain 'X', which is the required
        time component for the splinetime package. 
        Offending element(s) at indices %s: '%s'.",
        paste(which(grepl("X", designs) == FALSE), collapse = ", "),
        paste(missing_x_elements, collapse = "', '")
      ))
    }
    
    error_message <- paste(
      "'designs' must be a character vector with no element over 75 characters 
      and each must contain 'X'.",
      paste(error_message_parts, collapse=" Please shorten the variable names 
            or ensure 'X' is included.")
    )
    stop(error_message)
  }
  
  if (!is.character(modes) || !all(modes %in% c("isolated", "integrated"))) {
    stop("'modes' must be a character vector containing only 'isolated' or 
         'integrated'.")
  }
  
  if (!is.character(condition) && !(length(condition) == 1)) {
    stop("'condition' must be a single character")
  }
  
  
  if ("spline_type" %in% names(spline_params)) {
    if (!all(spline_params$spline_type %in% c("b", "n"))) {
      stop("Elements of spline_type must be either 'b' for B-splines or 'n'. for
           natural cubic splines")
    }
  } else {
    stop("spline_type is missing.")
  }
  
  # Check if degrees exists and is an integer vector
  if ("degrees" %in% names(spline_params)) {
    if (!all(spline_params$degrees == as.integer(spline_params$degrees))) {
      stop("degrees must be an integer vector.")
    }
  } else {
    stop("degrees is missing.")
  }
  
  # Check if DoFs exists and is an integer vector
  if ("DoFs" %in% names(spline_params)) {
    if (!all(spline_params$DoFs == as.integer(spline_params$DoFs))) {
      stop("DoFs must be an integer vector.")
    }
  }
  
  # Check if knots exists and is a list of numeric vectors
  if ("knots" %in% names(spline_params)) {
    if (!is.list(spline_params$knots) || 
        any(sapply(spline_params$knots, function(x) !is.numeric(x)))) {
      stop("knots must be a list of numeric vectors.")
    }
  }
  
  if (("DoFs" %in% names(spline_params)) && 
      ("knots" %in% names(spline_params))) {
    stop("Either DoFs or knots must be present, but not both.")
  } else if (!("DoFs" %in% names(spline_params)) && 
             !("knots" %in% names(spline_params))) {
    stop("At least one of DoFs or knots must be present.")
  }
  
  # Check if bknots exists and is a list of numeric vectors
  if ("bknots" %in% names(spline_params)) {
    if (!is.list(spline_params$bknots) || 
        any(sapply(spline_params$bknots, function(x) !is.numeric(x)))) {
      stop("bknots must be a list of numeric vectors.")
    }
  }
  
  check_spline_configs_element_len(spline_configs)
  
  
  if (!is.character(feature_names)) {
    stop("'feature_names' must be a character vector.")
  }
  
  if (!is.numeric(pthresholds) || 
      any(pthresholds <= 0) || 
      any(pthresholds >= 1)) {
    stop("'pthresholds' must be a numeric vector with 
         all elements > 0 and < 1.")
  }
  
  
  if (!is.character(padjust_method) || length(padjust_method) != 1) {
    stop("'padjust_method' must be a single character string.")
  }
  
  # Ensure that datas and metas have the same length, as well as designs and 
  # modes
  if(length(datas) != length(metas) || length(designs) != length(modes)) {
    stop("datas and metas, designs and modes must have the same length.")
  }
}


#' @importFrom tidyr expand_grid
#' @importFrom dplyr mutate
#' @importFrom purrr pmap
#' @importFrom purrr set_names
#' 
get_limma_combos_results <- function(datas, 
                                     metas, 
                                     designs, 
                                     modes, 
                                     condition, 
                                     spline_configs,
                                     feature_names, 
                                     pthresholds, 
                                     padjust_method) {
  combos <- tidyr::expand_grid(
    data_index = seq_along(datas),
    design_index = seq_along(designs),
    spline_config_index = seq_along(spline_configs$spline_type),
    pthreshold = pthresholds
  ) %>% 
    dplyr::mutate(id = paste0("Data_", data_index, 
                       "_Design_", design_index, 
                       "_SConfig_", spline_config_index, 
                       "_PThresh_", pthreshold))
  
  # Define the function to process each combination
  process_combo <- function(data_index, 
                            design_index, 
                            spline_config_index, 
                            pthreshold, 
                            ...) {
    data <- datas[[data_index]]
    meta <- metas[[data_index]]
    design <- designs[[design_index]]
    mode <- modes[[design_index]]
    
    spline_params <- create_spline_params(config_list = spline_configs, 
                                          index = spline_config_index, 
                                          meta = meta, 
                                          condition = condition, 
                                          mode = mode)
    
    result <- run_limma_splines(data = data, 
                                meta = meta, 
                                design = design, 
                                spline_params = spline_params, 
                                condition = condition,
                                feature_names = feature_names, 
                                mode = mode, 
                                padjust_method = padjust_method)
    
    result$top_tables
  }
  
  purrr::pmap(combos, process_combo) %>% 
    purrr::set_names(combos$id)
}


#' @importFrom stringr str_extract
#' @importFrom progress progress_bar
#' @importFrom purrr set_names
#' @importFrom purrr map
#' 
plot_limma_combos_results <- function(all_combos_top_tables,
                                      datas,
                                      metas) {
  
  names_extracted <- stringr::str_extract(names(all_combos_top_tables),
                                          "Data_\\d+_Design_\\d+")
  
  combos_separated <- lapply(unique(names_extracted), function(id) {
    all_combos_top_tables[names_extracted == id]
  })
  
  names(combos_separated) <- unique(names_extracted)
  
  combos <- names(combos_separated)
  combo_pairs <- combn(combos, 2, simplify = FALSE)
  
  print("Generating the plots for all pairwise hyperparams-combo comparisons")
  progress_ticks <- length(combo_pairs)
  pb <- progress::progress_bar$new(total = progress_ticks, 
                                   format = "[:bar] :percent")
  pb$tick(0)
  
  combo_pair_results <- purrr::set_names(
    purrr::map(combo_pairs, function(pair) {
      combo_pair <- combos_separated[pair]
      
      hitcomp <- gen_hitcomp_plots(combo_pair)
      
      composites <- purrr::map(combo_pair, function(combo) {
        composite <- gen_composite_spline_plots(combo,
                                                datas,
                                                metas)
      })
      pb$tick()
      list(hitcomp = hitcomp, composites = composites)
    }
    ), purrr::map(combo_pairs, function(pair) paste(pair[1], "vs", pair[2],
                                             sep = "_"))
  )
}


#' @importFrom progress progress_bar
#' @importFrom purrr imap
#' 
generate_reports <- function(combo_pair_plots, 
                             report_dir,
                             timestamp) {
  print("Building .html reports for all pairwise hyperparams-combo comparisons")
  progress_ticks <- length(combo_pair_plots)
  pb <- progress::progress_bar$new(total = progress_ticks, 
                                   format = "[:bar] :percent")
  
  result <- purrr::imap(combo_pair_plots, ~{
    process_combo_pair(.x, .y, report_dir, timestamp)
    pb$tick()  
  })
}


#' @importFrom here here
#' @importFrom knitr kable
#' @importFrom kableExtra kable_styling
#' 
generate_reports_meta <- function(datas_descr, 
                                  designs, 
                                  modes, 
                                  spline_configs,
                                  report_dir,
                                  timestamp) {
  formatted_spline_configs <- flatten_spline_configs(spline_configs)
  
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


check_spline_configs_element_len <- function(spline_configs) {
  expected_length <- length(spline_configs$spline_type)
  
  # Check if all elements have the correct length
  for (element in names(spline_configs)) {
    if (length(spline_configs[[element]]) != expected_length) {
      stop(paste("Error:", element, "length does not match 'spline_type' 
                 length."))
    }
  }
}


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


#' @importFrom tidyr expand_grid unnest_longer replace_na pivot_wider
#' @importFrom tibble enframe column_to_rownames
#' @importFrom dplyr mutate select left_join
#' @importFrom pheatmap pheatmap
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
    tidyr::unnest_longer(features) %>%
    dplyr::mutate(x1 = 1)
  
  df_2 <-
    hits_2 %>% 
    tibble::enframe("params", "features") %>% 
    tidyr::unnest_longer(features) %>%
    dplyr::mutate(x2 = 2)
  
  venn_matrix <- 
    df %>% 
    dplyr::left_join(df_1, by = c("features", "params")) %>% 
    dplyr::left_join(df_2, by = c("features", "params")) %>% 
    tidyr::replace_na(list(x1 = 0, x2 = 0)) %>% 
    dplyr::mutate(x = x1 + x2) %>% 
    dplyr::select(!c(x1, x2)) %>%
    tidyr::pivot_wider(names_from = params, values_from = x) %>% 
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


#' @importFrom ggplot2 geom_col geom_text facet_wrap scale_y_continuous 
#' @importFrom ggplot2 theme_minimal element_text element_blank expansion xlab
#' @importFrom ggplot2 theme
#' @importFrom purrr map_int set_names 
#' @importFrom tibble enframe
#' @importFrom dplyr bind_rows
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
  
  ggplot2::ggplot(plot_data, aes(x = params, y = n_hits)) +
    geom_col() +
    geom_text(
      aes(label = n_hits),
      vjust = -0.5,
      
    ) +
    xlab("Parameters") +
    scale_y_continuous(
      "number of significant features",
      expand = expansion(mult = c(0, .2))
    ) +
    facet_wrap(vars(condition), ncol = 1) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      panel.grid.major.x = element_blank()
    ) +
    NULL
}


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


#' One half of one condition comparison HTML 
#' (composite spline plots for one 'condition' inside one condition comparison)
#' @importFrom utils tail
#' 
gen_composite_spline_plots <- function(internal_combos, 
                                       datas, 
                                       metas) {
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
                                         top_table, 
                                         combo_name, 
                                         indices,
                                         type)
        
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


process_combo_pair <- function(combo_pair, 
                               combo_pair_name, 
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
  generate_report_html_new(plots, 
                           plots_len, 
                           report_dir, 
                           combo_pair_name,
                           timestamp)
}


create_spline_params <- function(config_list, 
                                 index, 
                                 meta, 
                                 condition, 
                                 mode) {
  num_levels <- length(unique(meta[[condition]]))  
  
  process_item <- function(item, index, num_levels, mode) {
    if (mode == "integrated") {
      # Just take the index element from each item or sub-item
      if (is.list(item)) {
        list(item[[index]])
      } else {
        item[index]
      }
    } else {     # mode = 'isolated'
      if (is.list(item)) {
        rep(list(item[[index]]), num_levels)
      } else {
        rep(item[index], num_levels)
      }
    } 
  }
  result <- lapply(config_list, process_item, index, num_levels, mode)
}


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


# Level 3 internal functions  ---------------------------------------------------


#' @importFrom dplyr pull
store_hits <- function(condition) {
  hits_cond <- list()
  
  for (item in condition) {
    df <- item$DataFrame
    adj_p_value_treshold <- item$Parameters$adj_p_value_threshold
    key_name <- item$Parameters$params_id
    
    # Filters df for only features with adjP < threshold
    hits_cond[[key_name]] <-
      df %>%
      subset(adj.P.Val < adj_p_value_treshold) %>% 
      dplyr::pull(feature_index) %>% 
      as.character()
  }
  
  return(hits_cond)
}


#' @importFrom splines ns
#' @importFrom ggplot2 ggplot geom_point geom_line theme_minimal 
#' @importFrom ggplot2 scale_x_continuous labs annotate theme
#' @importFrom patchwork wrap_plots plot_annotation
#' 
plot_composite_splines <- function(data, 
                                   meta, 
                                   top_table, 
                                   top_table_name, 
                                   indices,
                                   type) {
  plot_list <- list()
  DoF <- as.integer(sub(".*SConfig_([0-9]+)_.*", "\\1", top_table_name)) 
  
  # Generate all the individual plots
  for(index in indices) {
    smooth_timepoints <- seq(meta$Time[1], meta$Time[length(meta$Time)], 
                             length.out = 100)
    X <- splines::ns(smooth_timepoints, df = DoF, intercept = FALSE)
    
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
      geom_point(data = plot_data, aes(x = Time, y = Intensity), 
                 color = 'blue') +
      geom_line(data = plot_spline, aes(x = Time, y = Fitted), 
                color = 'red') +
      theme_minimal() +
      scale_x_continuous(limits = c(min(meta$Time), x_max + x_extension))
    labs(x = "Time [min]", y = "Intensity")
    
    title <- top_table$feature_names[index]  
    if (is.na(title)) {
      title <- paste("Feature:", index)
    }
    
    p <- p + labs(title = title, 
                  x = "Time [min]", y = "Intensity") +
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
