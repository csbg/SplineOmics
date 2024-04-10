# Library imports ---------------------------------------------------------
library(pheatmap)
library(ggplot2)
library(stats)
library(purrr)
library(dplyr)
library(patchwork)



# Internal functions level 3 ----------------------------------------------


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
      pull(feature_index) %>% 
      as.character()
  }
  
  return(hits_cond)
}


plot_composite_splines <- function(data, meta, top_table, top_table_name, 
                                   indices) {

  plot_list <- list()
  DoF <- as.integer(sub(".*DoF_([0-9]+)_.*", "\\1", top_table_name)) 
  
  # Generate all the individual plots
  for(index in indices) {
    smooth_timepoints <- seq(meta$Time[1], meta$Time[length(meta$Time)], 
                             length.out = 100)
    X <- ns(smooth_timepoints, df = DoF, intercept = FALSE)
    
    spline_coeffs <- 
      top_table[top_table$feature_index == index, paste0("X", 1:DoF)]
    
    intercept_value <- top_table$intercept[top_table$feature_index == index]
    
    fitted_values <- X %*% spline_coeffs + intercept
    
    plot_data <- data.frame(Time = meta$Time, 
                            Y = data[index, ])
    
    plot_spline <- data.frame(Time = smooth_timepoints,
                              Fitted = fitted_values)
    
    # Calculate the extension for the x-axis
    x_max <- max(meta$Time)
    x_extension <- x_max * 0.05  # Extend the x-axis by 5% of its maximum value
    
    p <- ggplot() +
      geom_point(data = plot_data, aes(x = Time, y = Y), color = 'blue') +
      geom_line(data = plot_spline, aes(x = Time, y = Fitted), 
                color = 'red') +
      theme_minimal() +
      scale_x_continuous(limits = c(min(timePoints), x_max + x_extension))
    labs(x = "Time [min]", y = "Intensity")
    

    title <- top_table$feature_name[index]  
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
    nrow <- ceiling(num_plots / ncol)
    
    composite_plot <- patchwork::wrap_plots(plot_list, ncol = 3) + 
      plot_annotation(title = paste(data_descriptor, "| DoF:", 
                                    DoF, 
                                    "| p-value threshold:", 
                                    p_value_threshold), 
                      theme = theme(plot.title = element_text(hjust = 0.5, 
                                                              size = 14)))
    return(list(composite_plot = composite_plot, composite_plot_len = nrow))
  } else {
    return(FALSE)
  }
}



# Internal functions level 2 ----------------------------------------------


hc_new <- function(cond1name = "Condition 1", cond2name = "Condition 2") {
  
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


hc_add <- function(hc_obj, top_table, params_id, condition = 1, threshold = 0.05){
  # Validate input
  if (!is.data.frame(top_table)) stop("top_table must be a dataframe.")
  if (!is.character(params_id) || nchar(params_id) > 50) stop("max len = 25")
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


hc_vennheatmap <- function(hc_obj) {
  hits_1 <- store_hits(hc_obj$data[[1]])
  hits_2 <- store_hits(hc_obj$data[[2]])
  
  df <- expand_grid(
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
    enframe("params", "features") %>% 
    unnest_longer(features) %>%
    mutate(x1 = 1)
  
  df_2 <-
    hits_2 %>% 
    enframe("params", "features") %>% 
    unnest_longer(features) %>%
    mutate(x2 = 2)

  venn_matrix <- 
    df %>% 
    left_join(df_1, by = c("features", "params")) %>% 
    left_join(df_2, by = c("features", "params")) %>% 
    replace_na(list(x1 = 0, x2 = 0)) %>% 
    mutate(x = x1 + x2) %>% 
    select(!c(x1, x2)) %>%
    pivot_wider(names_from = params, values_from = x) %>% 
    column_to_rownames("features") %>% 
    as.matrix()
  
  venn_matrix <- venn_matrix[, order(colnames(venn_matrix))]
  
  color_palette <- c("white", "blue", "yellow", "green")
  breaks <- c(-0.5, 0.5, 1.5, 2.5, 3.5)
  
  plot_title <- sprintf("0 -> none, 1 -> %s, 2 -> %s, 3 -> both", 
                        hc_obj$condition_names[[1]],
                        hc_obj$condition_names[[2]])
  
  # Generating the heatmap plot
  vennheatmap_plot <- pheatmap::pheatmap(venn_matrix, color = color_palette,
                                         breaks = breaks,
                                         cluster_cols = FALSE,
                                         cluster_rows = TRUE,
                                         show_rownames = TRUE,
                                         show_colnames = TRUE,
                                         border_color = NA,
                                         main = plot_title)
  
  return(list(vennheatmap = vennheatmap_plot, 
              nrhits = nrow(venn_matrix)))
}


hc_barplot <- function(hc_obj) {
  plot_data <- 
    list(
      store_hits(hc_obj$data[[1]]) %>% 
        map_int(length) %>% 
        enframe("params", "n_hits"),
      store_hits(hc_obj$data[[2]]) %>% 
        map_int(length) %>% 
        enframe("params", "n_hits")
    ) %>% 
    set_names(hc_obj$condition_names) %>% 
    bind_rows(.id = "condition")
  
  ggplot(plot_data, aes(x = params, y = n_hits)) +
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


gen_hitcomp_plots <- function(all_combos_top_tables, plots, plots_len) {
  hitcomp <- hc_new("combined", "individual")
  
  combo_names <- sapply(names(all_combos_top_tables), function(name) {
    # Using sub to extract everything from 'DoF' onwards
    match <- sub(".*(DoF.*)$", "\\1", name)
    return(match)
  })
  
  condition_nr <- 1L
  for (i in 1:length(all_combos_top_tables)) {
    combo_top_tables <- all_combos_top_tables[[i]]
    combo_name <- combo_names[i]
    # Extract the part where DoF and p-value threshold are written.
    pthreshold <- as.numeric(sub(".*_([^_]+)$", "\\1", combo_name))
    
    if (i == length(all_combos_top_tables)/2 + 1) {
      condition_nr <- 2L
    }
    
    for (top_table_name in names(combo_top_tables)) {
      top_table <- combo_top_tables[[top_table_name]]
      id <- paste(top_table_name, combo_name, sep = "_")
      hitcomp <- hc_add(hitcomp, top_table, id, condition_nr, pthreshold)
    }
  }
  
  result <- hc_vennheatmap(hitcomp)

  plots[[length(plots) + 1]] <- result$vennheatmap
  plots_len[[length(plots_len) + 1]] <- result$nrhits
  
  
  barplot <- hc_barplot(hitcomp)
  
  plots[[length(plots) + 1]] <- barplot
  plots_len[[length(plots_len) + 1]] <- 2
  
  return(list(plots = plots, plots_len = plots_len))
}


gen_composite_spline_plots <- function(all_combos_top_tables, datas, metas, 
                                      plots, plots_len) {
  for (combo_name in names(all_combos_top_tables)) {
    top_tables_levels <- all_combos_top_tables[[combo_name]]
    
    data <- datas[[as.integer(strsplit(combo_name, "_")[[1]][2])]]
    meta <- metas[[as.integer(strsplit(combo_name, "_")[[1]][2])]]
    
    ptresh <- as.numeric(tail(strsplit(combo_name, "_")[[1]], 1)) 
    
    for (top_table_name in names(top_tables_levels)) {
      top_table <- top_tables_levels[[top_table_name]]
      
      parts <- strsplit(top_table_name, "_")[[1]]
      factor <- parts[1]
      level <- parts[2]
      meta_level <- meta[meta[[column_name]] == level, ]
      
      for(type in c('significant', 'not_significant')) {
        
        if (type == "significant") {
          filtered_rows <- top_table[top_table$adj.P.Val < pthresh, ]
          selected_rows <- if(nrow(filtered_rows) > 30) {
            filtered_rows[sample(nrow(filtered_rows), 30), ]
          } else {
            filtered_rows
          }
          indices <- as.integer(selected_rows$feature_index)
        } else if (type == "not_significant") {
          filtered_rows <- top_table[top_table$adj.P.Val >= pthresh, ]
          selected_rows <- if(nrow(filtered_rows) > 30) {
            filtered_rows[sample(nrow(filtered_rows), 30), ]
          } else {
            filtered_rows
          }
          indices <- as.integer(selected_rows$feature_index)
        }
         
        result <- plot_composite_splines(data, meta_level, top_table, 
                                         combo_name, indices)
        
        if (is.list(result)) {
          plots[[length(plots) + 1]] <- result$composite_plot
          plots_len[[length(plots_len) + 1]] <- result$composite_plot_len
        }
      }
    }
  }
  
  return(list(plots = plots, plots_len = plots_len))
}



# Internal functions level 1 ----------------------------------------------


control_inputs_hyperpara_screen <- function(datas, metas, designs, modes, 
                                            factors, DoFs, feature_names, 
                                            pthresholds, padjust_method) {
  
  if (!is.list(datas) || any(!sapply(datas, is.data.frame))) {
    stop("'datas' must be a list of dataframes.")
  }
  
  if (!is.list(metas) || any(!sapply(metas, is.data.frame))) {
    stop("'metas' must be a list of dataframes.")
  }
  
  if (!is.character(designs)) {
    stop("'designs' must be a character vector.")
  }
  
  if (!is.character(modes)) {
    stop("'modes' must be a character vector.")
  }
  
  if (!is.character(factors)) {
    stop("'factors' must be a character vector.")
  }
  
  if (!is.integer(DoFs)) {
    stop("'DoFs' must be an integer vector.")
  }
  
  if (!is.character(feature_names)) {
    stop("'feature_names' must be a character vector.")
  }
  
  if (!is.numeric(pthresholds)) {
    stop("'pthresholds' must be a numeric vector.")
  }
  
  if (!is.character(padjust_method) || length(padjust_method) != 1) {
    stop("'padjust_method' must be a single character string.")
  }
  
  # Ensure that datas and metas have the same length, as well as designs and modes
  if(length(datas) != length(metas) || length(designs) != length(modes)) {
    stop("Datas and Metas, Designs and Modes must have the same length.")
  }
}


get_limma_combos_results <- function (datas, metas, designs, modes, factors,
                                      DoFs, feature_names, pthresholds,
                                      padjust_method) {
  top_tables_combos <- list()
  
  for(i in seq_along(datas)) {
    data <- datas[[i]]
    meta <- metas[[i]]
    
    for(j in seq_along(designs)) {
      design <- designs[[j]]
      mode <- modes[[j]]
      
      for(DoF in DoFs) {
        for(pthreshold in pthresholds) {
          result <- run_limma_splines(data, meta, design, DoF, factors,
                                      feature_names, mode, padjust_method)
          
          # Construct a unique id for the current combination of hyperparameters
          id <- paste0("DataMeta_", i, "_DesignMode_", j, "_DoF_", DoF, 
                       "_PThresh_", pthreshold)
          
          top_tables_combos[[id]] <- result$top_tables
        }
      }
    }
  }
  
  return(top_tables_combos)
}


plot_limma_combos_results <- function(all_combos_top_tables, datas, metas, 
                                      annotation) {
  plots <- list()
  plots_len <- list()
  
  result <- gen_hitcomp_plots(all_combos_top_tables, plots, plots_len)

  result <-
    gen_composite_spline_plots(all_combos_top_tables, datas, metas, 
                               result$plots, result$plots_len)
  
  return(list(plots = result$plots, plots_len = result$plots_len))
}



# Exported functions ----------------------------------------


#' Limma Hyperparameters Screening
#'
#' This function screens through various combinations of hyperparameters for limma analysis,
#' including designs, modes, and degrees of freedom. It validates inputs, generates results for all
#' combinations, and plots the outcomes. Finally, it may also be involved in generating an HTML report
#' as part of a larger analysis workflow.
#'
#' @param datas A list of data frames containing the datasets to be analyzed.
#' @param metas A list of data frames containing metadata for each dataset in `datas`.
#' @param designs A character vector of design formulas for the limma analysis.
#' @param modes A character vector indicating the mode of analysis for each design.
#' @param factors A character vector of factors to be considered in the analysis.
#' @param DoFs An integer vector of degrees of freedom to be used for each factor in `factors`.
#' @param feature_names A character vector of feature names to be analyzed.
#' @param pthresholds A numeric vector of p-value thresholds for significance determination.
#' @param padjust_method A character string specifying the method for p-value adjustment.
#'
#' @return Returns a list of plots generated from the limma analysis results.
#'         Each element in the list corresponds to a different combination of hyperparameters.
#'
#' @examples
#' \dontrun{
#'   results <- limma_hyperparams_screen(datas, metas, designs, modes, factors,
#'                                       DoFs = c(2L), feature_names, pthresholds = c(0.05),
#'                                       padjust_method = "BH")
#' }
#'
#' @export
#' @importFrom magrittr "%>%"
#' @importFrom limma lmFit eBayes topTable
#'
limma_hyperparams_screen <- function(datas, metas, designs, modes, 
                                     factors, DoFs = c(2L), feature_names, 
                                     pthresholds = c(0.05),
                                     padjust_method = "BH") {
  
  control_inputs_hyperpara_screen(datas, metas, designs, modes, factors, DoFs, 
                                  feature_names, pthresholds, padjust_method)
  
  top_tables_combos <- get_limma_combos_results(datas, metas, designs, modes, 
                                                factors, DoFs, feature_names, 
                                                pthresholds, padjust_method)
  
  # debug(plot_limma_combos_results)
  result <- plot_limma_combos_results(top_tables_combos, datas, metas)
  
  for (plot in result) {
    print(plot)
  }
  
  # generate_report()
  
}
