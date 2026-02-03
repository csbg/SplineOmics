#' Create PVC HTML report from statistical PVC results
#'
#' @description
#' Generates per-feature PVC plots and writes an HTML report to disk from the
#' statistical output of `find_pvc()`. This function does not recompute any
#' statistics; it consumes `pvc_results` (including its attributes) and the
#' original input data stored in `splineomics`.
#'
#' @param splineomics `list`: Preprocessed time-series input. Must include at
#'   least \code{data}, \code{meta}, \code{annotation}, \code{condition},
#'   \code{meta_batch_column}, \code{meta_batch2_column} (optional),
#'   \code{report_info}, and \code{feature_name_columns}.
#'
#' @param pvc_results `list`: Output from `find_pvc()`, a named list by
#'   condition level. Each level must contain \code{alpha} and
#'   \code{pvc_adj_pvals}. Attributes \code{padjust_method} and \code{support}
#'   are used for report header settings.
#'
#' @param plot_info `list`: Plot annotation options passed to `plot_pvc()`.
#'   See `find_pvc()` for the expected structure. Treatment annotations are
#'   applied per condition level if \code{treatment_labels} and
#'   \code{treatment_timepoints} are provided.
#'   
#' @param verbose Boolean flag indicating if messages should be shown.
#'
#' @param report_dir `character(1)`: Output directory for the HTML report and
#'   any associated files.
#'
#' @return The input \code{pvc_results} with an additional \code{plots} element
#'   added under each condition level. The report is written to
#'   \code{report_dir}.
#'
#' @details
#' For each condition level, the function subsets \code{splineomics$data} and
#' \code{splineomics$meta}, generates per-feature plots via `plot_pvc()`, then
#' calls `generate_report_html()` with the plots, metadata, and report header
#' settings derived from \code{pvc_results} and \code{splineomics}.
#'
#' @export
#'
create_pvc_report <- function(
        splineomics,
        pvc_results,
        plot_info = list(
            y_axis_label = "Value",
            time_unit = "min",
            treatment_labels = NA,
            treatment_timepoints = NA
        ),
        verbose = FALSE,
        report_dir = tempdir()
) {
    check_splineomics_elements(
        splineomics = splineomics,
        func_type = "find_peaks_valleys"
    )

    data <- splineomics[["data"]]
    meta <- splineomics[["meta"]]
    annotation <- splineomics[["annotation"]]
    condition <- splineomics[["condition"]]
    meta_batch_column <- splineomics[["meta_batch_column"]]
    meta_batch2_column <- splineomics[["meta_batch2_column"]]
    report_info <- splineomics[["report_info"]]
    feature_name_columns <- splineomics[["feature_name_columns"]]

    batch_effects <- c(meta_batch_column, meta_batch2_column)
    condition_levels <- names(pvc_results)
    alphas <- purrr::map(
        pvc_results,
        ~ .x[["alpha"]]
    )

    for (level in condition_levels) {
        level_chr <- as.character(level)
        selected_samples <- meta[[condition]] == level
        sub_data <- data[, selected_samples, drop = FALSE]
        sub_meta <- meta[selected_samples, , drop = FALSE]
        
        plots <- plot_pvc(
            pvc_pvals = pvc_results[[level_chr]][["pvc_adj_pvals"]],
            data = sub_data,
            meta = sub_meta,
            meta_batch_column = meta_batch_column,
            alpha = pvc_results[[level_chr]][["alpha"]],
            plot_info = plot_info,
            level = level_chr
        )
        
        pvc_results[[level_chr]][["plots"]] <- plots
    }
    
    padjust_method <- attr(
        pvc_results,
        "padjust_method",
        exact = TRUE
    )
    support <- attr(
        pvc_results,
        "support",
        exact = TRUE
    )
    
    level_headers_info <- list(
        pvc_settings = list(
            adj_value_thresh = alphas,
            padjust_method = padjust_method,
            support = support,
            batch_effects = batch_effects
        )
    )
    
    report_data <- bind_data_with_annotation(
        data,
        annotation
    )
    
    generate_report_html(
        plots = pvc_results,
        plots_sizes = NULL,
        report_info = report_info,
        data = report_data,
        meta = meta,
        level_headers_info = level_headers_info,
        report_type = "find_pvc",
        feature_name_columns = feature_name_columns,
        report_dir = report_dir
    )
    
    if (verbose) {
        print_info_message(
            message_prefix = "PVC report generation",
            report_dir = report_dir
        )
    }

    pvc_results
}


# Level 1 function definitions -------------------------------------------------


#' Plot Peaks and Valleys in Time-Series Omics Data
#'
#' @noRd
#'
#' @description
#' This function generates scatter plots for features that exhibit significant
#' local peaks or valleys (excursions) in time-series omics data. Excursion
#' points are highlighted in red, while normal points remain grey. If the
#' excursion is statistically significant
#' (based on the compound contrast test),
#' a significance star is shown directly above the excursion point.
#'
#' @param results A list returned from `detect_excursions()`, containing
#' `results_df` (excursion matrix) and `pairwise_pvals`
#' (one-tailed p-values for the excursion contrast).
#' @param data A numeric matrix, where rows correspond to features
#' (e.g., genes, proteins, metabolites) and columns correspond to samples.
#' @param meta A data frame containing metadata for the samples. Must include
#' a column named `"Time"` that specifies the timepoint for each sample.
#' @param meta_replicates_column A character string specifying the column name
#' in `meta` that indicates biological replicates.
#' @param alpha A numeric value specifying the significance threshold
#' for displaying stars above excursion points. Defaults to `0.05`.
#' @param plot_info List containing the elements y_axis_label (string),
#'                  time_unit (string), treatment_labels (character vector),
#'                  treatment_timepoints (integer vector). All can also be NA.
#'                  This list is used to add this info to the spline plots.
#'                  time_unit is used to label the x-axis, and treatment_labels
#'                  and -timepoints are used to create vertical dashed lines,
#'                  indicating the positions of the treatments (such as
#'                  feeding, temperature shift, etc.).
#'
#' @return A named list of ggplot objects, where each element corresponds to a
#' feature with at least one detected excursion. Each plot displays the
#' expression levels across timepoints, with replicates distinguished by shape.
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_shape_manual
#'             scale_color_manual geom_text labs theme_minimal
#' @importFrom rlang .data
#'
plot_pvc <- function(
        pvc_pvals,
        data,
        meta,
        meta_batch_column,
        alpha = 0.05,
        plot_info,
        level
        ) {
    peak_valley_flags <- ifelse(
        is.na(pvc_pvals),
        0,
        ifelse(pvc_pvals < alpha, 1, 0)
    )
    
    unique_timepoints <- sort(unique(meta$Time))
    num_timepoints <- length(unique_timepoints)
    
    unique_replicates <- unique(meta[[meta_batch_column]])
    num_replicates <- length(unique_replicates)
    
    symbols_available <- c(21, 22, 23, 24, 25, 7, 8, 10, 12)
    symbols <- symbols_available[seq_len(num_replicates)]
    
    plots <- list()
    
    for (protein_index in which(rowSums(peak_valley_flags) > 0)) {
        feature_name <- rownames(peak_valley_flags)[protein_index]
        protein_data <- data[feature_name, ]
        
        plot_data <- data.frame(
            Time = meta$Time,
            Feature_value = as.numeric(protein_data),
            Replicate = as.factor(meta[[meta_batch_column]])
        )
        
        excursion_flags <- as.numeric(peak_valley_flags[protein_index, ])
        
        # Step 1: define internal timepoints
        internal_timepoints <- unique_timepoints[2:(num_timepoints - 1)]
        
        # Step 2: map excursion flags to only internal timepoints
        named_flags <- excursion_flags
        names(named_flags) <- as.character(internal_timepoints)
        
        # Step 3: assign only to matching timepoints
        plot_data$Excursion <- named_flags[as.character(plot_data$Time)]
        plot_data$Excursion[is.na(plot_data$Excursion)] <- 0
        
        plot_data$Point_Type <- factor(
            ifelse(
                plot_data$Excursion == 1,
                "Excursion",
                "Normal"
            ),
            levels = c("Normal", "Excursion")
        )
        
        # For plotting significance stars directly above excursion points
        sig_df <- data.frame(
            Time = numeric(0),
            Label = character(0),
            y_pos = numeric(0)
        )
        
        get_stars <- function(p, alpha) {
            if (p < alpha / 500) {
                return("****")
            } else if (p < alpha / 50) {
                return("***")
            } else if (p < alpha / 5) {
                return("**")
            } else if (p < alpha) {
                return("*")
            } else {
                return("")
            }
        }
        
        # Loop over internal timepoints
        for (t in 2:(num_timepoints - 1)) {
            timepoint <- unique_timepoints[t]
            flag_index <- t - 1
            
            if (excursion_flags[flag_index] == 1) {
                p_val <- pvc_pvals[protein_index, flag_index]
                
                if (p_val < alpha) {
                    stars <- get_stars(p_val, alpha)
                    
                    max_value <- max(
                        plot_data$Feature_value[plot_data$Time == timepoint],
                        na.rm = TRUE
                    )
                    
                    sig_df <- rbind(
                        sig_df,
                        data.frame(
                            Time = timepoint,
                            Label = stars,
                            max_value = max_value,
                            PValue = p_val
                        )
                    )
                }
            }
        }
        
        # Build the p-value string if any are significant
        if (nrow(sig_df) > 0) {
            pval_str <- paste0(
                "adj.p-val: ",
                paste0(
                    "T=", sig_df$Time, " -> ",
                    formatC(sig_df$PValue, format = "fg", digits = 4),
                    collapse = "; "
                )
            )
            plot_title <- paste(feature_name, "\n", pval_str)
        } else {
            plot_title <- feature_name
        }
        
        p <- ggplot2::ggplot(
            plot_data,
            ggplot2::aes(
                x = !!rlang::sym("Time"),
                y = !!rlang::sym("Feature_value")
            )
        ) +
            ggplot2::geom_point(
                ggplot2::aes(
                    shape = !!rlang::sym("Replicate"),
                    color = !!rlang::sym("Point_Type")
                ),
                size = 3,
                stroke = 1.2,
                fill = "white",
                na.rm = TRUE
            ) +
            ggplot2::scale_shape_manual(values = symbols) +
            ggplot2::geom_text(
                data = sig_df,
                ggplot2::aes(
                    x = !!rlang::sym("Time"),
                    y = !!rlang::sym("max_value"),
                    label = !!rlang::sym("Label")
                ),
                size = 5,
                hjust = 0.5,
                vjust = -0.5
            ) +
            ggplot2::labs(
                title = plot_title,
                x = paste0("Time [", plot_info[["time_unit"]], "]"),
                y = plot_info[["y_axis_label"]],
                color = "Timepoints",
                shape = "Replicates"
            ) +
            ggplot2::theme_minimal()
        
        y_min <- min(plot_data$Feature_value, na.rm = TRUE)
        y_max <- max(plot_data$Feature_value, na.rm = TRUE)
        y_extension <- (y_max - y_min) * 0.1
        y_pos <- y_max + y_extension
        
        result <- maybe_add_dashed_lines(
            p = p,
            plot_info = plot_info,
            level = level,
            y_pos = y_pos
        )
        
        p <- result$p
        treatment_colors <- result$treatment_colors
        
        color_values <- c(
            "Normal" = "grey40",
            "Excursion" = "red",
            treatment_colors
        )
        
        p <- p + ggplot2::scale_color_manual(
            values = color_values,
            guide = ggplot2::guide_legend(
                override.aes = list(
                    size = c(
                        rep(1.5, 2),
                        rep(0.5, length(treatment_colors))
                    )
                )
            )
        )
        
        plots[[feature_name]] <- p
    }
    
    return(plots)
}


#' Build full PVC HTML report
#'
#' @noRd
#'
#' @description
#' This function assembles the complete HTML report for PVC analysis, including
#' section headers, significance metadata, subplot titles, plots, and a table
#'  of contents. It processes a nested list of plots and writes the final
#' report to the specified output file path.
#'
#' @param header_section A character string of HTML content that represents the
#'   top-level header section of the report (e.g., title, description, metadata)
#' @param plots A named list of plot objects, each containing a `plots` field
#'   (itself a named list of ggplot2 objects). The names define report sections.
#' @param level_headers_info A list of metadata associated with each report
#'   section, including p-value thresholds and padjust methods, structured
#'   under the `pvc_settings` key.
#' @param report_info A list containing metadata about the entire report (e.g.,
#'   parameters used, timestamps, version info) used in the report footer.
#' @param output_file_path File path (character) where the final HTML report
#'   should be written. 
#'
#' @return Used for side effects only. The function writes an HTML file to disk.
#'
build_pvc_report <- function(
        header_section,
        plots,
        level_headers_info,
        report_info,
        output_file_path
        ) {
    html_content <- paste(
        header_section,
        "<!--TOC-->",
        sep = "\n"
    )
    
    toc <- create_toc()
    
    styles <- define_html_styles()
    section_header_style <- styles$section_header_style
    toc_style <- styles$toc_style
    
    current_header_index <- 1
    
    # Flatten the nested structure into a list of (header_name, plot) pairs
    flattened_plots <- list()
    for (header_name in names(plots)) {
        subplots <- plots[[header_name]]$plots
        subplot_names <- names(subplots)
        
        for (i in seq_along(subplots)) {
            flattened_plots[[length(flattened_plots) + 1]] <- list(
                header_name = header_name,
                subplot_name = subplot_names[i],
                plot = subplots[[i]]
            )
        }
    }
    
    # Clean up header info
    levels <- names(plots)
    n_plots_per_header <- vapply(plots, function(x) length(x$plots), integer(1))
    
    # Track current section and plot offset
    current_header_index <- 1
    plots_processed_in_section <- 0
    
    # Create progress bar
    pb <- create_progress_bar(flattened_plots)
    
    for (index in seq_along(flattened_plots)) {
        if (plots_processed_in_section == 0) {
            # Insert section header
            level <- levels[[current_header_index]]
            
            nr_of_hits <- n_plots_per_header[[current_header_index]]
            
            html_content <- paste(
                html_content,
                generate_section_header_block(
                    level = levels[[current_header_index]],
                    index = index,
                    section_header_style = section_header_style,
                    level_headers_info = level_headers_info,
                    nr_of_hits = nr_of_hits
                ),
                sep = "\n"
            )
            
            hits_info <- "<p style='text-align: center; font-size: 30px;'>"
            html_content <- paste(
                html_content,
                hits_info,
                sep = "\n"
            )
            
            toc_entry <- sprintf(
                "<li style='%s'><a href='#section%d'>%s</a></li>",
                toc_style,
                index,
                level
            )
            toc <- paste(
                toc,
                toc_entry,
                sep = "\n"
            )
        }
        
        # Add subplot name as an HTML sub-header
        subplot_name <- flattened_plots[[index]]$subplot_name
        subplot_header <- sprintf(
            paste0(
                '<h4 style="text-align: center; font-size: 20px; ',
                'margin-top: 30px;">\n%s\n</h4>'
            ),
            subplot_name
        )
        
        # Inject it into the HTML content before the plot
        html_content <- paste(
            html_content,
            subplot_header,
            sep = "\n"
        )
        
        result <- process_plots(
            plots_element = flattened_plots[[index]]$plot,
            plots_size = 2,
            html_content = html_content,
            toc = toc,
            header_index = current_header_index
        )
        html_content <- result$html_content
        toc <- result$toc
        
        # Update section counters
        plots_processed_in_section <- plots_processed_in_section + 1
        if (plots_processed_in_section >=
            n_plots_per_header[[current_header_index]]) {
            current_header_index <- current_header_index + 1
            plots_processed_in_section <- 0
        }
        
        pb$tick()
    }
    
    generate_and_write_html(
        toc = toc,
        html_content = html_content,
        report_info = report_info,
        output_file_path = output_file_path
    )
}


# Level 2 function definitions -------------------------------------------------


#' Generate HTML header block for a condition level section
#'
#' @noRd
#'
#' @description
#' Constructs an HTML block that includes the section header and related
#' metadata
#' for a given condition level. This includes the adjusted p-value threshold,
#' the p-adjustment method, and a visual legend for significance stars based on
#' the threshold. The output is centered and styled with inline HTML for use in
#' reports or R Markdown rendering.
#'
#' @param level A character string representing the current condition level name
#' @param index An integer index used to create a unique HTML section ID.
#' @param section_header_style A character string of inline CSS styles to apply
#'   to the section header (`<h2>` element).
#' @param level_headers_info A list containing `pvc_settings`, with named
#'   entries for `adj_value_thresh` (a named list of alpha values by level) and
#'   `padjust_method` (a scalar string).
#'
#' @return A character string of HTML content to be included in the full HTML
#'         output.
#'
generate_section_header_block <- function(
        level,
        index,
        section_header_style,
        level_headers_info,
        nr_of_hits
        ) {
    # Access per-level alpha threshold
    alpha_thresh <-
        level_headers_info[["pvc_settings"]][["adj_value_thresh"]][[level]]
    
    # Access global padjust_method
    padjust_method <- level_headers_info[["pvc_settings"]][["padjust_method"]]
    support <- level_headers_info[["pvc_settings"]][["support"]]
    batch_effects <- level_headers_info[["pvc_settings"]][["batch_effects"]]
    
    # Format alpha threshold smartly
    format_alpha <- function(x) {
        if (x >= 0.0001) {
            # remove trailing zeros
            sub("\\.?0+$", "", format(round(x, 4), nsmall = 0))
        } else {
            format(
                x,
                scientific = TRUE,
                digits = 3
            )
        }
    }
    
    # Compute significance thresholds
    thresholds <- c(
        "*"    = alpha_thresh,
        "**"   = alpha_thresh / 5,
        "***"  = alpha_thresh / 50,
        "****" = alpha_thresh / 500
    )
    
    thresholds_formatted <- vapply(
        thresholds,
        format_alpha,
        character(1)
    )
    alpha_thresh_formatted <- format_alpha(alpha_thresh)
    
    # Generate section header
    section_header <- sprintf(
        "<h2 style='%s' id='section%d'>%s</h2>",
        section_header_style,
        index,
        level
    )
    
    batch_effects_formatted <- if (isFALSE(batch_effects)) {
        "None"
    } else {
        paste(batch_effects, collapse = ", ")
    }
    
    # Generate info block
    info_block <- sprintf(
        "<div style='text-align: center; font-size: 24px; margin-bottom: 20px;'>
       <p><strong>Adjusted p-value threshold:</strong> %s<br>
       <strong>p-adjustment method:</strong> %s<br>
       <strong>Support (min nr non-NA in all timepoints of hit):</strong> %s</p>
       <strong>Batch effects removed:</strong> %s</p>
       <strong>Number of hits:</strong> %s</p>
       <p style='margin-top: 15px; font-size: 22px;'>
         <strong>Significance stars:</strong><br>
         <span>*</span> : p &lt; %s<br>
         <span>**</span> : p &lt; %s<br>
         <span>***</span> : p &lt; %s<br>
         <span>****</span> : p &lt; %s
       </p>
       <hr style='border-top: 1px dashed #999; width: 100%%; margin-top: 15px;'>
     </div>",
        alpha_thresh_formatted,
        padjust_method,
        support,
        batch_effects_formatted,
        nr_of_hits,
        thresholds_formatted["*"],
        thresholds_formatted["**"],
        thresholds_formatted["***"],
        thresholds_formatted["****"]
    )
    
    # Combine and return the full HTML block
    return(paste(
        section_header,
        info_block,
        sep = "\n"
    ))
}