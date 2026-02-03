#' Generate an HTML clustering report from precomputed results
#'
#' Creates an HTML report and associated visualizations for clustering results
#' obtained from a prior SplineOmics analysis. This function performs no
#' statistical modeling or clustering itself; instead, it consumes a
#' precomputed \emph{report payload} and is solely responsible for report
#' generation, plotting, and writing outputs to disk.
#'
#' This strict separation between computation and reporting ensures that
#' analytical results can be inspected, tested, and reused independently of
#' visualization or file-system side effects.
#'
#' @param report_payload `list`: A structured report payload object as returned
#'   by upstream analysis and clustering functions (e.g.,
#'   \code{cluster_hits()}). The payload must contain all data and metadata
#'   required to generate the report, including (non-exhaustive):
#'   \itemize{
#'     \item \code{splineomics}: A `SplineOmics` object containing the original
#'       input data, metadata, design, spline parameters, and analysis settings.
#'     \item \code{all_levels_clustering}: Clustering assignments for all
#'       condition levels and analysis categories.
#'     \item \code{genes}: `character()` Vector of gene or feature names aligned
#'       to the rows of the input data.
#'     \item \code{predicted_timecurves}: Predicted spline trajectories on a
#'       common time grid.
#'     \item \code{category_2_and_3_hits}: Data describing significant features
#'       from average-difference and interaction analyses.
#'     \item Adjusted p-value thresholds and additional metadata required for
#'       report annotation.
#'   }
#'
#' @param plot_info `list`: List with optional elements used to annotate spline
#'   plots:
#'
#'   \itemize{
#'     \item \code{y_axis_label}: `character(1)` Label for the y-axis.
#'     \item \code{time_unit}: `character(1)` Unit used in the x-axis label.
#'     \item \code{treatment_labels}: `list(character(1))` Named list of labels
#'       describing experimental interventions.
#'     \item \code{treatment_timepoints}: `list(numeric(1))` Named list of time
#'       points at which interventions occur.
#'   }
#'
#'   If any treatment list is supplied, both must be present and must share
#'   identical names corresponding to levels of
#'   \code{meta[[condition]]}. Vertical dashed lines are drawn at the specified
#'   time points in the spline plots and annotated with the provided labels.
#'
#' @param plot_options `list`: Named list controlling optional plot
#'   customization. Supported entries include:
#'   \itemize{
#'     \item \code{cluster_heatmap_columns} `logical(1)`
#'       (default = \code{FALSE}):
#'       Whether to cluster columns in the heatmap.
#'
#'     \item \code{meta_replicate_column} `character(1)`
#'       (default = \code{NULL}):
#'       Column in \code{meta} encoding replicate information. When supplied,
#'       spline plot data points are colored by replicate.
#'
#'     \item \code{condition_colours} `list`
#'       (default = \code{NULL}):
#'       Optional named list mapping condition levels to colours. Names must
#'       correspond to values in the condition column of \code{meta}, and
#'       values must be valid R colour specifications accepted by
#'       \code{ggplot2}. When provided, these colours override the default
#'       colours for the corresponding conditions.
#'   }
#'
#' @param verbose `logical(1)`: Logical flag controlling the verbosity of status
#'   and progress messages during report generation.
#'
#' @param max_hit_number `integer(1)`: Maximum number of hits plotted per
#'   cluster.
#'   This parameter can be used to limit report size and computation time when
#'   many significant features are present.
#'
#' @param raw_data `matrix`: Optional matrix containing the raw (unimputed) data
#'   with missing values still present. When supplied, data points that were
#'   originally missing and subsequently imputed are highlighted in the spline
#'   plots.
#'
#' @param report_dir `character(1)`: Path to the directory where the HTML report
#'   and all associated output files will be written. The directory is created
#'   if it does not already exist.
#'
#' @return
#' A named list of plots corresponding to the visualizations included in the
#' generated HTML report. This includes:
#' \itemize{
#'   \item Cluster-specific spline plots for time effects.
#'   \item Heatmaps summarizing cluster assignments.
#'   \item Plots for average-difference and interaction (condition–time)
#'     clustering results, when available.
#' }
#'
#' The returned plots are invisibly written to disk as part of the report
#' generation process.
#'
#' @export
#' 
create_clustering_report <- function(
        report_payload,
        plot_info = list(
            y_axis_label = "Value",
            time_unit = "min",
            treatment_labels = NA,
            treatment_timepoints = NA
        ),
        plot_options = list(
            cluster_heatmap_columns = FALSE,
            meta_replicate_column = NULL,
            condition_colours = NULL
        ),
        verbose = FALSE,
        max_hit_number = 25,
        raw_data = NULL,
        report_dir = here::here()
        ) {
    check_inputs_clustering_report(
        max_hit_number = max_hit_number
    )
    args <- lapply(
        as.list(match.call()[-1]),
        eval,
        parent.frame()
    )
    meta <- splineomics[["meta"]]
    condition <- splineomics[["condition"]]
    args[["meta"]]      <- meta
    args[["condition"]] <- condition
    args[["verbose"]]   <- verbose
    check_null_elements(args)
    input_control <- InputControl$new(args)
    input_control$auto_validate()

    report_dir <- normalizePath(
        report_dir,
        mustWork = FALSE
    )
    
    splineomics <- report_payload[["splineomics"]]
    data <- splineomics[["data"]]
    mode <- splineomics[["mode"]]
    spline_params <- splineomics[["spline_params"]]
    report_info <- splineomics[["report_info"]]
    design <- splineomics[["design"]]
    meta_batch_column <- splineomics[["meta_batch_column"]]
    meta_batch2_column <- splineomics[["meta_batch2_column"]]
    feature_name_columns <- splineomics[["feature_name_columns"]]
    annotation <- splineomics[["annotation"]]
    
    spline_comp_plots <- generate_spline_comparisons(  
        splineomics = splineomics,
        data = data,
        meta = meta,
        condition = condition,
        replicate_column = plot_options[["meta_replicate_column"]],
        plot_info = plot_info,
        plot_options = plot_options,
        category_2_and_3_hits = report_payload[["category_2_and_3_hits"]],
        adj_pthresh_avrg_diff_conditions = 
            report_payload[["adj_pthresh_avrg_diff_conditions"]],
        adj_pthresh_interaction = 
            report_payload[["adj_pthresh_interaction_condition_time"]],
        min_effect_size = report_payload[["min_effect_size"]],
        raw_data = raw_data,
        predicted_timecurves = report_payload[["predicted_timecurves"]],
        max_hit_number = max_hit_number
    )
    
    # Put them in there under those names, so that the report generation fun
    # can access them directly like this.
    effects <- extract_effects(report_payload[["design"]])
    report_info[["Fixed effects"]] <- effects[["fixed_effects"]]
    report_info[["Random effects"]] <- effects[["random_effects"]]
    report_info[["meta_condition"]] <- c(condition)
    report_info[["plot_data_batch_correction"]] <- paste(
        meta_batch_column,
        meta_batch2_column,
        sep = ", "
    )
    report_info[["homosc_violation_result"]] <-
        splineomics[["homosc_violation_result"]]
    report_info[["min_effect_size"]] <- report_payload[["min_effect_size"]]
    
    if (!is.null(splineomics[["use_array_weights"]])) {
        report_info[["use_array_weights"]] <- splineomics[["use_array_weights"]]
        report_info[["heteroscedasticity"]] <- "not tested"
    } else {
        report_info[["use_array_weights"]] <- paste(
            "automatic (decided by Levene's test), array_weights only used ",
            "when heteroscedasticity is detected (% violating features >= 10)"
        )
        report_info[["heteroscedasticity"]] <- sprintf(
            "Heteroscedasticity detected: %s (%.1f%% of features violated the
      assumption of homoscedasticity)",
            ifelse(
                splineomics[["homosc_violation_result"]][["violation"]],
                "Yes",
                "No"
            ),
            splineomics[["homosc_violation_result"]][["percent_violated"]]
        )
    }
    
    plots <- make_clustering_report(
        all_levels_clustering = report_payload[["all_levels_clustering"]],
        condition = condition,
        data = data,
        meta = meta,
        annotation = annotation,
        genes = report_payload[["genes"]],
        spline_params = spline_params,
        adj_pthresh_time_effect = report_payload[["adj_pthresh_time_effect"]],
        adj_pthresh_avrg_diff_conditions =
            report_payload[["adj_pthresh_avrg_diff_conditions"]],
        adj_pthresh_interaction_condition_time =
            report_payload[["adj_pthresh_interaction_condition_time"]],
        category_2_and_3_hits = report_payload[["category_2_and_3_hits"]],
        report_dir = report_dir,
        mode = mode,
        report_info = report_info,
        predicted_timecurves = report_payload[["predicted_timecurves"]],
        design = design,
        meta_batch_column = meta_batch_column,
        meta_batch2_column = meta_batch2_column,
        plot_info = plot_info,
        plot_options = plot_options,
        feature_name_columns = feature_name_columns,
        spline_comp_plots = spline_comp_plots,
        raw_data = raw_data,
        max_hit_number = max_hit_number,
        verbose = verbose
    )
    print_info_message(
        message_prefix = "Clustering the hits",
        report_dir = report_dir
    )
    plots
}


# Level 1 internal functions ---------------------------------------------------


check_inputs_clustering_report <- function(
        max_hit_number
        ) {
    # validate max_hit_number
    if (!is.numeric(max_hit_number) ||
        length(max_hit_number) != 1L ||
        !(is.infinite(max_hit_number) ||
          (max_hit_number >= 1 &&
           max_hit_number == as.integer(max_hit_number)))) {
        stop_call_false(
            "`max_hit_number` must be a single positive integer or Inf."
        )
    }
    
    TRUE
}

#' Make Clustering Report
#'
#' @noRd
#'
#' @description
#' Generates a detailed clustering report including heatmaps, dendrograms,
#' curve plots, and consensus shapes for each level within a condition.
#'
#' @param splineomics A list containing the splineomics results, including
#' time effects, average difference between conditions, and interaction between
#' condition and time.
#' @param all_levels_clustering A list containing clustering results for each
#' level within a condition.
#' @param condition A character string specifying the condition.
#' @param data A matrix of data values.
#' @param meta A dataframe containing metadata.
#' @param annotation Dataframe containig the annotation info of the features,
#'                   such as gene and uniprotID, for example.
#' @param genes Character vector containing the genes of the features.
#' @param spline_params A list of spline parameters for the analysis.
#' @param adj_pthresh_time_effect `numeric(1)`: adj. p-value threshold
#' for the limma time effect results (category 1).
#' @param adj_pthresh_avrg_diff_conditions adj. p-value threshold
#' for the limma average difference between conditions results (category 2).
#' @param adj_pthresh_interaction_condition_time adj. p-value threshold
#' for the limma interaction condition time results (category 3).
#' @param category_2_and_3_hits List of dataframes, where each df is the part
#' of the toptable that contains the significant features of the respective
#' limma result category (2 or 3).
#' @param report_dir A character string specifying the report directory.
#' @param mode A character string specifying the mode
#' ('isolated' or 'integrated').
#' @param report_info A named list containing report information such as analyst
#'                    name, fixed and random effects, etc.
#' @param predicted_timecurves A list containing:
#'   \describe{
#'     \item{`time_grid`}{A numeric vector of dense time points used for
#'       evaluation.}
#'     \item{`predictions`}{A named list of matrices, one per condition level.
#'       Each matrix contains predicted values (rows = features, columns =
#'       timepoints).}
#'   }
#' @param fit Full fitted model returned by limma or variancePartition::dream.
#' @param design A string representing the limma design formula
#' @param meta_batch_column A character string specifying the meta batch column.
#' @param meta_batch2_column A character string specifying the second meta
#'                           batch column.
#' @param plot_info List containing the elements y_axis_label (string),
#'                  time_unit (string), treatment_labels (character vector),
#'                  treatment_timepoints (integer vector). All can also be NA.
#'                  This list is used to add this info to the spline plots.
#'                  time_unit is used to label the x-axis, and treatment_labels
#'                  and -timepoints are used to create vertical dashed lines,
#'                  indicating the positions of the treatments (such as
#'                  feeding, temperature shift, etc.).
#' @param plot_options `list`: Named list controlling optional plot
#'   customization. Supported entries include:
#'   \itemize{
#'     \item \code{cluster_heatmap_columns} `logical(1)`
#'       (default = \code{FALSE}):
#'       Whether to cluster columns in the heatmap.
#'
#'     \item \code{meta_replicate_column} `character(1)`
#'       (default = \code{NULL}):
#'       Column in \code{meta} encoding replicate information. When supplied,
#'       spline plot data points are colored by replicate.
#'
#'     \item \code{condition_colours} `list`
#'       (default = \code{NULL}):
#'       Optional named list mapping condition levels to colours. Names must
#'       correspond to values in the condition column of \code{meta}, and
#'       values must be valid R colour specifications accepted by
#'       \code{ggplot2}. When provided, these colours override the default
#'       colours for the corresponding conditions.
#'   }
#' @param feature_name_columns Character vector containing the column names of
#'                             the annotation info that describe the features.
#'                             This argument is used to specify in the HTML
#'                             report how exactly the feature names displayed
#'                             above each individual spline plot have been
#'                             created. Use the same vector that was used to
#'                             create the row headers for the data matrix!
#' @param spline_comp_plots List containing the list of lists with all
#' the plots for all the pairwise comparisons of the condition in terms of
#' average spline diff and interaction condition time, and another list of lists
#' where the respective names of each plot are stored.
#' @param raw_data Optional. Data matrix with the raw (unimputed) data, still
#' containing NA values. When provided, it highlights the datapoints in the
#' spline plots that originally where NA and that were imputed.
#' @param max_hit_number Maximum number of hits which are plotted within each
#' cluster. This can be used to limit the computation time and size of
#' the HTML report in the case of many hits.
#' @param verbose Boolean flag controlling the display of messages.
#'
#' @return No return value, called for side effects.
#'
#' @seealso
#' \code{\link{removeBatchEffect}}, \code{\link{plot_heatmap}},
#' \code{\link{plot_cluster_mean_splines}}, \code{\link{plot_splines}},
#' \code{\link{generate_report_html}}
#'
#' @importFrom limma removeBatchEffect
#' @importFrom dplyr filter
#' @importFrom stats na.omit
#' @importFrom rlang .data
#'
make_clustering_report <- function(
        all_levels_clustering,
        condition,
        data,
        meta,
        annotation,
        genes,
        spline_params,
        adj_pthresh_time_effect,
        adj_pthresh_avrg_diff_conditions,
        adj_pthresh_interaction_condition_time,
        category_2_and_3_hits,
        report_dir,
        mode,
        report_info,
        predicted_timecurves,
        design,
        meta_batch_column,
        meta_batch2_column,
        plot_info,
        plot_options,
        feature_name_columns,
        spline_comp_plots,
        raw_data,
        max_hit_number,
        verbose
        ) {
    design <- gsub("Time", "X", design)
    effects <- extract_effects(design)
    datas <- split_data_by_condition(
        data = data,
        meta = meta,
        condition = condition,
        mode = mode
    )
    
    # To extract the stored value for the potential auto cluster decision.
    # collect k values only from valid list entries
    clusters <- integer(0)
    # Normalize: replace string placeholders with NULL
    all_levels_clustering <- lapply(all_levels_clustering, function(x) {
        if (is.character(x)) {
            return(NA)
        }
        x
    })
    
    for (i in seq_along(all_levels_clustering)) {
        x <- all_levels_clustering[[i]]
        
        # skip non-lists (e.g., "No result ...") and NULL/NA
        if (!is.list(x) || is.null(x) || all(is.na(x))) next
        
        # collect k if present and valid; then remove the field
        if ("clusters" %in% names(x)) {
            k <- x$clusters
            if (is.numeric(k) && length(k) == 1L && !is.na(k)) {
                clusters <- c(clusters, as.integer(k))
            }
            all_levels_clustering[[i]]$clusters <- NULL
        }
    }
    
    # ensure the report dir exists
    if (!dir.exists(report_dir)) {
        dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    time_unit_label <- paste0("[", plot_info$time_unit, "]")
    
    if (isTRUE(verbose)) message("Generating heatmap...")
    heatmaps <- plot_heatmap(
        datas = datas,
        meta = meta,
        mode = mode,
        condition = condition,
        all_levels_clustering = all_levels_clustering,
        time_unit_label = time_unit_label,
        cluster_heatmap_columns = plot_options[["cluster_heatmap_columns"]],
        max_hit_number = max_hit_number
    )
    
    cluster_quality_plots <- lapply(
        all_levels_clustering,
        plot_cluster_quality
    )
    
    level_headers_info <- list()
    plots <- list()
    plots_sizes <- list()
    q <- 0
    levels <- sub(
        paste0("^", condition, "_"),
        "",
        names(all_levels_clustering)
    )
    
    for (i in seq_along(all_levels_clustering)) {
        # When a level has < 2 hits
        if (is.null(all_levels_clustering[[i]]) ||
            all(is.na(all_levels_clustering[[i]]))) {
            next
        } else {
            q <- q + 1
        }
        level_clustering <- all_levels_clustering[[i]]
        
        if (length(levels) >= i) {
            level <- as.character(levels[i])
            
            # Get indices of columns in meta that match the given level
            # (condition)
            condition_indices <- which(meta[[condition]] == level)
            
            # Subset raw_data to only include these columns (keeping all rows)
            raw_data_level <- raw_data[, condition_indices, drop = FALSE]
            
            # Construct header name
            header_name <- level
            
            nr_hits <- nrow(level_clustering$clustered_hits)
            
            header_info <- list(
                header_name = header_name,
                nr_hits = nr_hits,
                adj_pvalue_threshold = adj_pthresh_time_effect
            )
            
            level_headers_info[[i]] <- header_info
        }
        
        curve_values <- level_clustering$curve_values
        
        p_curves <- plot_all_mean_splines(
            curve_values = curve_values,
            plot_info = plot_info,
            level = level
        )
        
        if (verbose) {
            message(paste("Generating cluster mean splines for level: ", level))
        }
        cluster_mean_splines <- plot_cluster_mean_splines(
            curve_values = curve_values,
            plot_info = plot_info,
            level = level
        )
        
        top_table <- level_clustering$top_table
        col_indices <- which(meta[[condition]] == levels[i])
        
        if (mode == "integrated") {
            data_level <- datas[[i]][, col_indices]
        } else { # mode == "isolated"
            data_level <- datas[[i]]
        }
        
        meta_level <- meta |> dplyr::filter(.data[[condition]] == levels[i])
        
        clusters_spline_plots <- list()
        
        if (isTRUE(verbose)) message("Generating spline plots...")
        for (nr_cluster in sort(unique(stats::na.omit(top_table$cluster)))) {
            nr_of_hits <- sum(
                level_clustering$clustered_hits$cluster == nr_cluster,
                na.rm = TRUE
            )
            main_title <- paste(
                "Cluster",
                nr_cluster,
                " | Hits:",
                nr_of_hits,
                sep = " "
            )
            
            top_table_cluster <- top_table |>
                dplyr::filter(!!rlang::sym("cluster") == nr_cluster)
            
            X <- level_clustering$X
            
            spline_plots <- plot_splines(
                top_table = top_table_cluster,
                data = data_level,
                meta = meta_level,
                predicted_timecurves = predicted_timecurves,
                time_unit_label = time_unit_label,
                plot_info = plot_info,
                plot_options = plot_options,
                adj_pthreshold = adj_pthresh_time_effect,
                replicate_column = plot_options[["meta_replicate_column"]],
                level = level,
                raw_data = raw_data_level,
                report_info = report_info,
                max_hit_number = max_hit_number,
                all_levels_clustering = all_levels_clustering,
                condition = condition
            )
            
            clusters_spline_plots[[length(clusters_spline_plots) + 1]] <- list(
                spline_plots = spline_plots,
                cluster_main_title = main_title
            )
        }
        
        plots <- c(
            plots,
            new_level = "level_header", # is the signal for the plotting code
            p_curves = list(p_curves),
            cluster_mean_splines = list(cluster_mean_splines),
            cluster_quality_plots = list(cluster_quality_plots[[i]]),
            heatmap = heatmaps[[i]],
            individual_spline_plots = clusters_spline_plots 
        )
        
        # For every plot in plots, this determines the size in the HTML
        plots_sizes <- c(
            plots_sizes,
            999, # dummy size for "next_level" signal
            1.5,
            1,
            1.5,
            1.5,
            rep(1, length(clusters_spline_plots))
        )
    }
    
    topTables <- list()
    
    # Loop over each element in all_levels_clustering
    for (i in seq_along(all_levels_clustering)) {
        if (is.logical(all_levels_clustering[[i]])) next
        
        # Get the current element, which is a list
        current_element <- all_levels_clustering[[i]]
        
        # Extract the top_table element
        top_table_element <- current_element$top_table
        
        # Get the name of the outer list element
        element_name <- names(all_levels_clustering)[i]
        
        # Trim the name to 30 characters if necessary
        if (nchar(element_name) > 30) {
            element_name <- substr(element_name, 1, 30)
        }
        
        topTables[[element_name]] <- top_table_element
    }
    
    topTables <- transfer_sr2cc(
        topTables = topTables,
        all_levels_clustering = all_levels_clustering
    )
    
    if (!is.null(genes)) {
        enrichr_format <- prepare_gene_lists_for_enrichr(
            all_levels_clustering,
            genes
        )
    } else {
        enrichr_format <- NA
    }
    
    all_levels_clustering <- merge_annotation_all_levels_clustering(
        all_levels_clustering = all_levels_clustering,
        annotation = annotation
    )
    
    if (isTRUE(verbose)) message("Generating report. This takes a few seconds.")
    report_info[["max_hit_number"]] <- max_hit_number
    
    generate_report_html(
        plots = plots,
        limma_result_2_and_3_plots = spline_comp_plots,
        plots_sizes = plots_sizes,
        level_headers_info = level_headers_info,
        spline_params = spline_params,
        report_info = report_info,
        data = bind_data_with_annotation(data, annotation),
        meta = meta,
        topTables = topTables,
        category_2_and_3_hits = category_2_and_3_hits,
        enrichr_format = enrichr_format,
        adj_pthresh_time_effect = adj_pthresh_time_effect,
        adj_pthresh_avrg_diff_conditions = adj_pthresh_avrg_diff_conditions,
        adj_pthresh_interaction_condition_time =
            adj_pthresh_interaction_condition_time,
        report_type = "cluster_hits",
        feature_name_columns = feature_name_columns,
        mode = mode,
        filename = "report_clustered_hits",
        report_dir = report_dir
    )
    
    return(plots)
}


# Level 2 internal functions ---------------------------------------------------


#' Prepare Gene Lists for Enrichr and Return as String
#'
#' @noRd
#'
#' @description
#' This function processes the clustered hits in each element of
#' `all_levels_clustering`, formats the gene names for easy copy-pasting into
#' Enrichr, and returns the formatted gene lists as a string.
#'
#' @param all_levels_clustering A list where each element contains a dataframe
#' `clustered_hits` with columns `feature` and `cluster`.
#' @param genes A vector of gene names corresponding to the feature indices.
#'
#' @return A character vector with the formatted gene lists for each cluster.
#'
prepare_gene_lists_for_enrichr <- function(
        all_levels_clustering,
        genes) {
    formatted_gene_lists <- list()
    
    for (i in seq_along(all_levels_clustering)) {
        if (is.logical(all_levels_clustering[[i]])) next
        
        level_name <- names(all_levels_clustering)[i]
        clustered_hits <- all_levels_clustering[[i]]$clustered_hits
        
        # Process each cluster
        clusters <- split(
            clustered_hits$feature,
            clustered_hits$cluster
        )
        
        level_gene_lists <- list()
        
        for (cluster_id in names(clusters)) {
            cluster_genes <- clusters[[cluster_id]]
            
            gene_list <- genes[cluster_genes]
            gene_list <- na.omit(gene_list) # Remove NAs if any
            
            if (length(gene_list) > 0) {
                level_gene_lists[[paste0("Cluster ", cluster_id)]] <-
                    paste(gene_list, collapse = "\n")
            }
        }
        
        formatted_gene_lists[[level_name]] <- level_gene_lists
    }
    
    # Prepare the background genes list using preprocessed genes
    background_gene_list <- paste(
        na.omit(genes),
        collapse = "\n"
    )
    
    return(list(
        gene_lists = formatted_gene_lists,
        background = background_gene_list
    ))
}


#' Split or duplicate a data matrix by condition level
#'
#' @noRd
#'
#' @description
#' This function returns a list of data matrices, either split by
#' condition level (if mode is "isolated") or duplicated for each level
#' (if mode is "integrated"). This is used to prepare per-condition data
#' subsets for downstream processing or plotting.
#'
#' @param data A numeric matrix with features as rows and samples as
#'   columns.
#' @param meta A data.frame containing sample-level metadata. Each row
#'   corresponds to a column in the data matrix.
#' @param condition A character string specifying the column in `meta`
#'   that defines the condition or group variable.
#' @param mode A string, either "isolated" or "integrated", determining
#'   whether to split the data or return full copies per level.
#'
#' @return A named list of data matrices, one for each condition level.
#'
split_data_by_condition <- function(
        data,
        meta,
        condition,
        mode) {
    datas <- list()
    levels <- unique(meta[[condition]])
    
    for (level in levels) {
        if (mode == "isolated") {
            cols <- which(meta[[condition]] == level)
            datas[[level]] <- data[, cols, drop = FALSE]
        } else {
            datas[[level]] <- data # same full data for all levels
        }
    }
    
    return(datas)
}


#' Merge Annotation with All Top Tables
#'
#' @noRd
#'
#' @description
#' This function merges annotation information into the `top_table` of each
#' non-logical element in a list.
#'
#' @param all_levels_clustering A list where each element contains a `top_table`
#' dataframe with a `feature_nr` column. Some elements may be logical values.
#' @param annotation A dataframe containing the annotation information.
#'
#' @return A list with updated `top_table` dataframes containing merged
#' annotation information.
#'
merge_annotation_all_levels_clustering <- function(
        all_levels_clustering,
        annotation = NULL) {
    all_levels_clustering <- lapply(
        all_levels_clustering,
        function(x) {
            # Check if x is not logical and annotation is not NULL
            if (!is.logical(x) && !is.null(annotation)) {
                x$top_table <- merge_top_table_with_annotation(
                    x$top_table,
                    annotation
                )
            }
            return(x)
        }
    )
    
    return(all_levels_clustering)
}


#' Transfer per-member cluster quality scores into top tables
#'
#' @noRd
#'
#' @description
#' For each condition level present in both `topTables` and
#' `all_levels_clustering`, transfer the per-feature values from
#' `cluster_quality$per_member` (in `all_levels_clustering`) into the
#' corresponding tibble in `topTables`.
#'
#' The mapping is done by matching feature identifiers:
#' - `feature` in `clustered_hits` is aligned positionally with entries
#'   in `per_member` (row i of `clustered_hits` corresponds to element i
#'   of `per_member`).
#' - `feature_nr` in each tibble of `topTables` is then matched against
#'   `feature` in `clustered_hits` to retrieve the correct per-member
#'   score.
#'
#' A new column `sr2cc` is added to each tibble in `topTables`:
#' - Values are numeric, taken from the corresponding `per_member`
#'   element after ID-based matching.
#' - Unmatched rows are assigned `NA`.
#' - Any existing `sr2cc` column is removed and replaced.
#'
#' @param topTables
#'   A named list of tibbles. Each tibble must contain a column
#'   `feature_nr` with numeric feature IDs used for matching.
#' @param all_levels_clustering
#'   A named list of condition-level results. Each element must contain:
#'   - `cluster_quality$per_member`: numeric vector of per-feature
#'     scores (aligned to rows of `clustered_hits`).
#'   - `clustered_hits`: data frame with a numeric `feature` column.
#'
#' @return
#'   A list of tibbles (same structure as `topTables`), each with a new
#'   column `sr2cc` holding per-member cluster quality scores or `NA`
#'   if no match was found.
#'
transfer_sr2cc <- function(
        topTables,
        all_levels_clustering
) {
    out <- topTables
    lvl_names <- intersect(names(all_levels_clustering), names(out))
    if (length(lvl_names) == 0L) {
        return(out)
    }
    
    for (lvl in lvl_names) {
        alc <- all_levels_clustering[[lvl]]
        tt <- out[[lvl]]
        
        if (is.null(alc$cluster_quality$per_member) ||
            is.null(alc$clustered_hits) ||
            !is.data.frame(alc$clustered_hits) ||
            is.null(tt) || !is.data.frame(tt) ||
            !("feature_nr" %in% names(tt))) {
            next
        }
        
        per_member <- alc$cluster_quality$per_member
        hits_df <- alc$clustered_hits
        
        # Coerce IDs to numeric for robust matching
        f_ids <- as.numeric(hits_df$feature)
        t_ids <- as.numeric(tt$feature_nr)
        
        # Match feature_nr (tt) to feature (hits_df)
        idx <- match(t_ids, f_ids)
        
        # Pull per_member by the matched row position in hits_df
        sr2cc <- rep(NA_real_, nrow(tt))
        ok <- !is.na(idx) & idx >= 1 & idx <= length(per_member)
        sr2cc[ok] <- as.numeric(per_member[idx[ok]])
        
        tt$sr2cc <- sr2cc
        out[[lvl]] <- tt
    }
    
    out
}


#' Plot per-cluster signed r^2 quality distributions (one category block)
#'
#' @noRd
#'
#' @description
#' For a single result/category component, generate a density-style histogram
#' of member \emph{signed} r\eqn{^2} (variance explained with the sign of the
#' correlation) for each cluster. Clusters without valid quality values are
#' skipped. The returned list is named and sorted as \code{cluster_1},
#' \code{cluster_2}, ... by numeric cluster id.
#'
#' @param category_result A list-like object for one result/category component.
#'   Must contain:
#'   \describe{
#'     \item{\code{cluster_quality}}{A list with \code{per_member}, the vector
#'       of signed r\eqn{^2} scores aligned to rows in \code{clustered_hits}.}
#'     \item{\code{clustered_hits}}{A data.frame with a \code{cluster} column
#'       giving the cluster id for each row/feature.}
#'   }
#'
#' @return A named list of \code{ggplot} objects (one per cluster), sorted by
#'   increasing cluster id. Returns \code{NULL} if no valid plots can be made.
#'
#' @details
#' This expects that \code{cluster_quality$per_member} came from
#' \code{compute_cluster_fits()} (signed r\eqn{^2} to centroid) and that its
#' order matches \code{clustered_hits}.
#'
plot_cluster_quality <- function(category_result) {
    # Basic structure checks
    if (!is.list(category_result)) {
        return(NULL)
    }
    cq <- category_result$cluster_quality
    ch <- category_result$clustered_hits
    if (is.null(cq) || is.null(cq$per_member) || is.null(ch)) {
        return(NULL)
    }
    if (!("cluster" %in% names(ch))) {
        return(NULL)
    }
    
    # Collect cluster ids, keep numeric/finite only, sort
    cl_ids_num <- sort(as.integer(unique(ch$cluster)))
    cl_ids_num <- cl_ids_num[is.finite(cl_ids_num)]
    if (length(cl_ids_num) == 0) {
        return(NULL)
    }
    
    # Build plots in sorted order; drop empty/NA-only
    plots <- lapply(cl_ids_num, function(clid) {
        idx <- which(as.numeric(ch$cluster) == clid)
        sr2 <- cq$per_member[idx] # signed r^2 per member
        if (length(sr2) == 0 || all(is.na(sr2))) {
            return(NULL)
        }
        plot_cluster_quality_distribution(sr2, clid)
    })
    
    # Drop NULLs
    keep <- !vapply(plots, is.null, logical(1))
    plots <- plots[keep]
    if (length(plots) == 0) {
        return(NULL)
    }
    
    # Name by sorted cluster id (cluster_1, cluster_2, ...)
    names(plots) <- paste0("cluster_", cl_ids_num[keep])
    plots
}


#' Plot All Mean Splines
#'
#' @noRd
#'
#' @description
#' Generates a plot of average curves for each cluster, showing z-score
#' normalized intensities over time.
#'
#' @param curve_values A dataframe containing curve values and cluster
#' assignments.
#' @param plot_info List containing the elements y_axis_label (string),
#'                  time_unit (string), treatment_labels (character vector),
#'                  treatment_timepoints (integer vector). All can also be NA.
#'                  This list is used to add this info to the spline plots.
#'                  time_unit is used to label the x-axis, and treatment_labels
#'                  and -timepoints are used to create vertical dashed lines,
#'                  indicating the positions of the treatments (such as
#'                  feeding, temperature shift, etc.).
#' @param level One of the unique values of the meta condition column. This is
#'              a factor that separates the experiment.
#'
#' @return A ggplot object representing the average curves by cluster.
#'
#' @importFrom scales hue_pal
#' @importFrom ggplot2 ggplot geom_line ggtitle xlab ylab scale_color_brewer
#'                     theme_minimal aes element_text
#' @importFrom rlang .data
#'
plot_all_mean_splines <- function(
        curve_values,
        plot_info,
        level
        ) {
    time <- as.numeric(colnames(curve_values)[-length(colnames(curve_values))])
    
    clusters <- unique(curve_values$cluster)
    average_curves <- data.frame()
    
    # Loop through each unique cluster value to calculate the average curve
    for (current_cluster in clusters) {
        # Filter rows for the current cluster
        subset_hits <- curve_values[curve_values$cluster == current_cluster, ]
        last_timepoint <- (which(names(curve_values) == "cluster")) - 1
        average_curve <- colMeans(subset_hits[, seq_len(last_timepoint)])
        
        # Create a data frame for the average curve with an additional 'Cluster'
        # column
        curve_df <- data.frame(
            Time = time, Value = average_curve,
            cluster = as.factor(current_cluster)
        )
        
        # Bind the curve data frame to the cumulative data frame
        average_curves <- rbind(
            average_curves,
            curve_df
        )
    }
    
    average_curves$cluster <- factor(
        average_curves$cluster,
        levels = sort(
            unique(as.numeric(average_curves$cluster))
        )
    )
    
    time_unit_label <- paste0("[", plot_info$time_unit, "]")
    
    cluster_colors <- get_cluster_colors(curve_values)
    
    if (length(cluster_colors) > length(unique(average_curves$cluster))) {
        cluster_colors <-
            cluster_colors[seq_len(length(unique(average_curves$cluster)))]
    }
    names(cluster_colors) <- paste(
        "Cluster",
        levels(average_curves$cluster)
    )
    
    color_values <- c(cluster_colors)
    distinct_colors <- c()
    
    # Create the base plot
    p_curves <- ggplot2::ggplot(
        average_curves,
        ggplot2::aes(
            x = !!rlang::sym("Time"),
            y = !!rlang::sym("Value"),
            color = paste("Cluster", factor(!!rlang::sym("cluster")))
        )
    ) +
        ggplot2::geom_line() +
        ggplot2::ggtitle(
            sprintf("Cluster Centroid (average spline) - %s", level)
        ) +
        ggplot2::xlab(paste("Time", time_unit_label)) +
        ggplot2::ylab(paste("z-score norm.", plot_info$y_axis_label)) +
        ggplot2::theme_minimal() +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    
    # Call the wrapper function to conditionally add dashed lines and get
    # treatment colors
    result <- maybe_add_dashed_lines(
        p = p_curves,
        plot_info = plot_info,
        level = level
    )
    
    p_curves <- result$p
    treatment_colors <- result$treatment_colors
    
    # Combine cluster colors and treatment colors for a single color scale
    all_colors <- c(cluster_colors, treatment_colors)
    
    # Finalize color scale and theme adjustments
    p_curves <- p_curves +
        ggplot2::scale_color_manual(
            values = all_colors, # Combine both cluster and treatment colors
            name = NULL # No legend title
        ) +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(size = 12),
            axis.text.y = ggplot2::element_text(size = 12),
            axis.title.x = ggplot2::element_text(size = 14),
            axis.title.y = ggplot2::element_text(size = 14),
            legend.text = ggplot2::element_text(size = 12),
            legend.key.size = grid::unit(0.9, "cm"),
            legend.key.height = grid::unit(0.6, "cm"),
            plot.title = ggplot2::element_text(size = 16)
        )
    
    
    return(p_curves)
}


#' Generate spline comparison plots for all condition pairs
#'
#' @noRd
#'
#' @description
#' Generates the "double spline plots" (limma result categories 2 & 3)
#' for all pairwise combinations of condition levels in the metadata.
#' For each condition pair, the function:
#' \itemize{
#'   \item extracts the spline-based time effects for both conditions,
#'   \item subsets the corresponding average-difference (Category 2) and
#'     interaction (Category 3) top tables for that pair,
#'   \item calls \code{plot_spline_comparisons()} to plot data points and
#'     overlay fitted spline curves.
#' }
#' Plots are only generated for features that pass the adjusted p-value
#' thresholds and additional effect-size criteria applied upstream in
#' \code{get_category_2_and_3_hits()}.
#'
#' @param splineomics A list containing the splineomics results, including
#'   per-condition time effects and pairwise contrast tables produced by
#'   the spline-based limma/dream workflow.
#'   
#' @param data The data matrix containing the measurements.
#' 
#' @param meta The metadata associated with the measurements, which includes
#'   the condition column.
#'   
#' @param condition Column name of \code{meta} that contains the experimental
#'   condition levels.
#'   
#' @param replicate_column Column name of the \code{meta} column that
#'   specifies the replicates per timepoint. For example \code{Reactor} with
#'   values like \code{"ReactorE16"}, \code{"ReactorE17"}, etc., indicating
#'   multiple bioreactors contributing samples at each timepoint.
#'   
#' @param plot_info A list containing plotting information such as time unit
#'   and axis labels.
#'   
#' @param plot_options `list`: Named list controlling optional plot
#'   customization. Supported entries include:
#'   \itemize{
#'     \item \code{cluster_heatmap_columns} `logical(1)`
#'       (default = \code{FALSE}):
#'       Whether to cluster columns in the heatmap.
#'
#'     \item \code{meta_replicate_column} `character(1)`
#'       (default = \code{NULL}):
#'       Column in \code{meta} encoding replicate information. When supplied,
#'       spline plot data points are colored by replicate.
#'
#'     \item \code{condition_colours} `list`
#'       (default = \code{NULL}):
#'       Optional named list mapping condition levels to colours. Names must
#'       correspond to values in the condition column of \code{meta}, and
#'       values must be valid R colour specifications accepted by
#'       \code{ggplot2}. When provided, these colours override the default
#'       colours for the corresponding conditions.
#'   }
#'   
#' @param raw_data Optional data matrix with the raw (unimputed) data,
#'   still containing \code{NA} values. When provided, datapoints that were
#'   originally \code{NA} and later imputed can be highlighted in the plots.
#'   
#' @param predicted_timecurves A list containing:
#'   \describe{
#'     \item{time_grid}{Numeric vector of dense time points used for
#'       evaluation.}
#'     \item{predictions}{Named list of matrices, one per condition level.
#'       Each matrix contains predicted values (rows = features, columns =
#'       time points).}
#'   }
#'   
#' @param max_hit_number Maximum number of hits per condition pair for which
#'   individual spline plots are shown. This can be used to limit computation
#'   time and the size of the HTML report when many hits are present.
#'   
#' @param category_2_and_3_hits List returned by
#'   \code{get_category_2_and_3_hits()}, containing per-contrast Category 2
#'   and Category 3 hit tables (one element per pairwise comparison).
#'   
#' @param adj_pthresh_avrg_diff_conditions Adjusted p-value threshold for
#'   the average difference between conditions (Category 2). Passed through
#'   to plotting helpers for annotation or additional filtering.
#'   
#' @param adj_pthresh_interaction Adjusted p-value threshold for the
#'   condition × time interaction (Category 3). Passed through to plotting
#'   helpers for annotation or additional filtering.
#'   
#' @param min_effect_size `list`: A named list that specifies the minimum 
#' effect size thresholds to consider a feature as biologically meaningful, in 
#' addition to statistical significance. This allows users to filter out 
#' "trivial" hits that pass adjusted p-value cutoffs but show negligible effect 
#' sizes.
#'
#'   The list must contain the following elements:
#'   - `time_effect`: `numeric(1)` Minimum cumulative travel for time effects 
#'     (Category 1). Features with a smaller travel will be ignored even if 
#'     significant.
#'   - `avg_diff_cond`: `numeric(1)` Minimum absolute effect size for average 
#'     differences between conditions (Category 2). Ensures that only contrasts 
#'     with a relevant magnitude are reported.
#'   - `interaction_cond_time`: `numeric(1)` Minimum effect size for the 
#'     interaction between condition and time (Category 3). This controls how 
#'     large the differential curve travel must be across conditions to count 
#'     as a hit.
#'
#'   Values should be numeric scalars (typically >0). For example: 
#'   `min_effect_size = list(time_effect = 1, avg_diff_cond = 1, 
#'   interaction_cond_time = 2)` will only keep features with cumulative 
#'   travels or condition-time differences above those cutoffs. Use smaller 
#'   values (e.g., 0.1) for permissive filtering, or larger values for more 
#'   conservative thresholds.
#'
#'   The default is 0 for all three elements.
#'
#' @return A named list of lists. Each top-level element corresponds to a
#'   condition pair (e.g. \code{"A_vs_B"}) and contains the comparison plots
#'   and associated feature names for that pair.
#'
generate_spline_comparisons <- function(
        splineomics,
        data,
        meta,
        condition,
        replicate_column,
        plot_info,
        plot_options,
        raw_data,
        predicted_timecurves,
        max_hit_number,
        category_2_and_3_hits,
        adj_pthresh_avrg_diff_conditions,
        adj_pthresh_interaction,
        min_effect_size
) {
    # Ensure `condition` column is character
    meta[[condition]] <- as.character(meta[[condition]])
    
    limma_res <- splineomics[["limma_splines_result"]]
    te_list <- limma_res[["time_effect"]]
    avrg_list <- limma_res[["avrg_diff_conditions"]]
    inter_list <- limma_res[["interaction_condition_time"]]
    
    if (is.null(avrg_list) && is.null(inter_list)) {
        return(list())
    }
    
    # Derive contrast suffixes like "A_vs_B" from list names
    suffix_from_avrg <- if (!is.null(avrg_list)) {
        sub("^avrg_diff_", "", names(avrg_list))
    } else {
        character(0)
    }
    
    suffix_from_inter <- if (!is.null(inter_list)) {
        sub("^time_interaction_", "", names(inter_list))
    } else {
        character(0)
    }
    
    suffixes <- sort(unique(c(suffix_from_avrg, suffix_from_inter)))
    suffixes <- suffixes[nzchar(suffixes)]
    
    if (length(suffixes) == 0L) {
        return(list())
    }
    
    comparison_plots <- vector("list", length(suffixes))
    names(comparison_plots) <- suffixes
    
    empty_tbl <- tibble::tibble(
        feature_nr = numeric(0),
        feature_names = character(0)
    )
    
    for (suffix in suffixes) {
        # suffix is e.g. "A_vs_B"
        parts <- strsplit(suffix, "_vs_")[[1L]]
        if (length(parts) != 2L) {
            next
        }
        c1 <- parts[1L]
        c2 <- parts[2L]
        
        # time_effect names like "condition_A", "condition_B", ...
        te_name_1 <- paste0(condition, "_", c1)
        te_name_2 <- paste0(condition, "_", c2)
        time_effect_1 <- te_list[[te_name_1]]
        time_effect_2 <- te_list[[te_name_2]]
        
        # per-contrast limma tables
        avrg_name <- paste0("avrg_diff_", suffix)
        inter_name <- paste0("time_interaction_", suffix)
        
        avrg_diff_df <- if (!is.null(avrg_list) &&
                            avrg_name %in% names(avrg_list)) {
            
            avrg_list[[avrg_name]]
        } else {
            empty_tbl
        }
        
        interaction_df <- if (!is.null(inter_list) &&
                              inter_name %in% names(inter_list)) {
            
            inter_list[[inter_name]]
        } else {
            empty_tbl
        }
        
        plots_and_feature_names <- plot_spline_comparisons(
            time_effect_1 = time_effect_1,
            condition_1 = c1,
            time_effect_2 = time_effect_2,
            condition_2 = c2,
            avrg_diff_conditions = avrg_diff_df,
            interaction_condition_time = interaction_df,
            data = data,
            meta = meta,
            condition = condition,
            replicate_column = replicate_column,
            predicted_timecurves = predicted_timecurves,
            adj_pthresh_avrg_diff_conditions =
                adj_pthresh_avrg_diff_conditions,
            adj_pthresh_interaction = adj_pthresh_interaction,
            min_effect_size = min_effect_size,
            plot_info = plot_info,
            plot_options = plot_options,
            raw_data = raw_data,
            max_hit_number = max_hit_number
        )
        comparison_plots[[suffix]] <- plots_and_feature_names
    }
    comparison_plots
}


#' Plot Splines for Features Based on Top Table Information
#'
#' @noRd
#'
#' @description This function generates plots for each feature listed in the
#' top table using spline
#' interpolation for fitted values. It creates individual plots for each feature
#' and combines them into a single composite plot.
#'
#' @param top_table A dataframe containing the indices and names of features,
#' along with their
#'                  statistical metrics such as intercepts and spline
#'                  coefficients.
#' @param data A matrix or dataframe containing the raw data values for each
#' feature.
#' @param meta A dataframe containing metadata for the data, including time
#' points.
#' @param predicted_timecurves A list containing:
#'   \describe{
#'     \item{`time_grid`}{A numeric vector of dense time points used for
#'       evaluation.}
#'     \item{`predictions`}{A named list of matrices, one per condition level.
#'       Each matrix contains predicted values (rows = features, columns =
#'       timepoints).}
#'   }
#' @param time_unit_label A string shown in the plots as the unit for the time,
#' such as min or hours.
#' @param plot_info List containing the elements y_axis_label (string),
#'                  time_unit (string), treatment_labels (character vector),
#'                  treatment_timepoints (integer vector). All can also be NA.
#'                  This list is used to add this info to the spline plots.
#'                  time_unit is used to label the x-axis, and treatment_labels
#'                  and -timepoints are used to create vertical dashed lines,
#'                  indicating the positions of the treatments (such as
#'                  feeding, temperature shift, etc.).
#' @param plot_options `list`: Named list controlling optional plot
#'   customization. Supported entries include:
#'   \itemize{
#'     \item \code{cluster_heatmap_columns} `logical(1)`
#'       (default = \code{FALSE}):
#'       Whether to cluster columns in the heatmap.
#'
#'     \item \code{meta_replicate_column} `character(1)`
#'       (default = \code{NULL}):
#'       Column in \code{meta} encoding replicate information. When supplied,
#'       spline plot data points are colored by replicate.
#'
#'     \item \code{condition_colours} `list`
#'       (default = \code{NULL}):
#'       Optional named list mapping condition levels to colours. Names must
#'       correspond to values in the condition column of \code{meta}, and
#'       values must be valid R colour specifications accepted by
#'       \code{ggplot2}. When provided, these colours override the default
#'       colours for the corresponding conditions.
#'   }
#' @param adj_pthreshold Double > 0 and < 1 specifying the adj. p-val threshold.
#' @param replicate_column String specifying the column of the meta dataframe
#' that contains the labels of the replicate measurents. When that is not
#' given, this argument is NULL.
#' @param level Unique value of the meta condition column, such as 'treatment'
#' or 'control'.
#' @param raw_data Optional. Data matrix with the raw (unimputed) data, still
#' containing NA values. When provided, it highlights the datapoints in the
#' spline plots that originally where NA and that were imputed.
#' @param report_info A named list containing report information such as analyst
#'                    name, fixed and random effects, etc.
#'
#' @return A list containing the composite plot and the number of rows used in
#' the plot layout.
#'
#' @importFrom splines ns
#' @importFrom ggplot2 ggplot geom_point geom_line theme_minimal labs theme
#'                     scale_x_continuous annotate
#' @importFrom scales hue_pal
#' @importFrom rlang .data
#'
plot_splines <- function(
        top_table,
        data,
        meta,
        predicted_timecurves,
        time_unit_label,
        plot_info,
        plot_options,
        adj_pthreshold,
        replicate_column,
        level,
        raw_data,
        report_info,
        max_hit_number,
        all_levels_clustering,
        condition
        ) {
    # Sort so that HTML reports are easier to read and comparisons are easier.
    top_table <- top_table |> dplyr::arrange(.data$feature_names)
    smooth_timepoints <- predicted_timecurves$time_grid
    pred_mat_level <- predicted_timecurves$predictions[[level]]
    
    # pick the right clustering sub-result for this level
    level_key <- if (!is.null(condition)) {
        paste0(condition, "_", level)
    } else {
        level
    }
    level_result <- NULL
    if (is.list(all_levels_clustering)) {
        if (!is.null(all_levels_clustering[[level_key]])) {
            level_result <- all_levels_clustering[[level_key]]
        } else if (!is.null(all_levels_clustering[[level]])) {
            # fallback if names don't carry the condition_ prefix
            level_result <- all_levels_clustering[[level]]
        }
    }
    
    # Helper to fetch similarity + cluster id for a feature_nr
    .get_sim_for_feature <- function(level_result, feature_nr) {
        out <- list(sim = NA_real_, cl = NA)
        if (!is.list(level_result)) {
            return(out)
        }
        cq <- level_result$cluster_quality
        ch <- level_result$clustered_hits
        if (is.null(cq) || is.null(cq$per_member) || is.null(ch)) {
            return(out)
        }
        if (!("feature" %in% names(ch)) || !("cluster" %in% names(ch))) {
            return(out)
        }
        idx <- which(ch$feature == as.integer(feature_nr))
        if (length(idx) == 0) {
            return(out)
        }
        out$sim <- cq$per_member[idx[1]]
        out$cl <- ch$cluster[idx[1]]
        out
    }
    
    DoF <- which(names(top_table) == "AveExpr") - 1
    time_points <- meta[["Time"]]
    
    titles <- data.frame(
        FeatureID = top_table$feature_nr,
        feature_names = top_table$feature_names
    )
    
    shape_values <- c( # 16 = circle, 17 = triangle
        "Measured" = 16,
        "Imputed" = 17
    )
    
    plot_list <- list()
    n_hits <- min(
        max_hit_number,
        nrow(top_table)
    )
    
    for (hit in seq_len(n_hits)) {
        hit_index <- as.numeric(top_table$feature_nr[hit])
        feature_name <- top_table$feature_names[hit]
        sim_info <- .get_sim_for_feature(level_result, hit_index)
        sim_str <- if (
            is.finite(sim_info$sim)
        ) {
            sprintf(
                " | sr<sup>2</sup><sub>cc</sub>: %.2f (cl %s)",
                sim_info$sim, as.character(sim_info$cl)
            )
        } else {
            ""
        }
        # cumulative travel (effect size) for this feature in this level
        cum_travel_val <- NA_real_
        es_vec <- predicted_timecurves$time_effect_effect_size[[level]]
        if (!is.null(es_vec)) {
            tmp <- unname(es_vec[feature_name])
            if (length(tmp)) cum_travel_val <- tmp[1]
        }
        fitted_values <- as.numeric(
            pred_mat_level[feature_name, ]
        )
        y_values <- data[hit_index, ]
        
        homosc_result <- report_info[["homosc_violation_result"]][["bp_df"]]
        heteroscedasticity <- homosc_result$violation_flag[hit_index]
        high_var_group <- homosc_result$max_var_group[hit_index]
        
        plot_data <- data.frame(
            Time = time_points,
            Y = y_values
        )
        
        # Mark original NA values from raw_data if available
        if (!is.null(raw_data)) {
            # Identify NA positions and mark them as "Imputed" in `plot_data`
            na_indices <- which(is.na(raw_data[hit_index, ]))
            plot_data$IsNA <- "Measured"
            plot_data$IsNA[na_indices] <- "Imputed"
        } else {
            plot_data$IsNA <- "Measured"
        }
        
        # If replicate_column is specified (i.e., a string), use replicate info
        if (!is.null(replicate_column) && is.character(replicate_column)) {
            replicates <- meta[[replicate_column]] # Get the replicate info
            plot_data$Replicate <- replicates # Add replicate info to plot data
            
            # Create color palette for replicates
            replicate_colors <- scales::hue_pal()(length(unique(replicates)))
            names(replicate_colors) <- unique(replicates)
            
            color_values <- c(
                "Spline" = "red",
                replicate_colors
            )
        } else {
            color_values <- c(
                "Data" = "blue",
                "Spline" = "red"
            )
        }
        
        cc <- plot_options[["condition_colours"]]
        
        if (!is.null(cc) &&
            is.list(cc) &&
            level %in% names(cc) &&
            is.character(cc[[level]]) &&
            length(cc[[level]]) == 1L) {
            color_values["Spline"] <- cc[[level]]
        }
        
        # Get adjusted p-value and significance stars
        adj_p_value <- as.numeric(top_table[hit, "adj.P.Val"])
        significance_stars <- ifelse(
            adj_p_value < adj_pthreshold / 500,
            "****",
            ifelse(
                adj_p_value < adj_pthreshold / 50,
                "***",
                ifelse(
                    adj_p_value < adj_pthreshold / 5,
                    "**",
                    ifelse(
                        adj_p_value < adj_pthreshold,
                        "*",
                        ""
                    )
                )
            )
        )
        
        avg_cv <- calc_cv(
            time_values = time_points,
            response_values = y_values
        )
        
        # Use local environment to avoid unwanted updating dynamic legend label.
        p <- local({
            plot_spline <- data.frame(
                Time = smooth_timepoints,
                Fitted = fitted_values
            )
            
            x_min <- min(time_points)
            x_max <- max(time_points)
            x_extension <- (x_max - x_min) * 0.001 # Etxtension on each side
            
            # Define color column outside aes()
            color_column_values <- if (!is.null(replicate_column) &&
                                       is.character(replicate_column)) {
                plot_data$Replicate # Use replicate column if it exists
            } else {
                rep("Data", nrow(plot_data))
            }
            
            plot_data$color_column <- factor(color_column_values)
            
            y_max <- max(c(y_values, fitted_values), na.rm = TRUE)
            y_min <- min(c(y_values, fitted_values), na.rm = TRUE)
            y_extension <- (y_max - y_min) * 0.1
            
            p <- ggplot2::ggplot() +
                ggplot2::geom_point(
                    data = dplyr::filter(plot_data, !is.na(.data$Y)),
                    ggplot2::aes(
                        x = .data$Time,
                        y = .data$Y,
                        color = .data$color_column,
                        shape = factor(.data$IsNA) 
                    ),
                    alpha = 0.5 # 50% transparent data dots
                ) +
                ggplot2::geom_line(
                    data = plot_spline,
                    ggplot2::aes(
                        x = .data$Time,
                        y = .data$Fitted,
                        color = "Spline"
                    )
                ) +
                ggplot2::scale_shape_manual(values = shape_values) +
                ggplot2::theme_minimal() +
                ggplot2::scale_x_continuous(
                    limits = c(x_min - x_extension, x_max + x_extension),
                    labels = scales::label_number_auto()
                ) +
                ggplot2::guides(x = ggplot2::guide_axis(check.overlap = TRUE)) +
                ggplot2::coord_cartesian(ylim = c(y_min, y_max + y_extension)) +
                ggplot2::labs(
                    x = paste0("Time ", time_unit_label),
                    y = plot_info$y_axis_label
                ) +
                ggplot2::guides(
                    color = ggplot2::guide_legend(title = NULL),
                    shape = if (any(plot_data$IsNA == "Imputed")) {
                        ggplot2::guide_legend(title = NULL)
                    } else {
                        "none" 
                    }
                )
            
            y_pos_label <- y_max + y_extension * 0.5
            
            result <- maybe_add_dashed_lines(
                p = p,
                plot_info = plot_info,
                level = level,
                y_pos = y_pos_label,
                horizontal_labels = TRUE
            )
            
            p <- result$p # Updated plot with dashed lines
            treatment_colors <- result$treatment_colors 
            
            color_values <- c(
                color_values,
                treatment_colors
            )
            
            # Add title and annotations
            matched_row <- dplyr::filter(
                titles,
                !!rlang::sym("FeatureID") == hit_index
            )
            
            title <- as.character(matched_row$feature_name)
            
            if (isTRUE(heteroscedasticity)) {
                if (!is.na(high_var_group)) {
                    title_prefix <- paste0(
                        "\u26A0 (",
                        high_var_group,
                        " \u2191) | "
                    )
                } else {
                    title_prefix <- "\u26A0\uFE0F "
                }
            } else {
                title_prefix <- ""
            }
            
            title <- paste0(
                title_prefix,
                title
            )
            
            if (nchar(title) > 100) {
                title_before <- title
                title <- paste0(substr(title, 1, 100), " ...")
                message(paste(
                    "The feature ID", title_before, "is > 100 characters.",
                    "Truncating it to 100 chars:", title
                ))
            }
            
            if (is.na(title)) {
                title <- paste("feature:", hit_index)
            }
            
            p <- p +
                ggplot2::scale_colour_manual(
                    values = color_values,
                ) +
                ggplot2::labs(
                    title = paste(
                        "<b>", title, "</b>",
                        "<br>",
                        "cT:",
                        ifelse(
                            is.na(cum_travel_val),
                            "NA", signif(cum_travel_val, 3)
                        ),
                        "  |  avg CV: ", round(avg_cv, 2), "%",
                        "  |  adj. p-val: ", signif(adj_p_value, digits = 2),
                        " ", significance_stars,
                        sim_str
                    ),
                    x = paste("Time", time_unit_label),
                    y = paste(plot_info$y_axis_label)
                ) +
                ggplot2::theme(
                    plot.title = ggplot2::element_text(size = 6),
                    axis.title.x = ggplot2::element_text(size = 14),
                    axis.title.y = ggplot2::element_text(size = 14),
                    legend.key.size = grid::unit(0.8, "cm"),
                    legend.key.height = grid::unit(0.5, "cm"),
                    legend.title = ggplot2::element_text(size = 8),
                    legend.text = ggplot2::element_text(size = 12),
                    axis.text.x = ggplot2::element_text(size = 12),
                    axis.text.y = ggplot2::element_text(size = 12)
                )
            
            p
        })
        
        plot_list[[hit]] <- p
    }
    
    return(plot_list)
}


#' Plot Consensus Shapes
#'
#' @noRd
#'
#' @description
#' Generates composite plots of single and consensus shapes for each cluster
#' of curve values.
#'
#' @param curve_values A dataframe containing curve values and cluster
#'  assignments.
#' @param plot_info List containing the elements y_axis_label (string),
#'                  time_unit (string), treatment_labels (character vector),
#'                  treatment_timepoints (integer vector). All can also be NA.
#'                  This list is used to add this info to the spline plots.
#'                  time_unit is used to label the x-axis, and treatment_labels
#'                  and -timepoints are used to create vertical dashed lines,
#'                  indicating the positions of the treatments (such as
#'                  feeding, temperature shift, etc.).
#' @param level Unique value within the condition.
#' @param max_hit_number Maximum number of hits which are plotted within each
#' cluster. This can be used to limit the computation time and size of
#' the HTML report in the case of many hits.
#'
#' @return A list containing a plot for every cluster
#'
#' @seealso
#' \code{\link{plot_single_and_mean_splines}}
#'
plot_cluster_mean_splines <- function(
        curve_values,
        plot_info,
        level
        ) {
    clusters <- sort(unique(curve_values$cluster))
    plots <- list()
    cluster_colors <- get_cluster_colors(curve_values)
    
    for (current_cluster in clusters) {
        subset_df <- subset(
            curve_values,
            curve_values$cluster == current_cluster
        )
        
        nr_of_hits <- nrow(subset_df)
        
        subset_df$cluster <- NULL
        current_title <- paste(
            "Cluster",
            current_cluster,
            "| Hits:",
            nr_of_hits,
            "|",
            level,
            sep = " "
        )
        
        plots[[length(plots) + 1]] <-
            plot_single_and_mean_splines(
                subset_df,
                current_title,
                plot_info = plot_info,
                level,
                cluster_color = 
                    cluster_colors[[paste("Cluster", current_cluster)]]
            )
    }
    return(plots)
}


#' Plot Heatmap
#'
#' @noRd
#'
#' @description
#' Generates heatmaps for each level within a condition, showing z-scores of
#' log2 intensity values, split by clusters.
#'
#' @param datas A matrix of data values.
#' @param meta A dataframe containing metadata.
#' @param mode A character vector with length 1, specifying the type of limma
#'             design formula (integrated for formulas with interaction effects
#'             between the levels, isolated for formulas where each level is
#'             analysed in isolation (no interaction effects))
#' @param condition A character string specifying the condition.
#' @param all_levels_clustering A list containing clustering results for each
#' level within the condition.
#' @param time_unit_label A character string specifying the time unit label.
#' @param cluster_heatmap_columns Boolean specifying wether to cluster the
#' columns of the heatmap or not.
#' @param max_hit_number Maximum number of hits which are plotted within each
#' cluster. This can be used to limit the computation time and size of
#' the HTML report in the case of many hits.
#'
#' @return A list of ComplexHeatmap heatmap objects for (one for each level).
#'
#' @seealso
#' \link[ComplexHeatmap]{Heatmap}, \link[dplyr]{arrange}
#'
#' @importFrom dplyr arrange mutate group_by summarize
#' @importFrom tidyr pivot_longer separate
#' @importFrom ComplexHeatmap Heatmap draw ht_opt
#' @importFrom ggplot2 ggplot geom_line facet_wrap geom_vline ylab theme unit
#' @importFrom ggplot2 theme_bw scale_x_continuous
#' @importFrom grid gpar
#' @importFrom rlang .data
#'
plot_heatmap <- function(
        datas,
        meta,
        mode,
        condition,
        all_levels_clustering,
        time_unit_label,
        cluster_heatmap_columns,
        max_hit_number) {
    BASE_TEXT_SIZE_PT <- 5
    
    ht_opt(
        simple_anno_size = unit(1.5, "mm"),
        COLUMN_ANNO_PADDING = unit(1, "pt"),
        DENDROGRAM_PADDING = unit(1, "pt"),
        HEATMAP_LEGEND_PADDING = unit(1, "mm"),
        ROW_ANNO_PADDING = unit(1, "pt"),
        TITLE_PADDING = unit(2, "mm"),
        heatmap_row_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
        heatmap_row_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
        heatmap_column_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
        heatmap_column_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
        legend_labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
        legend_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
        legend_border = FALSE
    )
    
    ht_opt$message <- FALSE
    
    levels <- unique(meta[[condition]])
    heatmaps <- list()
    
    # Generate a heatmap for every level
    for (i in seq_along(all_levels_clustering)) {
        # Skip when no result: NULL, all NA, or informative string
        if (is.null(all_levels_clustering[[i]]) ||
            all(is.na(all_levels_clustering[[i]])) ||
            is.character(all_levels_clustering[[i]])) {
            heatmaps[[length(heatmaps) + 1]] <- NA
            next
        }
        
        level_clustering <- all_levels_clustering[[i]]
        
        clustered_hits <- level_clustering$clustered_hits
        clusters <- clustered_hits |> dplyr::arrange(!!rlang::sym("cluster"))
        
        if (!is.infinite(max_hit_number)) {
            clusters <- clusters |>
                dplyr::group_by(!!rlang::sym("cluster")) |>
                dplyr::slice_head(n = max_hit_number) |>
                dplyr::ungroup()
        }
        
        level <- levels[[i]]
        level_indices <- which(meta[[condition]] == level)
        
        if (mode == "integrated") {
            data_level <- datas[[i]][, level_indices]
        } else { # mode == "isolated"
            data_level <- datas[[i]]
        }
        
        data_level <- data_level[as.numeric(clusters$feature), ]
        z_score <- t(scale(t(data_level)))
        
        meta_level <- meta[level_indices, ]
        
        row_labels <- truncate_row_names(rownames(data_level))
        
        if (is.null(cluster_heatmap_columns)) { # set default value
            cluster_heatmap_columns <- FALSE
        }
        
        ht <-
            ComplexHeatmap::Heatmap(
                z_score,
                name = paste0(
                    "left-labels = cluster,",
                    "top-labels = time"
                ),
                use_raster = TRUE,
                column_split = meta_level$Time,
                cluster_columns = cluster_heatmap_columns,
                row_split = clusters$cluster,
                cluster_rows = FALSE,
                heatmap_legend_param = list(
                    title = "z-score of log2 values",
                    title_position = "lefttop-rot"
                ),
                row_gap = unit(2, "pt"),
                column_gap = unit(2, "pt"),
                show_row_names = TRUE,
                row_labels = row_labels,
                show_column_names = TRUE,
                column_names_rot = 70,
                column_names_gp = gpar(fontsize = 5)
            )
        
        heatmaps[[length(heatmaps) + 1]] <- ht
    }
    # name the list elements by condition level
    names(heatmaps) <- levels[seq_along(heatmaps)]
    heatmaps
}


#' Build Cluster Hits Report
#'
#' @noRd
#'
#' @description
#' Generates an HTML report for clustered hits, including plots and
#' spline parameter details, with a table of contents.
#'
#' @param header_section A character string containing the HTML header section.
#' @param plots A list of ggplot2 plot objects.
#' @param limma_result_2_and_3_plots List containing the list of lists with all
#' the plots for all the pairwise comparisons of the condition in terms of
#' average spline diff and interaction condition time, and another list of lists
#' where the respective names of each plot are stored.
#' @param plots_sizes A list of integers specifying the size of each plot.
#' @param level_headers_info A list of header information for each level.
#' @param spline_params A list of spline parameters.
#' @param adj_pthresh_time_effect `numeric(1)`: adj. p-value threshold
#' for the limma time effect results (category 1).
#' @param adj_pthresh_avrg_diff_conditions Float
#' @param adj_pthresh_interaction_condition_time Float
#' @param row_counts_dict A nested list containing row counts for dataframes
#'   from `category_2_and_3_hits`. The outer keys are `"category_2"` and
#'   `"category_3"`, representing the two sublists. The inner keys are derived
#'   from the portion of each dataframe name after the second underscore (`_`),
#'   or the full name if fewer than two underscores exist. The values are
#'   integers representing the number of rows in each dataframe.
#' @param mode A character string specifying the mode
#'            ('isolated' or 'integrated').
#' @param report_info A named list containg the report info fields. Here used
#'                    for the email hotkey functionality.
#' @param output_file_path A character string specifying the path to save the
#'                         HTML report.
#'
#' @return No return value, called for side effects.
#'
#' @seealso
#' \code{\link{plot2base64}}, \code{\link{create_progress_bar}}
#'
build_cluster_hits_report <- function(
        header_section,
        plots,
        limma_result_2_and_3_plots,
        plots_sizes,
        level_headers_info,
        spline_params,
        adj_pthresh_time_effect,
        adj_pthresh_avrg_diff_conditions,
        adj_pthresh_interaction_condition_time,
        category_2_and_3_hit_counts,
        mode,
        report_info,
        output_file_path) {
    html_content <- paste(header_section, "<!--TOC-->", sep = "\n")
    
    toc <- create_toc()
    
    styles <- define_html_styles()
    section_header_style <- styles$section_header_style
    toc_style <- styles$toc_style
    
    current_header_index <- 1
    j <- 0
    level_headers_info <- Filter(
        Negate(is.null),
        level_headers_info
    )
    
    pb <- create_progress_bar(plots)
    
    header_index <- 0
    level_index <- 0
    
    # Generate the sections and plots
    for (index in seq_along(plots)) {
        header_index <- header_index + 1
        
        if (current_header_index <= length(level_headers_info)) {
            header_info <- level_headers_info[[current_header_index]]
            nr_hits <- header_info$nr_hits
            adj_pvalue_threshold <- header_info$adj_pvalue_threshold
            
            # means this is the section of a new level
            # The very first level is also a new level
            if (names(plots)[index] == "new_level") {
                level_index <- level_index + 1
                
                time_effect_section_header <- paste(
                    "Time Effect in Condition:",
                    header_info$header_name
                )
                
                section_header <- sprintf(
                    "<h2 style='%s' id='section%d'>%s</h2>",
                    section_header_style,
                    header_index,
                    time_effect_section_header
                )
                
                html_content <- paste(
                    html_content,
                    section_header,
                    sep = "\n"
                )
                
                if (mode == "integrated") {
                    j <- 1
                } else { # mode == "isolated" or mode == NA
                    j <- j + 1
                }
                
                spline_params_info <-
                    get_spline_params_info(
                        spline_params = spline_params,
                        j = j
                    )
                
                html_content <- paste(
                    html_content,
                    spline_params_info,
                    sep = "\n"
                )
                
                hits_info <- sprintf(
                    paste0(
                        "<p style='text-align: center; font-size: 30px;'>",
                        "adj.p-value threshold: %.4g</p>",
                        "<p style='text-align: center; font-size: 30px;'>",
                        "Number of hits: %d</p>",
                        "<div style='text-align: center; font-size:",
                        "30px;'>%s</div>",
                        "<hr>"
                    ),
                    adj_pvalue_threshold,
                    nr_hits,
                    generate_asterisks_definition(
                        adj_pvalue_threshold = adj_pvalue_threshold,
                        include_ns = FALSE
                    )
                )
                
                html_content <- paste(
                    html_content,
                    hits_info,
                    sep = "\n"
                )
                
                toc_entry <- sprintf(
                    "<li style='%s'><a href='#section%d'>%s</a></li>",
                    toc_style,
                    header_index,
                    time_effect_section_header
                )
                toc <- paste(
                    toc,
                    toc_entry,
                    sep = "\n"
                )
                
                current_header_index <- current_header_index + 1
                
                pb$tick()
                next
            }
        }
        
        element_name <- names(plots)[index]
        
        header_levels <- c(
            "dendrogram",
            "cluster_mean_splines",
            "cluster_quality_plots",
            "heatmap",
            "individual_spline_plots"
        )
        
        if (element_name %in% header_levels) {
            if (element_name == "dendrogram") {
                header_text <- "Overall Clustering"
            } else if (element_name == "cluster_mean_splines") {
                header_text <- "Z-score normalized individual and mean splines"
            } else if (element_name == "cluster_quality_plots") {
                header_text <-
                    "Variance-explained-by-cluster-centroid distribution plots"
                cquality_description <- paste(
                    "<div style='text-align: center; font-size: 1.5em;'>",
                    "Variance explained by cluster centroid:",
                    "<br>",
                    "<ul style='list-style-position: inside; text-align: left;",
                    "display: inline-block;'>",
                    "<li>0.90-1.00 = excellent</li>",
                    "<li>0.80-0.89 = very good</li>",
                    "<li>0.70-0.79 = good</li>",
                    "<li>0.60-0.69 = borderline</li>",
                    "<li>0.50-0.59 = poor</li>",
                    "<li>0.00-0.49 = very poor</li>",
                    "<li>&lt;0.00 = anti-pattern</li>",
                    "</ul>",
                    "<br>",
                    "sr<sup>2</sup><sub>cc</sub> = signed r<sup>2</sup>",
                    "by cluster centroid,",
                    "i.e. how well a gene's spline fits the centroid of",
                    "its assigned cluster.",
                    "<br><hr>",
                    "</div>"
                )
            } else if (element_name == "heatmap") {
                header_text <- "Z-Score of log2 Value Heatmap"
                
                heatmap_description <- paste(
                    "<div style='text-align: center; font-size: 1.5em;'>",
                    "Rows = features (labels on the right, cluster labels",
                    "on the left),",
                    "columns = timepoints; Blue = down, red = up, --> compared",
                    "to the rest of the row;",
                    "</div>"
                )
            } else { # element_name == "individual_spline_plots"
                adjusted_p_val <- adj_pthresh_time_effect
                header_text <- "Individual Significant Features (Hits) Splines"
                asterisks_definition <- 
                    generate_asterisks_definition(
                        adj_pvalue_threshold = adj_pvalue_threshold,
                        include_ns = FALSE
                    )
            }
            
            # Add the main title as a section title with an anchor
            # before the first plot
            header <- paste0(
                "<h2 id='section",
                header_index,
                "' style='text-align: center; font-size: 3.5em;'>",
                header_text,
                "</h2>",
                if (exists("heatmap_description")) {
                    heatmap_description
                } else if (exists("cquality_description")) {
                    cquality_description
                } else {
                    ""
                }
            )
            
            if (exists("cquality_description")) rm(cquality_description)
            if (exists("heatmap_description")) rm(heatmap_description)
            
            # Add the asterisks definition if it exists
            if (exists("asterisks_definition")) {
                header <- paste0(
                    header,
                    "<div style='text-align: center;",
                    "font-size: 1.5em;'>",
                    asterisks_definition,
                    "</div>"
                )
                # Otherwise, the next level has it everywhere
                rm(asterisks_definition) 
            }
            
            html_content <- paste(
                html_content,
                header,
                sep = "\n"
            )
            
            toc_entry <- paste0(
                "<li style='margin-left: 30px; font-size: 30px;'>",
                "<a href='#section",
                header_index,
                "'>",
                header_text,
                "</a></li>"
            )
            
            toc <- paste(toc, toc_entry, sep = "\n")
        }
        
        header_index <- header_index + 1
        
        result <- process_plots(
            plots_element = plots[[index]],
            element_name = names(plots)[index],
            plots_size = plots_sizes[[index]],
            html_content = html_content,
            toc = toc,
            header_index = header_index
        )
        
        html_content <- result$html_content
        toc <- result$toc
        
        pb$tick()
    }
    pb$terminate()
    
    # Avrg diff conditions & interaction condition time (per contrast)
    if (length(limma_result_2_and_3_plots) > 0) {
        # Main header for the whole limma section
        header_index <- header_index + 1
        
        limma_main_header <- sprintf(
            "<h2 style='%s' id='section%d'>%s</h2>",
            section_header_style,
            header_index,
            "Avrg diff conditions & interaction condition time"
        )
        
        html_content <- paste(
            html_content,
            limma_main_header,
            sep = "\n"
        )
        
        asterisks_definition_avrg_diff <- paste(
            "<div style='text-align:center; margin-bottom: 20px;'>",
            generate_asterisks_definition(
                adj_pvalue_threshold =
                    adj_pthresh_avrg_diff_conditions,
                header_html = paste(
                    "<b><span style='font-size:24pt;'>",
                    "Asterisks definition (Average Diff Conditions):",
                    "</span></b>",
                    sep = ""
                )
            ),
            "</div>",
            sep = "\n"
        )
        
        asterisks_definition_interaction <- paste(
            "<div style='text-align:center; margin-bottom: 40px;'>",
            generate_asterisks_definition(
                adj_pvalue_threshold =
                    adj_pthresh_interaction_condition_time,
                header_html = paste(
                    "<b><span style='font-size:24pt;'>",
                    "Asterisks definition (Interaction):",
                    "</span></b>",
                    sep = ""
                )
            ),
            "</div>",
            sep = "\n"
        )
        
        html_content <- paste(
            html_content,
            asterisks_definition_avrg_diff,
            asterisks_definition_interaction,
            sep = "\n"
        )
        
        # TOC entry for the whole limma section
        toc_entry <- sprintf(
            "<li style='%s'><a href='#section%d'>%s</a></li>",
            toc_style,
            header_index,
            "Avrg diff conditions & interaction condition time"
        )
        toc <- paste(
            toc,
            toc_entry,
            sep = "\n"
        )
        
        # Per-contrast subsections
        c2_counts <- category_2_and_3_hit_counts[["category_2"]]
        c3_counts <- category_2_and_3_hit_counts[["category_3"]]
        
        for (comparison_name in names(limma_result_2_and_3_plots)) {
            # Build title including comparison name
            header_index <- header_index + 1
            subsection_title <- comparison_name
            
            subheader <- sprintf(
                "<h3 style='font-size: 3.5em; color: #001F3F;
                text-align: center;' id='section%d'>%s</h3>",
                header_index,
                subsection_title
            )
            
            # Map suffix to count names
            avrg_key <- paste0("avrg_diff_", comparison_name)
            inter_key <- paste0("time_interaction_", comparison_name)
            
            avrg_diff_hits <- if (!is.null(c2_counts) &&
                                  avrg_key %in% names(c2_counts)) {
                c2_counts[[avrg_key]]
            } else 0L
            
            interaction_hits <- if (!is.null(c3_counts) &&
                                    inter_key %in% names(c3_counts)) {
                c3_counts[[inter_key]]
            } else 0L
            
            hits_info <- sprintf(
                paste0(
                    "<p style='font-size: 2em; text-align: center;'>",
                    "Avrg diff conditions hits: %d</p>",
                    "<p style='font-size: 2em; text-align: center;'>",
                    "Interaction condition time hits: %d</p>",
                    "<hr>"
                ),
                avrg_diff_hits,
                interaction_hits
            )
            
            html_content <- paste(
                html_content,
                subheader,
                hits_info,
                sep = "\n"
            )
            
            # TOC entry for this comparison
            toc_entry <- paste0(
                "<li style='margin-left: 30px; font-size: 30px;'>",
                "<a href='#section", header_index, "'>",
                subsection_title,
                "</a></li>"
            )
            toc <- paste(
                toc,
                toc_entry,
                sep = "\n"
            )
            
            # Extract plots + feature names for this comparison
            comparison <- limma_result_2_and_3_plots[[comparison_name]]
            comparison_plots <- comparison$plots
            comparison_feature_names <- comparison$feature_names
            
            for (i in seq_along(comparison_plots)) {
                
                feature_name_div <- sprintf(
                    '<div style="text-align: center;
                    font-size: 36px; margin-bottom: 10px;">%s</div>',
                    comparison_feature_names[[i]]
                )
                
                # Corrected: build the new HTML first, then pass it in
                updated_html <- paste(
                    html_content,
                    feature_name_div,
                    sep = "\n"
                )
                
                result <- process_plots(
                    plots_element = comparison_plots[[i]],
                    plots_size    = 1.5,
                    html_content  = updated_html,
                    toc           = toc,
                    header_index  = header_index,
                    element_name  = ""
                )
                
                html_content <- result$html_content
                toc <- result$toc
            }
        }
    }
    
    generate_and_write_html(
        toc = toc,
        html_content = html_content,
        report_info = report_info,
        output_file_path = output_file_path
    )
}


# Level 3 internal functions ---------------------------------------------------


#' Plot the signed r^2 distribution for one cluster
#'
#' @noRd
#'
#' @description
#' Draws a histogram of per-member \emph{signed} r\eqn{^2} values
#' (variance explained with sign of the correlation) to the cluster centroid
#' for a single cluster. A red dashed vertical line marks the mean signed
#' r\eqn{^2} of that cluster, and a dotted vertical line at 0 highlights the
#' boundary between in-phase (\eqn{>0}) and inverted (\eqn{<0}) shapes.
#'
#' @param r2 Numeric vector of \emph{signed} r\eqn{^2} for members
#'   of one cluster (computed from Pearson r to that cluster's centroid as
#'   \code{sign(r) * r^2}). Values are in \code{[-1, 1]}; \code{NA}s allowed.
#' @param cluster_id Integer or character identifier shown in the plot title.
#'
#' @return A ggplot object showing the signed r\eqn{^2} histogram with mean
#'   and reference lines.
#'
#' @details
#' The histogram uses fixed bin width \code{0.1} over \code{[-1, 1]} (20 bins)
#' and displays counts on the y-axis. The mean line is drawn at
#' \code{mean(r2, na.rm = TRUE)}.
#'
#' @importFrom ggplot2 ggplot aes geom_histogram after_stat geom_vline theme
#'                     scale_color_manual labs coord_cartesian theme_minimal
#' @importFrom rlang .data
#'
plot_cluster_quality_distribution <- function(
        r2,
        cluster_id) {
    df <- data.frame(sr2 = r2)
    mean_sr2 <- mean(df$sr2, na.rm = TRUE)
    n_valid <- sum(is.finite(df$sr2))
    
    p <- ggplot2::ggplot(df, ggplot2::aes(x = !!rlang::sym("sr2")))
    
    # histogram with fixed 0.1 bins from -1 to 1
    if (n_valid >= 1) {
        p <- p + ggplot2::geom_histogram(
            binwidth = 0.05,
            boundary = -1,
            closed   = "right",
            fill     = "steelblue",
            color    = NA,
            na.rm    = TRUE
        )
    }
    
    # reference line at 0 (in-phase vs inverted boundary)
    p <- p + ggplot2::geom_vline(
        xintercept = 0,
        linetype = "dotted",
        color = "grey50"
    )
    
    # mean line + legend only if the mean is finite
    if (is.finite(mean_sr2)) {
        p <- p + ggplot2::geom_vline(
            ggplot2::aes(xintercept = mean_sr2, color = "mean"),
            linetype = "dashed"
        ) +
            ggplot2::scale_color_manual(
                values = c("mean" = "red"),
                name = NULL
            )
    }
    
    p +
        ggplot2::labs(
            title = paste0(
                "Cluster ", cluster_id,
                " variances explained by centroid (mean = ",
                round(mean_sr2, 3), ")"
            ),
            x = bquote("Signed " ~ r^2 ~ "(variance explained)"),
            y = "Hit Count"
        ) +
        ggplot2::coord_cartesian(xlim = c(-1, 1)) +
        ggplot2::theme_minimal(base_size = 13) +
        ggplot2::theme(
            legend.position = "right",
            legend.key.size = grid::unit(0.9, "cm"),
            legend.key.height = grid::unit(0.6, "cm"),
            aspect.ratio = 0.5,
            axis.title.x = ggplot2::element_text(size = 14),
            axis.title.y = ggplot2::element_text(size = 14),
            axis.text.x = ggplot2::element_text(size = 12),
            axis.text.y = ggplot2::element_text(size = 12),
            legend.text = ggplot2::element_text(size = 12)
        )
}


#' Merge Annotation with a Single Top Table
#'
#' @noRd
#'
#' @description
#' This function merges annotation information into a single `top_table`
#' dataframe based on the `feature_nr` column.
#'
#' @param top_table A dataframe containing the `top_table` with a `feature_nr`
#'                  column.
#' @param annotation A dataframe containing the annotation information.
#'
#' @return A dataframe with updated `top_table` containing merged annotation
#'         information.
#'
merge_top_table_with_annotation <- function(
        top_table,
        annotation) {
    if (!"feature_nr" %in% names(top_table)) {
        top_table$feature_nr <- seq_len(nrow(top_table))
    } else {
        top_table$feature_nr <- as.integer(as.character(top_table$feature_nr))
    }
    annotation_rows <- annotation[top_table$feature_nr, ]
    top_table <- cbind(top_table, annotation_rows)
}


#' Create spline comparison plots for two conditions
#'
#' @noRd
#'
#' @description
#' This function generates comparison plots for spline fits of two conditions
#' over time. It compares the time effects of two conditions, plots the data
#' points, and overlays the fitted spline curves. The function checks if the
#' adjusted p-values for the average difference between conditions and the
#' interaction between condition and time are below the specified thresholds
#' before generating plots. (this function generates the double spline plots).
#'
#' @param time_effect_1 A data frame containing the time effects for the first
#'  condition.
#' @param condition_1 The name of the first condition.
#' @param time_effect_2 A data frame containing the time effects for the second
#'  condition.
#' @param condition_2 The name of the second condition.
#' @param avrg_diff_conditions A data frame with the adjusted p-values for the
#'  average difference
#' between conditions.
#' @param interaction_condition_time A data frame with the adjusted p-values
#'  for the interaction between
#' condition and time.
#' @param data The data matrix containing the measurements.
#' @param meta The metadata associated with the measurements.
#' @param condition Column name of meta that contains the levels of the
#' experiment.
#' @param replicate_column Column name of the meta column that specifies the
#' replicates per timepoint. For example Reactor with the unique values:
#' 'ReactorE16', 'ReactorE17', ... which means that multiple bioreactors where
#' running this experiment and each timepoint has one sample from each reactor.
#' @param predicted_timecurves A list containing:
#'   \describe{
#'     \item{`time_grid`}{A numeric vector of dense time points used for
#'       evaluation.}
#'     \item{`predictions`}{A named list of matrices, one per condition level.
#'       Each matrix contains predicted values (rows = features, columns =
#'       timepoints).}
#'   }
#' @param plot_info A list containing plotting information such as time unit
#' and axis labels.
#' @param plot_options `list`: Named list controlling optional plot
#'   customization. Supported entries include:
#'   \itemize{
#'     \item \code{cluster_heatmap_columns} `logical(1)`
#'       (default = \code{FALSE}):
#'       Whether to cluster columns in the heatmap.
#'
#'     \item \code{meta_replicate_column} `character(1)`
#'       (default = \code{NULL}):
#'       Column in \code{meta} encoding replicate information. When supplied,
#'       spline plot data points are colored by replicate.
#'
#'     \item \code{condition_colours} `list`
#'       (default = \code{NULL}):
#'       Optional named list mapping condition levels to colours. Names must
#'       correspond to values in the condition column of \code{meta}, and
#'       values must be valid R colour specifications accepted by
#'       \code{ggplot2}. When provided, these colours override the default
#'       colours for the corresponding conditions.
#'   }
#' @param adj_pthresh_avrg_diff_conditions The adjusted p-value threshold for
#' the average difference
#' between conditions.
#' @param adj_pthresh_interaction The adjusted p-value threshold for the
#' interaction between
#' condition and time.
#' @param min_effect_size `list`: A named list that specifies the minimum 
#' effect size thresholds to consider a feature as biologically meaningful, in 
#' addition to statistical significance. This allows users to filter out 
#' "trivial" hits that pass adjusted p-value cutoffs but show negligible effect 
#' sizes.
#'
#'   The list must contain the following elements:
#'   - `time_effect`: `numeric(1)` Minimum cumulative travel for time effects 
#'     (Category 1). Features with a smaller travel will be ignored even if 
#'     significant.
#'   - `avg_diff_cond`: `numeric(1)` Minimum absolute effect size for average 
#'     differences between conditions (Category 2). Ensures that only contrasts 
#'     with a relevant magnitude are reported.
#'   - `interaction_cond_time`: `numeric(1)` Minimum effect size for the 
#'     interaction between condition and time (Category 3). This controls how 
#'     large the differential curve travel must be across conditions to count 
#'     as a hit.
#'
#'   Values should be numeric scalars (typically >0). For example: 
#'   `min_effect_size = list(time_effect = 1, avg_diff_cond = 1, 
#'   interaction_cond_time = 2)` will only keep features with cumulative 
#'   travels or condition-time differences above those cutoffs. Use smaller 
#'   values (e.g., 0.1) for permissive filtering, or larger values for more 
#'   conservative thresholds.
#'
#'   The default is 0 for all three elements.
#' @param raw_data Optional. Data matrix with the raw (unimputed) data, still
#' containing NA values. When provided, it highlights the datapoints in the
#' spline plots that originally where NA and that were imputed.
#' @param max_hit_number Maximum number of hits for which the individual spline
#' plots are shown. This can be used to limit the computation time and size of
#' the HTML report in the case of many hits.
#'
#' @return A list containing:
#' \describe{
#'   \item{plots}{A list of ggplot2 plots comparing the two conditions.}
#'   \item{feature_names}{A list of feature names for the plotted features.}
#' }
#'
#' @importFrom rlang .data
#'
plot_spline_comparisons <- function(
        time_effect_1,
        condition_1,
        time_effect_2,
        condition_2,
        avrg_diff_conditions,
        interaction_condition_time,
        data,
        meta,
        condition,
        replicate_column,
        predicted_timecurves,
        plot_info,
        plot_options,
        adj_pthresh_avrg_diff_conditions,
        adj_pthresh_interaction,
        min_effect_size,
        raw_data,
        max_hit_number
) {
    replicate_mapping <- .make_replicate_mapping(
        meta,
        replicate_column
    )
    
    shape_mapping <- .make_shape_mapping(
        meta,
        replicate_column
    )
    
    sorted <- .sort_inputs_for_plotting(
        time_effect_1,
        time_effect_2,
        avrg_diff_conditions,
        interaction_condition_time
    )
    time_effect_1 <- sorted$time_effect_1
    time_effect_2 <- sorted$time_effect_2
    avrg_diff_conditions <- sorted$avrg_diff_conditions
    interaction_condition_time <- sorted$interaction_condition_time
    
    mats <- .get_prediction_mats(
        predicted_timecurves,
        condition_1,
        condition_2
    )
    smooth_timepoints <- mats$smooth_timepoints
    pred_mat_1 <- mats$pred_mat_1
    pred_mat_2 <- mats$pred_mat_2
    
    features_to_plot <- .select_features_to_plot(
        avrg_diff_conditions,
        interaction_condition_time,
        max_hit_number,
        adj_pthresh_avrg_diff_conditions,
        adj_pthresh_interaction
    )
    
    features_to_plot <- .filter_features_in_prediction_mats(
        features_to_plot,
        pred_mat_1,
        pred_mat_2
    )
    
    plot_list <- list()
    feature_names_list <- list()
    
    for (i in seq_len(nrow(features_to_plot))) {
        hit_index <- as.numeric(features_to_plot$feature_nr[i])
        feature_name <- features_to_plot$feature_names[i]
        
        cat2_eff <- .get_cat2_effect(
            avrg_diff_conditions,
            feature_name
        )
        
        cat3 <- .get_cat3_effects(
            predicted_timecurves,
            condition_1,
            condition_2,
            feature_name
        )
        es1 <- cat3$es1
        es2 <- cat3$es2
        diff_es <- cat3$diff_es
        
        avrg_diff_pval <- safe_pull_pval(
            avrg_diff_conditions,
            feature_name,
            "adj.P.Val"
        )
        
        interaction_pval <- safe_pull_pval(
            interaction_condition_time,
            feature_name,
            "adj.P.Val"
        )
        
        cat2_ok <- .is_sig_and_rel(
            avrg_diff_pval,
            adj_pthresh_avrg_diff_conditions,
            cat2_eff,
            min_effect_size[["avg_diff_cond"]]
        )
        
        cat3_ok <- .is_sig_and_rel(
            interaction_pval,
            adj_pthresh_interaction,
            diff_es,
            min_effect_size[["interaction_cond_time"]]
        )
        
        if (!(cat2_ok || cat3_ok)) {
            next
        }
        
        pd <- .build_plot_data(
            data,
            raw_data,
            meta,
            condition,
            replicate_column,
            replicate_mapping,
            hit_index,
            condition_1,
            condition_2
        )
        plot_data <- pd$plot_data
        has_imputed_1 <- pd$has_imputed_1
        has_imputed_2 <- pd$has_imputed_2
        
        fitted_values_1 <- as.numeric(
            pred_mat_1[feature_name, ]
        )
        fitted_values_2 <- as.numeric(
            pred_mat_2[feature_name, ]
        )
        
        cv_1 <- calc_cv(
            time_values = plot_data$Time,
            response_values = plot_data$Y1
        )
        cv_2 <- calc_cv(
            time_values = plot_data$Time,
            response_values = plot_data$Y2
        )
        
        title_lines <- .build_title_lines(
            feature_name = feature_name,
            avrg_diff_pval = avrg_diff_pval,
            interaction_pval = interaction_pval,
            adj_pthresh_avrg_diff_conditions =
                adj_pthresh_avrg_diff_conditions,
            adj_pthresh_interaction =
                adj_pthresh_interaction,
            cat2_eff = cat2_eff,
            es1 = es1,
            es2 = es2,
            diff_es = diff_es,
            condition_1 = condition_1,
            condition_2 = condition_2,
            min_effect_size = min_effect_size,
            cv_1 = cv_1,
            cv_2 = cv_2
        )
        
        p <- .build_comparison_plot(
            plot_data = plot_data,
            smooth_timepoints = smooth_timepoints,
            fitted_values_1 = fitted_values_1,
            fitted_values_2 = fitted_values_2,
            condition_1 = condition_1,
            condition_2 = condition_2,
            replicate_column = replicate_column,
            shape_mapping = shape_mapping,
            plot_info = plot_info,
            plot_options = plot_options,
            title_lines = title_lines,
            has_imputed_1 = has_imputed_1,
            has_imputed_2 = has_imputed_2
        )
        
        plot_list[[length(plot_list) + 1]] <- p
        feature_names_list[[length(feature_names_list) + 1]] <-
            feature_name
    }
    
    list(
        plots = plot_list,
        feature_names = feature_names_list
    )
}


#' Generate consistent colors for clusters
#'
#' @noRd
#'
#' @description
#' This internal helper function assigns a distinct color to each
#' cluster in a dataset, ensuring that plots use a reproducible and
#' consistent color mapping across different functions.
#'
#' @param curve_values A data frame containing a column named
#'   `cluster` with cluster assignments (numeric or coercible to numeric).
#'
#' @return
#' A named character vector of hex color codes. Names are of the
#' form `"Cluster 1"`, `"Cluster 2"`, etc., corresponding to the
#' sorted cluster IDs.
#'
#' @details
#' Colors are generated with [scales::hue_pal()], which cycles
#' through distinct hues on the color wheel. The number of colors
#' is determined by the number of unique clusters present.
#'
#' @importFrom scales hue_pal
#'
get_cluster_colors <- function(curve_values) {
    lvls <- sort(unique(as.numeric(curve_values$cluster)))
    base <- scales::hue_pal()(length(lvls))
    names(base) <- paste("Cluster", lvls)
    base
}


#' Plot Single and Mean Splines
#'
#' @noRd
#'
#' @description
#' Generates a plot showing individual time series shapes and their consensus
#' (mean) shape.
#'
#' @param time_series_data A dataframe or matrix with time series data.
#' @param title A character string specifying the title of the plot.
#' @param plot_info List containing the elements y_axis_label (string),
#'                  time_unit (string), treatment_labels (character vector),
#'                  treatment_timepoints (integer vector). All can also be NA.
#'                  This list is used to add this info to the spline plots.
#'                  time_unit is used to label the x-axis, and treatment_labels
#'                  and -timepoints are used to create vertical dashed lines,
#'                  indicating the positions of the treatments (such as
#'                  feeding, temperature shift, etc.).
#' @param level Level of the condition, which is a factor (categorical
#'              predictor of the linear model)
#' @param cluster_color Color to be used for the splines of this cluster.
#'
#' @return A ggplot object representing the single and consensus shapes.
#'
#' @seealso
#' \code{\link{ggplot2}}
#'
#' @importFrom dplyr arrange mutate
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot geom_line scale_colour_manual theme_minimal
#'                     ggtitle aes labs element_rect
#' @importFrom rlang sym .data
#' @importFrom scales hue_pal
#'
plot_single_and_mean_splines <- function(
        time_series_data,
        title,
        plot_info,
        level,
        cluster_color
        ) {
    time_col <- rlang::sym("time")
    feature_col <- rlang::sym("feature")
    
    # Convert data to long format
    df_long <- as.data.frame(t(time_series_data)) |>
        tibble::rownames_to_column(var = "time") |>
        tidyr::pivot_longer(
            cols = -!!time_col,
            names_to = "feature",
            values_to = "intensity"
        ) |>
        dplyr::arrange(!!feature_col) |>
        dplyr::mutate(time = as.numeric(.data$time))
    
    max_lines <- 100L
    
    features <- unique(df_long$feature)
    
    if (length(features) > max_lines) {
        keep_features <- sample(features, size = max_lines)
        df_long <- df_long[df_long$feature %in% keep_features, , drop = FALSE]
    }
    
    # Compute consensus (mean of each column)
    consensus <- colMeans(time_series_data, na.rm = TRUE)
    
    consensus_df <- data.frame(
        time = as.numeric(colnames(time_series_data)),
        consensus = consensus
    )
    
    time_unit_label <- paste0("[", plot_info$time_unit, "]")
    
    color_values <- c(
        "Mean"   = "black",
        "Spline" = cluster_color
    )
    
    # Draw only the individual splines here
    p <- ggplot2::ggplot() +
        ggplot2::geom_line(
            data = df_long,
            ggplot2::aes(
                x = !!rlang::sym("time"),
                y = !!rlang::sym("intensity"),
                group = !!rlang::sym("feature"),
                colour = "Spline"
            ),
            alpha = 0.4, linewidth = 0.5
        )
    
    treatment_labels <- NA
    
    result <- maybe_add_dashed_lines(
        p = p,
        plot_info = plot_info,
        level = level
    )
    
    p <- result$p
    treatment_colors <- result$treatment_colors
    
    # Combine the original colors with the treatment colors
    color_values <- c(color_values, treatment_colors)
    
    # Add the final scale for colors and adjust legend
    p <- result$p
    treatment_colors <- result$treatment_colors
    
    color_values <- c(color_values, treatment_colors)
    
    p <- p +
        ggplot2::scale_colour_manual(
            name = "",
            values = color_values
        ) +
        ggplot2::coord_cartesian(clip = "off") +
        ggplot2::theme_minimal() +
        ggplot2::labs(
            title = title,
            x = paste("Time", time_unit_label),
            y = paste("z-score norm.", plot_info$y_axis_label)
        ) +
        ggplot2::theme(
            plot.margin = grid::unit(c(1, 1, 1.5, 1), "lines"),
            legend.position = "right",
            legend.box = "vertical",
            legend.background = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(hjust = 0.5),
            legend.key.size = grid::unit(0.9, "cm"),
            legend.key.height = grid::unit(0.6, "cm"),
            axis.title.x = ggplot2::element_text(size = 14),
            axis.title.y = ggplot2::element_text(size = 8),
            axis.text.x = ggplot2::element_text(size = 12),
            axis.text.y = ggplot2::element_text(size = 12),
            legend.text = ggplot2::element_text(size = 12)
        )
    
    # Add mean line last so that it is always on top.
    p <- p +
        ggplot2::geom_line(
            data = consensus_df,
            ggplot2::aes(
                x = !!rlang::sym("time"),
                y = consensus,
                colour = "Mean"
            ),
            linewidth = 1.5
        )
    
    return(p)
}


#' Get Spline Parameters Info
#'
#' @noRd
#'
#' @description
#' This function retrieves the spline parameters information for a given index.
#' It ensures the spline parameters are valid and constructs an HTML string
#' describing the spline parameters.
#'
#' @param spline_params A list containing the spline parameters. The list should
#'                      include elements: `spline_type`, `degree`, and `dof`.
#' @param j An integer specifying the index of the spline parameters to
#' retrieve.
#'
#' @details
#' The function checks if the spline parameters are not `NULL` and have a length
#' greater than or equal to the specified index `j`. If a parameter is
#' invalid or
#' missing, it sets the parameter to `NA`. It then constructs an HTML string
#' describing the spline parameters, including spline type, degree, degrees of
#' freedom (DoF).
#'
#' @return A character string containing HTML-formatted information about the
#'         spline parameters at the specified index.
#'
get_spline_params_info <- function(
        spline_params,
        j) {
    if (!is.null(spline_params$spline_type) &&
        length(spline_params$spline_type) >= j) {
        spline_params$spline_type[j] <- spline_params$spline_type[j]
    } else {
        spline_params$spline_type[j] <- NA
    }
    
    if (!is.null(spline_params$degree) &&
        length(spline_params$degree) >= j) {
        spline_params$degree[j] <- spline_params$degree[j]
    } else {
        spline_params$degree[j] <- NA
    }
    
    if (!is.null(spline_params$dof) &&
        length(spline_params$dof) >= j) {
        spline_params$dof[j] <- spline_params$dof[j]
    } else {
        spline_params$dof[j] <- NA
    }
    
    if (spline_params$spline_type[j] == "b") {
        spline_params_info <- sprintf(
            "
    <p style='text-align: center; font-size: 30px;'>
        <span style='color: blue;'>Spline-type:</span> B-spline<br>
        <span style='color: blue;'>Degree:</span> %s<br>
        <span style='color: blue;'>DoF:</span> %s<br>
    </p>",
            spline_params$degree[j], spline_params$dof[j]
        )
    } else { # spline_type == "n"
        spline_params_info <- sprintf(
            "
    <p style='text-align: center; font-size: 30px;'>
        <span style='color: blue;'>Spline-type:</span> Natural cubic spline<br>
        <span style='color: blue;'>DoF:</span> %s<br>
    </p>",
            spline_params$dof[j]
        )
    }
    return(spline_params_info)
}


#' Truncate Row Names
#'
#' @noRd
#'
#' @description
#' This function truncates row names that exceed a specified maximum length.
#' If the row name length exceeds the maximum length, it appends " ..."
#' to indicate truncation.
#'
#' @param names A character vector of row names.
#' @param max_length An integer specifying the maximum length of the row names.
#' Default is 40.
#'
#' @return A character vector of truncated row names.
#'
truncate_row_names <- function(
        names,
        max_length = 40) {
    vapply(names, function(x) {
        if (nchar(x) > max_length) {
            return(paste0(substr(x, 1, max_length - 3), " ..."))
        } else {
            return(x)
        }
    }, character(1))
}


#' Calculate average CV across unique time points
#'
#' @noRd
#'
#' @description
#' This function calculates the coefficient of variation (CV) for each unique
#' time point based on the provided time values and response values. It then
#' returns the average CV across all time points. The CV is only calculated if
#' there are more than one valid (non-NA) values for a given time point and
#' the mean of the values is non-zero.
#'
#' @param time_values A numeric vector containing the time points. Time points
#' may repeat across replicates.
#' @param response_values A numeric vector of response values corresponding to
#' the time points.
#'
#' @return The average coefficient of variation (CV) across all time points.
#' Returns NA if all CVs are NA.
#'
calc_cv <- function(
        time_values,
        response_values) {
    time_data <- data.frame(
        Time = time_values,
        Response = response_values
    )
    
    unique_times <- unique(time_data$Time)
    
    cvs <- vapply(
        unique_times,
        function(t) {
            # Subset for the specific time point
            values_at_time <- time_data$Response[time_data$Time == t]
            if (mean(values_at_time, na.rm = TRUE) != 0 &&
                sum(!is.na(values_at_time)) > 1) {
                (sd(
                    values_at_time,
                    na.rm = TRUE
                ) /
                    mean(
                        values_at_time,
                        na.rm = TRUE
                    )) * 100
            } else {
                NA # Return NA for CV when mean is 0 or insufficient data points
            }
        },
        numeric(1)
    )
    # Return the average CV across time points
    return(mean(
        cvs,
        na.rm = TRUE
    ))
}


#' Generate Asterisks Definition HTML
#'
#' @noRd
#'
#' @description
#' Generate an HTML string defining the asterisk notation based on an adjusted
#' p-value threshold. Allows an optional title/header string and optional
#' inclusion of the "ns" (not significant) definition line.
#'
#' @param adj_pvalue_threshold Numeric. The adjusted p-value threshold.
#' @param header_html Character. Optional HTML shown above the definition lines
#'   (e.g., a <div> with centered title). If NULL, a default header is used.
#' @param include_ns Logical. Whether to include the "Adj. p-value > threshold
#'   --> ns" line. Default TRUE.
#' @param ns_label Character. Label appended after "ns". Default "ns".
#'
#' @return A character string containing the HTML definition for the asterisks.
#'
generate_asterisks_definition <- function(
        adj_pvalue_threshold,
        header_html = NULL,
        include_ns = TRUE,
        ns_label = "ns"
) {
    if (is.null(header_html)) {
        header_html <- paste(
            "<b><span style='font-size:20pt; margin-bottom: 0;'>",
            "Asterisks definition time effect:</span></b>",
            sep = ""
        )
    }
    
    lines <- character(0)
    
    if (isTRUE(include_ns)) {
        lines <- c(
            lines,
            paste(
                "Adj. p-value >",
                adj_pvalue_threshold,
                "-->",
                ns_label,
                sep = " "
            )
        )
    }
    
    lines <- c(
        lines,
        paste(
            "Adj. p-value <",
            adj_pvalue_threshold,
            "--> *",
            sep = " "
        ),
        paste(
            "Adj. p-value <",
            adj_pvalue_threshold / 5,
            "--> **",
            sep = " "
        ),
        paste(
            "Adj. p-value <",
            adj_pvalue_threshold / 50,
            "--> ***",
            sep = " "
        ),
        paste(
            "Adj. p-value <",
            adj_pvalue_threshold / 500,
            "--> ****",
            sep = " "
        )
    )
    
    paste(
        c(header_html, lines),
        collapse = "<br>"
    )
}


# Level 4 internal functions ---------------------------------------------------


#' Build title lines for a spline comparison plot
#'
#' @noRd
#'
#' @description
#' This internal helper function assembles the multi-line plot title used in
#' spline comparison visualizations. The title summarizes statistical
#' significance, effect sizes, relevance annotations, and average coefficient
#' of variation (CV) for the two compared conditions.
#'
#' @param feature_name A character scalar giving the feature identifier.
#' @param avrg_diff_pval A numeric scalar giving the adjusted p-value for the
#'   category-2 (average difference) test.
#' @param interaction_pval A numeric scalar giving the adjusted p-value for
#'   the category-3 (interaction) test.
#' @param adj_pthresh_avrg_diff_conditions Numeric. Adjusted p-value threshold
#'   for category-2 significance.
#' @param adj_pthresh_interaction Numeric. Adjusted p-value threshold for
#'   category-3 significance.
#' @param cat2_eff A numeric scalar giving the category-2 effect size, or
#'   `NA_real_` if unavailable.
#' @param es1 A numeric scalar giving the per-condition time effect size for
#'   `condition_1`, or `NA_real_`.
#' @param es2 A numeric scalar giving the per-condition time effect size for
#'   `condition_2`, or `NA_real_`.
#' @param diff_es A numeric scalar giving the category-3 interaction effect
#'   size between the two conditions, or `NA_real_`.
#' @param condition_1 A character scalar naming the first condition.
#' @param condition_2 A character scalar naming the second condition.
#' @param min_effect_size A named list of minimum effect-size thresholds,
#'   used to label effects as not relevant (`"nr"`).
#' @param cv_1 Numeric. Average coefficient of variation for `condition_1`.
#' @param cv_2 Numeric. Average coefficient of variation for `condition_2`.
#'
#' @return
#' A character vector of title lines to be collapsed with newline separators
#' when added to a ggplot object.
#'
#' @details
#' The title includes:
#' \itemize{
#'   \item The feature name.
#'   \item Adjusted p-values for category-2 and category-3 tests with
#'   asterisk or `"ns"` annotations.
#'   \item Category-2 and category-3 effect sizes, annotated with `"nr"`
#'   when below the corresponding effect-size thresholds.
#'   \item Average CV values for both conditions.
#' }
#' All numeric values are formatted for concise display, and missing values
#' are rendered as `"NA"`.
#'
.build_title_lines <- function(
        feature_name,
        avrg_diff_pval,
        interaction_pval,
        adj_pthresh_avrg_diff_conditions,
        adj_pthresh_interaction,
        cat2_eff,
        es1,
        es2,
        diff_es,
        condition_1,
        condition_2,
        min_effect_size,
        cv_1,
        cv_2
) {
    avrg_diff_stars <- .stars_from(
        avrg_diff_pval,
        adj_pthresh_avrg_diff_conditions
    )
    interaction_stars <- .stars_from(
        interaction_pval,
        adj_pthresh_interaction
    )
    
    title_lines <- c(
        feature_name,
        paste(
            "adj.P.Val avrg_diff_conditions:",
            .fmt_p_for_title(avrg_diff_pval),
            avrg_diff_stars
        ),
        paste(
            "adj.P.Val interaction_condition_time:",
            .fmt_p_for_title(interaction_pval),
            interaction_stars
        )
    )
    
    if (!is.na(cat2_eff)) {
        cat2_lab <- signif(cat2_eff, 3)
        cat2_nr <- !is.null(min_effect_size[["avg_diff_cond"]]) &&
            !is.na(min_effect_size[["avg_diff_cond"]]) &&
            abs(cat2_eff) < min_effect_size[["avg_diff_cond"]]
        
        title_lines <- c(
            title_lines,
            paste0(
                "Avrg diff conditions: ",
                cat2_lab,
                if (cat2_nr) " nr" else ""
            )
        )
    }
    
    if (!is.na(es1) || !is.na(es2) || !is.na(diff_es)) {
        ct1 <- if (is.na(es1)) {
            "NA"
        } else {
            as.character(signif(es1, 3))
        }
        
        ct2 <- if (is.na(es2)) {
            "NA"
        } else {
            as.character(signif(es2, 3))
        }
        
        cdt <- if (is.na(diff_es)) {
            "NA"
        } else {
            as.character(signif(diff_es, 3))
        }
        
        cdt_nr <- !is.null(min_effect_size[["interaction_cond_time"]]) &&
            !is.na(min_effect_size[["interaction_cond_time"]]) &&
            !is.na(diff_es) &&
            abs(diff_es) < min_effect_size[["interaction_cond_time"]]
        
        title_lines <- c(
            title_lines,
            paste0(
                "cT: ",
                condition_1,
                "=",
                ct1,
                " | ",
                condition_2,
                "=",
                ct2,
                " | cDT: ",
                cdt,
                if (cdt_nr) " nr" else ""
            )
        )
    }
    
    title_lines <- c(
        title_lines,
        paste0(
            "avg CV ",
            condition_1,
            ": ",
            round(cv_1, 2),
            "% | avg CV ",
            condition_2,
            ": ",
            round(cv_2, 2),
            "%"
        )
    )
    
    title_lines
}


#' Build a two-condition spline comparison ggplot
#'
#' @noRd
#'
#' @description
#' This internal helper function constructs the ggplot object for a spline
#' comparison between two conditions. It overlays observed data points and
#' fitted spline curves, optionally encodes replicate identity via point
#' shape, and configures legend entries for measured versus imputed data.
#'
#' @param plot_data A data frame produced by `.build_plot_data()`, containing
#'   `Time`, condition-specific response columns (`Y1`, `Y2`), and imputation
#'   flags (`IsImputed1`, `IsImputed2`). If `replicate_column` is used, it
#'   should also contain `Replicate`.
#'   
#' @param smooth_timepoints A numeric vector of time points used to draw
#'   fitted spline curves.
#'   
#' @param fitted_values_1 A numeric vector of fitted spline values for
#'   `condition_1`, aligned to `smooth_timepoints`.
#'   
#' @param fitted_values_2 A numeric vector of fitted spline values for
#'   `condition_2`, aligned to `smooth_timepoints`.
#'   
#' @param condition_1 A character scalar naming the first condition.
#' 
#' @param condition_2 A character scalar naming the second condition.
#' 
#' @param replicate_column Optional. A character scalar naming the replicate
#'   column in `plot_data`. If `NULL`, replicate shapes are not mapped.
#'   
#' @param shape_mapping Optional. A named integer vector mapping replicate
#'   identifiers to ggplot2 shape codes, typically produced by
#'   `.make_shape_mapping()`.
#'   
#' @param plot_info A named list of plotting metadata, including `time_unit`
#'   and `y_axis_label`.
#'   
#' @param plot_options `list`: Named list controlling optional plot
#'   customization. Supported entries include:
#'   \itemize{
#'     \item \code{cluster_heatmap_columns} `logical(1)`
#'       (default = \code{FALSE}):
#'       Whether to cluster columns in the heatmap.
#'
#'     \item \code{meta_replicate_column} `character(1)`
#'       (default = \code{NULL}):
#'       Column in \code{meta} encoding replicate information. When supplied,
#'       spline plot data points are colored by replicate.
#'
#'     \item \code{condition_colours} `list`
#'       (default = \code{NULL}):
#'       Optional named list mapping condition levels to colours. Names must
#'       correspond to values in the condition column of \code{meta}, and
#'       values must be valid R colour specifications accepted by
#'       \code{ggplot2}. When provided, these colours override the default
#'       colours for the corresponding conditions.
#'   }

#' @param title_lines A character vector of title lines, typically produced
#'   by `.build_title_lines()`, which will be collapsed with newlines.
#'   
#' @param has_imputed_1 Logical. Provided for compatibility; recomputed
#'   internally based on non-`NA` plotted points for `condition_1`.
#'   
#' @param has_imputed_2 Logical. Provided for compatibility; recomputed
#'   internally based on non-`NA` plotted points for `condition_2`.
#'
#' @return
#' A ggplot object showing observed points and fitted spline curves for both
#' conditions, with legends for data/spline elements and optional imputation
#' labels.
#'
#' @details
#' To prevent legend entries for imputed data when no imputed points are
#' actually plotted, the function subsets the input into condition-specific
#' point data (`plot_data_1`, `plot_data_2`) by removing rows with `NA`
#' responses. Imputation presence is then evaluated on these subsets and
#' used to control which legend labels are included in the manual color
#' scale.
#'
.build_comparison_plot <- function(
        plot_data,
        smooth_timepoints,
        fitted_values_1,
        fitted_values_2,
        condition_1,
        condition_2,
        replicate_column,
        shape_mapping,
        plot_info,
        plot_options,
        title_lines,
        has_imputed_1,
        has_imputed_2
) {
    plot_data$ColorLabel1 <- ifelse(
        plot_data$IsImputed1 == "Imputed",
        paste("Imputed data", condition_1),
        paste("Data", condition_1)
    )
    plot_data$ColorLabel2 <- ifelse(
        plot_data$IsImputed2 == "Imputed",
        paste("Imputed data", condition_2),
        paste("Data", condition_2)
    )
    
    plot_data_1 <- plot_data[!is.na(plot_data$Y1), , drop = FALSE]
    plot_data_2 <- plot_data[!is.na(plot_data$Y2), , drop = FALSE]
    
    has_imputed_1 <- any(plot_data_1$IsImputed1 == "Imputed")
    has_imputed_2 <- any(plot_data_2$IsImputed2 == "Imputed")
    
    p <- ggplot2::ggplot() +
        ggplot2::geom_point(
            data = plot_data_1,
            ggplot2::aes(
                x = .data$Time,
                y = .data$Y1,
                color = .data$ColorLabel1,
                shape = if (!is.null(replicate_column)) {
                    .data$Replicate
                } else {
                    NULL
                }
            ),
            na.rm = TRUE,
            alpha = 0.5
        ) +
        ggplot2::geom_line(
            data = data.frame(
                Time = smooth_timepoints,
                Fitted = fitted_values_1
            ),
            ggplot2::aes(
                x = .data$Time,
                y = .data$Fitted,
                color = paste("Spline", condition_1)
            )
        ) +
        ggplot2::geom_point(
            data = plot_data_2,
            ggplot2::aes(
                x = .data$Time,
                y = .data$Y2,
                color = .data$ColorLabel2,
                shape = if (!is.null(replicate_column)) {
                    .data$Replicate
                } else {
                    NULL
                }
            ),
            na.rm = TRUE,
            alpha = 0.5
        ) +
        ggplot2::geom_line(
            data = data.frame(
                Time = smooth_timepoints,
                Fitted = fitted_values_2
            ),
            ggplot2::aes(
                x = .data$Time,
                y = .data$Fitted,
                color = paste("Spline", condition_2)
            )
        ) +
        ggplot2::guides(
            color = ggplot2::guide_legend(title = NULL),
            shape = ggplot2::guide_legend(title = "Replicate")
        ) +
        ggplot2::scale_x_continuous(
            labels = scales::label_number_auto()
        ) +
        ggplot2::guides(
            x = ggplot2::guide_axis(check.overlap = TRUE)
        ) +
        ggplot2::labs(
            title = paste(title_lines, collapse = "\n"),
            x = paste0(
                "Time [",
                plot_info[["time_unit"]],
                "]"
            ),
            y = plot_info[["y_axis_label"]]
        )
    
    if (!is.null(replicate_column)) {
        p <- p + ggplot2::scale_shape_manual(
            values = shape_mapping,
            name = "Replicate"
        )
    }
    
    p <- p + ggplot2::theme_minimal() +
        ggplot2::theme(
            legend.position = "right",
            legend.title = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(size = 7),
            legend.text = ggplot2::element_text(size = 8),
            legend.key.height = ggplot2::unit(0.4, "cm"),
            legend.key.width = ggplot2::unit(0.8, "cm"),
            axis.title.x = ggplot2::element_text(size = 14),
            axis.title.y = ggplot2::element_text(size = 14),
            axis.text.x = ggplot2::element_text(size = 12),
            axis.text.y = ggplot2::element_text(size = 12)
        )
    
    y_combined <- c(plot_data$Y1, plot_data$Y2)
    y_max <- max(y_combined, na.rm = TRUE)
    y_min <- min(y_combined, na.rm = TRUE)
    y_extension <- (y_max - y_min) * 0.1
    y_pos_label <- y_max + y_extension * 0.5
    
    res <- maybe_add_dashed_lines(
        p = p,
        plot_info = plot_info,
        condition_1 = condition_1,
        condition_2 = condition_2,
        y_pos = y_pos_label,
        horizontal_labels = TRUE
    )
    p <- res$p
    treatment_colors <- res$treatment_colors
    
    color_values <- setNames(
        c(
            "orange",
            "orange",
            "purple",
            "purple",
            "red",
            "dodgerblue"
        ),
        c(
            paste("Data", condition_1),
            paste("Spline", condition_1),
            paste("Data", condition_2),
            paste("Spline", condition_2),
            paste("Imputed data", condition_1),
            paste("Imputed data", condition_2)
        )
    )

    cc <- plot_options[["condition_colours"]]
    
    if (!is.null(cc) &&
        is.list(cc) &&
        condition_1 %in% names(cc) &&
        is.character(cc[[condition_1]]) &&
        length(cc[[condition_1]]) == 1L) {
        color_values[paste("Data", condition_1)] <- cc[[condition_1]]
        color_values[paste("Spline", condition_1)] <- cc[[condition_1]]
    }
    
    if (!is.null(cc) &&
        is.list(cc) &&
        condition_2 %in% names(cc) &&
        is.character(cc[[condition_2]]) &&
        length(cc[[condition_2]]) == 1L) {
        color_values[paste("Data", condition_2)] <- cc[[condition_2]]
        color_values[paste("Spline", condition_2)] <- cc[[condition_2]]
    }
    
    filtered_labels <- c(
        paste("Data", condition_1),
        paste("Spline", condition_1),
        paste("Data", condition_2),
        paste("Spline", condition_2)
    )
    
    if (has_imputed_1) {
        filtered_labels <- c(
            filtered_labels,
            paste("Imputed data", condition_1)
        )
    }
    
    if (has_imputed_2) {
        filtered_labels <- c(
            filtered_labels,
            paste("Imputed data", condition_2)
        )
    }
    
    color_values <- c(
        color_values[
            names(color_values) %in% filtered_labels
        ],
        treatment_colors
    )
    
    p + ggplot2::scale_color_manual(
        values = color_values,
        drop = TRUE
    )
}


#' Conditionally add dashed lines for treatment timepoints
#'
#' @noRd
#'
#' @description
#' This internal function checks whether there are valid treatment
#' timepoints and labels in the `plot_info` list and, if so, adds
#' dashed vertical lines (and corresponding legend entries) to a
#' ggplot object.
#'
#' There are two modes of operation:
#' \itemize{
#'   \item \strong{Level-specific mode}: if `level` is not `NA`, the
#'     function looks for entries named `level` in
#'     `plot_info$treatment_labels` and `plot_info$treatment_timepoints`.
#'     If both are present and valid, a single set of dashed lines is
#'     added for that level.
#'
#'   \item \strong{Pairwise mode}: if `level` is `NA`, the function
#'     instead looks for entries named `condition_1` and `condition_2`
#'     (if provided) in `plot_info$treatment_labels` and
#'     `plot_info$treatment_timepoints`. For each condition that is
#'     present in both lists, its treatment lines are added. This can
#'     result in up to two sets of treatment lines (one per condition).
#' }
#'
#' In all cases, treatment labels and timepoints must be non-empty and
#' non-`NA`. When multiple sets are combined (pairwise mode), all
#' labels and timepoints are concatenated and passed together to
#' `add_dashed_lines()`. A color is assigned to each treatment label
#' using `scales::hue_pal()`, and the resulting named color vector is
#' returned for use in the plot's legend.
#'
#' @param p A ggplot object to which dashed lines and labels may be
#'   added.
#' @param plot_info A list containing treatment metadata. It is expected
#'   to contain:
#'   \itemize{
#'     \item `treatment_labels`: a named list or vector of labels.
#'     \item `treatment_timepoints`: a named list or vector of numeric
#'       timepoints.
#'   }
#'   When used in level-specific or pairwise mode, `treatment_labels`
#'   and `treatment_timepoints` must be named, with names matching
#'   `level`, `condition_1`, and/or `condition_2`.
#' @param level A character string or `NA`. When not `NA`, used to look
#'   up level-specific treatment labels and timepoints. When `NA`, the
#'   function instead checks `condition_1` and `condition_2`.
#' @param condition_1 A character string or `NA`. When `level` is `NA`,
#'   this is used as a name to look up treatment labels and timepoints
#'   for the first condition.
#' @param condition_2 A character string or `NA`. When `level` is `NA`,
#'   this is used as a name to look up treatment labels and timepoints
#'   for the second condition.
#' @param y_pos A numeric value specifying the y-axis position where
#'   the text labels should be placed.
#' @param horizontal_labels Logical flag indicating whether to place
#'   the treatment labels horizontally (`TRUE`) or vertically
#'   (`FALSE`, default).
#'
#' @return A list with components:
#'   \itemize{
#'     \item `p`: The ggplot object with possibly added dashed lines
#'       and labels.
#'     \item `treatment_colors`: A named character vector of colors
#'       used for the treatment labels. Names correspond to the label
#'       strings.
#'   }
#'
#' @importFrom scales hue_pal
#'
maybe_add_dashed_lines <- function(
        p,
        plot_info,
        level = NA,
        condition_1 = NA,
        condition_2 = NA,
        y_pos = 1,
        horizontal_labels = FALSE
        ) {
    treatment_colors <- c()
    
    labs <- plot_info$treatment_labels
    tps  <- plot_info$treatment_timepoints
    
    # If there is no treatment info at all, nothing to do
    if (is.null(labs) || is.null(tps)) {
        return(list(p = p, treatment_colors = treatment_colors))
    }
    
    # Helper to extract labels/timepoints for a single name
    extract_for_name <- function(nm) {
        if (is.null(nm) || is.na(nm)) {
            return(list(labels = character(0), tps = numeric(0)))
        }
        if (is.null(names(labs)) || is.null(names(tps))) {
            return(list(labels = character(0), tps = numeric(0)))
        }
        if (!(nm %in% names(labs)) || !(nm %in% names(tps))) {
            return(list(labels = character(0), tps = numeric(0)))
        }
        
        lbls <- as.character(labs[[nm]])
        tpts <- as.numeric(tps[[nm]])
        
        if (length(lbls) == 0L ||
            length(tpts) == 0L ||
            anyNA(lbls) ||
            anyNA(tpts)) {
            return(list(labels = character(0), tps = numeric(0)))
        }
        
        list(labels = lbls, tps = tpts)
    }
    
    chosen_labels <- character(0)
    chosen_tps    <- numeric(0)
    
    if (!is.na(level)) {
        # Level-specific mode: use only `level`
        res <- extract_for_name(level)
        chosen_labels <- res$labels
        chosen_tps    <- res$tps
    } else {
        # Pairwise mode: look at condition_1 and condition_2
        res1 <- extract_for_name(condition_1)
        res2 <- extract_for_name(condition_2)
        
        chosen_labels <- c(res1$labels, res2$labels)
        chosen_tps    <- c(res1$tps,    res2$tps)
    }
    
    # If nothing valid was found, return unchanged
    if (length(chosen_labels) == 0L ||
        length(chosen_tps) == 0L) {
        return(list(p = p, treatment_colors = treatment_colors))
    }
    
    treatment_colors <- scales::hue_pal()(length(chosen_labels))
    names(treatment_colors) <- chosen_labels
    
    p <- add_dashed_lines(
        p = p,
        treatment_timepoints = chosen_tps,
        treatment_labels = chosen_labels,
        y_pos = y_pos,
        horizontal_labels = horizontal_labels
    )
    
    list(
        p = p,
        treatment_colors = treatment_colors
    )
}



#' Test whether a result is statistically significant and relevant
#'
#' @noRd
#'
#' @description
#' This internal helper function evaluates whether a result satisfies both
#' statistical significance and effect-size relevance criteria. It is used
#' to gate plotting and reporting to features that pass both thresholds.
#'
#' @param pval A numeric scalar giving the adjusted p-value.
#' @param pthr A numeric scalar specifying the adjusted p-value significance
#'   threshold.
#' @param eff A numeric scalar giving the effect size.
#' @param effthr A numeric scalar specifying the minimum effect-size
#'   threshold required for relevance.
#'
#' @return
#' A logical scalar: `TRUE` if the result is both statistically significant
#' (`pval < pthr`) and relevant (`abs(eff) >= effthr`), otherwise `FALSE`.
#'
#' @details
#' Missing values in any of the inputs cause the function to return `FALSE`,
#' ensuring conservative behavior when gating downstream plots or summaries.
#'
.is_sig_and_rel <- function(
        pval,
        pthr,
        eff,
        effthr
) {
    isTRUE(!is.na(pval) && pval < pthr) &&
        isTRUE(
            !is.na(eff) &&
                !is.na(effthr) &&
                abs(eff) >= effthr
        )
}


#' Build plotting data for a two-condition spline comparison
#'
#' @noRd
#'
#' @description
#' This internal helper function constructs the long-format data frame used
#' for plotting observed values for two conditions across time. It also
#' marks imputed versus measured observations (when `raw_data` is provided)
#' and optionally attaches replicate identifiers for shape mapping.
#'
#' @param data A numeric matrix or data frame of processed values used for
#'   plotting, with features in rows and samples in columns.
#' @param raw_data Optional. A numeric matrix or data frame of raw (pre-
#'   imputation) values with the same shape as `data`. Missing entries are
#'   used to flag imputed observations. If `NULL`, all observations are
#'   treated as measured.
#' @param meta A data frame of sample metadata containing a numeric `Time`
#'   column and the grouping column specified by `condition`.
#' @param condition A character scalar naming the condition column in `meta`
#'   used to split samples into `condition_1` and `condition_2`.
#' @param replicate_column Optional. A character scalar naming the replicate
#'   column in `meta`. If `NULL`, replicate fields are not added.
#' @param replicate_mapping Optional. A named integer vector mapping
#'   replicate identifiers to indices, typically produced by
#'   `.make_replicate_mapping()`.
#' @param hit_index Integer. Row index of the feature in `data`/`raw_data`
#'   to extract for plotting.
#' @param condition_1 A character scalar giving the first condition level.
#' @param condition_2 A character scalar giving the second condition level.
#'
#' @return
#' A named list with three elements:
#' \describe{
#'   \item{plot_data}{A data frame containing columns `Time`, `Y1`, `Y2`,
#'   `IsImputed1`, and `IsImputed2`, and optionally `Replicate` and
#'   `ReplicateLabel`.}
#'   \item{has_imputed_1}{Logical indicating whether any imputed values are
#'   present for `condition_1` in the extracted row.}
#'   \item{has_imputed_2}{Logical indicating whether any imputed values are
#'   present for `condition_2` in the extracted row.}
#' }
#'
#' @details
#' The function creates `Y1` and `Y2` vectors by assigning the feature's
#' values to the matching condition and `NA` otherwise. When `raw_data` is
#' provided, imputation flags are derived from `NA` entries in the raw
#' matrix for the corresponding condition-specific sample columns.
#'
.build_plot_data <- function(
        data,
        raw_data,
        meta,
        condition,
        replicate_column,
        replicate_mapping,
        hit_index,
        condition_1,
        condition_2
) {
    time_points <- meta$Time
    row_values <- data[hit_index, ]
    
    if (!is.null(raw_data)) {
        columns_condition_1 <- which(
            meta[[condition]] == condition_1
        )
        columns_condition_2 <- which(
            meta[[condition]] == condition_2
        )
        
        na_indices_cond1 <- columns_condition_1[
            which(
                is.na(
                    raw_data[hit_index, columns_condition_1]
                )
            )
        ]
        na_indices_cond2 <- columns_condition_2[
            which(
                is.na(
                    raw_data[hit_index, columns_condition_2]
                )
            )
        ]
        
        plot_data <- data.frame(
            Time = time_points,
            Y1 = ifelse(
                meta[[condition]] == condition_1,
                row_values,
                NA
            ),
            Y2 = ifelse(
                meta[[condition]] == condition_2,
                row_values,
                NA
            ),
            IsImputed1 = ifelse(
                seq_along(row_values) %in% na_indices_cond1,
                "Imputed",
                "Measured"
            ),
            IsImputed2 = ifelse(
                seq_along(row_values) %in% na_indices_cond2,
                "Imputed",
                "Measured"
            )
        )
    } else {
        plot_data <- data.frame(
            Time = time_points,
            Y1 = ifelse(
                meta[[condition]] == condition_1,
                row_values,
                NA
            ),
            Y2 = ifelse(
                meta[[condition]] == condition_2,
                row_values,
                NA
            ),
            IsImputed1 = "Measured",
            IsImputed2 = "Measured"
        )
    }
    
    has_imputed_1 <- any(plot_data$IsImputed1 == "Imputed")
    has_imputed_2 <- any(plot_data$IsImputed2 == "Imputed")
    
    if (!is.null(replicate_column)) {
        plot_data$Replicate <- meta[[replicate_column]]
        plot_data$ReplicateLabel <- replicate_mapping[
            meta[[replicate_column]]
        ]
    }
    
    list(
        plot_data = plot_data,
        has_imputed_1 = has_imputed_1,
        has_imputed_2 = has_imputed_2
    )
}


#' Create a replicate-to-index mapping for plotting
#'
#' @noRd
#'
#' @description
#' This internal helper function generates a named integer vector that maps
#' replicate identifiers to consecutive numeric indices. The mapping is
#' used to provide stable shape assignments when plotting replicate-level
#' observations.
#'
#' @param meta A data frame containing sample-level metadata.
#' @param replicate_column A character scalar giving the name of the column
#'   in `meta` that identifies biological or technical replicates. If
#'   `NULL`, no mapping is created.
#'
#' @return
#' A named integer vector mapping replicate identifiers to indices, or
#' `NULL` if `replicate_column` is `NULL`.
#'
#' @details
#' Replicate identifiers are taken in the order returned by
#' `unique(meta[[replicate_column]])`. The mapping is intended for use with
#' discrete aesthetics such as point shapes in ggplot2.
#'
.make_replicate_mapping <- function(
        meta,
        replicate_column
) {
    if (is.null(replicate_column)) {
        return(NULL)
    }
    
    setNames(
        seq_along(unique(meta[[replicate_column]])),
        unique(meta[[replicate_column]])
    )
}


#' Create a shape mapping for replicate identifiers
#'
#' @noRd
#'
#' @description
#' This internal helper function assigns point shapes to replicate
#' identifiers for use in ggplot2 visualizations. A set of distinct shapes
#' is used first, followed by a fallback shape if the number of replicates
#' exceeds the distinct set.
#'
#' @param meta A data frame containing sample-level metadata.
#' @param replicate_column A character scalar giving the name of the column
#'   in `meta` that identifies biological or technical replicates. If
#'   `NULL`, no shape mapping is created.
#'
#' @return
#' A named integer vector mapping replicate identifiers to ggplot2 shape
#' codes, or `NULL` if `replicate_column` is `NULL`.
#'
#' @details
#' The function assigns shapes in the order of
#' `unique(meta[[replicate_column]])`. If the number of replicates exceeds
#' the predefined set of distinct shapes, the remaining replicates are
#' assigned a fallback shape to avoid errors.
#'
.make_shape_mapping <- function(
        meta,
        replicate_column
) {
    if (is.null(replicate_column)) {
        return(NULL)
    }
    
    distinct_shapes <- c(21, 22, 23, 24, 25, 3, 4, 8)
    fallback_shapes <- rep(1, 100)
    uniq_rep <- unique(meta[[replicate_column]])
    
    setNames(
        c(distinct_shapes, fallback_shapes)[seq_along(uniq_rep)],
        uniq_rep
    )
}


#' Extract category-2 (average difference) effect size for a feature
#'
#' @noRd
#'
#' @description
#' This internal helper function retrieves the category-2 effect size
#' (average difference between conditions) for a given feature from the
#' corresponding results table. If the feature is not present, a numeric
#' `NA` is returned.
#'
#' @param avrg_diff_conditions A data frame containing category-2 results,
#'   including a `feature_names` column and an effect-size column in the
#'   first position.
#' @param feature_name A character scalar giving the feature identifier
#'   to extract.
#'
#' @return
#' A single numeric value: the category-2 effect size for the feature, or
#' `NA_real_` if the feature is not present in `avrg_diff_conditions`.
#'
#' @details
#' The function assumes that the first column of `avrg_diff_conditions`
#' contains the effect size of interest. No assumptions are made about the
#' column name, allowing flexibility in upstream result formatting.
#'
.get_cat2_effect <- function(
        avrg_diff_conditions,
        feature_name
) {
    if (!feature_name %in% avrg_diff_conditions$feature_names) {
        return(NA_real_)
    }
    
    row_cat2 <- avrg_diff_conditions[
        avrg_diff_conditions$feature_names == feature_name,
        ,
        drop = FALSE
    ]
    
    if (nrow(row_cat2) == 0) {
        return(NA_real_)
    }
    
    as.numeric(row_cat2[[1]])
}


#' Extract category-3 (interaction) effect sizes for a feature
#'
#' @noRd
#'
#' @description
#' This internal helper function retrieves category-3 effect size information
#' for a given feature, including per-condition time effects (category 1)
#' and the interaction effect between two conditions. Missing values are
#' returned as numeric `NA`s to ensure safe downstream handling.
#'
#' @param predicted_timecurves A list-like object containing precomputed
#'   effect sizes, including `time_effect_effect_size` (per condition) and
#'   `interaction_effect_size` (per condition pair or a single vector).
#' @param condition_1 A character scalar giving the name of the first
#'   condition.
#' @param condition_2 A character scalar giving the name of the second
#'   condition.
#' @param feature_name A character scalar giving the feature identifier
#'   for which effect sizes should be extracted.
#'
#' @return
#' A named list with three numeric elements:
#' \describe{
#'   \item{es1}{Per-condition time effect size for `condition_1`, or
#'   `NA_real_` if unavailable.}
#'   \item{es2}{Per-condition time effect size for `condition_2`, or
#'   `NA_real_` if unavailable.}
#'   \item{diff_es}{Interaction (category-3) effect size between
#'   `condition_1` and `condition_2`, or `NA_real_` if unavailable.}
#' }
#'
#' @details
#' The function supports two interaction effect representations:
#' a named list indexed by condition-pair strings of the form
#' `"<cond1>_vs_<cond2>"` (with automatic fallback to the reversed order),
#' or a single unnamed vector for backward compatibility. All extracted
#' values are returned as numerics.
#'
.get_cat3_effects <- function(
        predicted_timecurves,
        condition_1,
        condition_2,
        feature_name
) {
    es1 <- NA_real_
    es2 <- NA_real_
    diff_es <- NA_real_
    
    es_list <- predicted_timecurves$time_effect_effect_size
    
    if (!is.null(es_list[[condition_1]])) {
        es1 <- unname(es_list[[condition_1]][feature_name])
    }
    
    if (!is.null(es_list[[condition_2]])) {
        es2 <- unname(es_list[[condition_2]][feature_name])
    }
    
    ies <- predicted_timecurves$interaction_effect_size
    
    if (!is.null(ies)) {
        if (is.list(ies)) {
            pair_name <- paste0(
                condition_1,
                "_vs_",
                condition_2
            )
            ies_vec <- ies[[pair_name]]
            
            if (is.null(ies_vec)) {
                pair_name_rev <- paste0(
                    condition_2,
                    "_vs_",
                    condition_1
                )
                ies_vec <- ies[[pair_name_rev]]
            }
        } else {
            ies_vec <- ies
        }
        
        if (!is.null(ies_vec)) {
            tmp <- unname(ies_vec[feature_name])
            
            if (length(tmp)) {
                diff_es <- tmp[1]
            }
        }
    }
    
    list(
        es1 = es1,
        es2 = es2,
        diff_es = diff_es
    )
}


#' Safely extract a numeric p-value for a feature
#'
#' @noRd
#'
#' @description
#' This internal helper function retrieves the adjusted p-value for a given
#' feature from a results table. It ensures that the returned value is always
#' numeric (`NA_real_` if unavailable), preventing downstream errors when
#' formatting or comparing p-values.
#'
#' @param tbl A data frame containing a column `feature_names` and a column
#'   specified by `col` (default `"adj.P.Val"`).
#' @param feature_name A character scalar giving the feature identifier
#'   to match against the `feature_names` column.
#' @param col A character scalar specifying the name of the column to pull
#'   (defaults to `"adj.P.Val"`).
#'
#' @return
#' A single numeric value: the p-value corresponding to the feature, or
#' `NA_real_` if the feature is not present or the value is missing.
#'
#' @details
#' The function performs a filtered lookup on `tbl` and extracts the first
#' matching value from column `col`. It explicitly converts the result to
#' numeric and returns `NA_real_` if the feature is absent or the value
#' is `NA`. This ensures safe use of p-values in subsequent numerical
#' operations such as significance testing, thresholding, or annotation.
#'
#' @importFrom dplyr filter pull
#'
safe_pull_pval <- function(
        tbl,
        feature_name,
        col = "adj.P.Val") {
    if (is.null(tbl) ||
        !is.data.frame(tbl) ||
        !(col %in% names(tbl))) {
        return(NA_real_)
    }
    v <- tbl |>
        dplyr::filter(.data$feature_names == feature_name) |>
        dplyr::pull({{ col }})
    if (length(v) == 0 || is.na(v[1])) NA_real_ else as.numeric(v[1])
}


#' Sort result tables by feature name for stable plotting
#'
#' @noRd
#'
#' @description
#' This internal helper function sorts the provided result tables by
#' `feature_names` to ensure stable and reproducible downstream behavior
#' (e.g., consistent plot ordering and matching across tables).
#'
#' @param time_effect_1 A data frame for the first condition containing a
#'   `feature_names` column.
#' @param time_effect_2 A data frame for the second condition containing a
#'   `feature_names` column.
#' @param avrg_diff_conditions A data frame containing `feature_names` for
#'   average-difference (category 2) results.
#' @param interaction_condition_time A data frame containing `feature_names`
#'   for interaction (category 3) results.
#'
#' @return
#' A named list with four elements: `time_effect_1`, `time_effect_2`,
#' `avrg_diff_conditions`, and `interaction_condition_time`, each sorted by
#' `feature_names`.
#'
#' @importFrom dplyr arrange
#'
.sort_inputs_for_plotting <- function(
        time_effect_1,
        time_effect_2,
        avrg_diff_conditions,
        interaction_condition_time
) {
    list(
        time_effect_1 = dplyr::arrange(
            time_effect_1,
            .data$feature_names
        ),
        time_effect_2 = dplyr::arrange(
            time_effect_2,
            .data$feature_names
        ),
        avrg_diff_conditions = dplyr::arrange(
            avrg_diff_conditions,
            .data$feature_names
        ),
        interaction_condition_time = dplyr::arrange(
            interaction_condition_time,
            .data$feature_names
        )
    )
}


#' Extract prediction matrices and time grid for two conditions
#'
#' @noRd
#'
#' @description
#' This internal helper function retrieves the spline prediction matrices
#' and the corresponding time grid for two specified conditions from a
#' `predicted_timecurves` object. The returned components are used for
#' plotting fitted spline curves alongside observed data.
#'
#' @param predicted_timecurves A list-like object containing a numeric vector
#'   `time_grid` and a named list `predictions`, where each element is a
#'   matrix of predicted values indexed by feature and time.
#' @param condition_1 A character scalar giving the name of the first
#'   condition to extract from `predicted_timecurves$predictions`.
#' @param condition_2 A character scalar giving the name of the second
#'   condition to extract from `predicted_timecurves$predictions`.
#'
#' @return
#' A named list with three elements: `smooth_timepoints` (the shared time
#' grid), `pred_mat_1` (prediction matrix for `condition_1`), and
#' `pred_mat_2` (prediction matrix for `condition_2`).
#'
#' @details
#' The function assumes that `condition_1` and `condition_2` are valid names
#' in `predicted_timecurves$predictions`. No additional validation is
#' performed; missing entries will propagate as `NULL` and should be
#' handled by the calling code.
#'
.get_prediction_mats <- function(
        predicted_timecurves,
        condition_1,
        condition_2
) {
    list(
        smooth_timepoints = predicted_timecurves$time_grid,
        pred_mat_1 = predicted_timecurves$predictions[[condition_1]],
        pred_mat_2 = predicted_timecurves$predictions[[condition_2]]
    )
}


#' Select a balanced set of features for spline comparison plots
#'
#' @noRd
#'
#' @description
#' This internal helper function delegates to `select_balanced_hits()` to
#' choose a limited, balanced set of features for plotting based on
#' category-2 (average difference) and category-3 (interaction) results.
#' Selection is driven by adjusted p-values and user-defined thresholds.
#'
#' @param avrg_diff_conditions A data frame containing category-2 results,
#'   including columns `feature_nr`, `feature_names`, and optionally
#'   `adj.P.Val`.
#' @param interaction_condition_time A data frame containing category-3
#'   results, including columns `feature_nr`, `feature_names`, and optionally
#'   `adj.P.Val`.
#' @param max_hit_number Integer. The maximum number of features to select
#'   for plotting.
#' @param adj_pthresh_avrg_diff_conditions Numeric. Adjusted p-value
#'   threshold used for category-2 results.
#' @param adj_pthresh_interaction Numeric. Adjusted p-value threshold used
#'   for category-3 interaction results.
#'
#' @return
#' A data frame with at least the columns `feature_nr` and `feature_names`,
#' containing the features selected for plotting.
#'
#' @details
#' The function performs no additional filtering beyond delegating to
#' `select_balanced_hits()`. It ensures that only the required columns are
#' passed downstream, providing a stable interface for the plotting
#' pipeline.
#'
#' @importFrom dplyr select any_of
#'
.select_features_to_plot <- function(
        avrg_diff_conditions,
        interaction_condition_time,
        max_hit_number,
        adj_pthresh_avrg_diff_conditions,
        adj_pthresh_interaction
) {
    select_balanced_hits(
        avrg_df = dplyr::select(
            avrg_diff_conditions,
            "feature_nr",
            "feature_names",
            dplyr::any_of("adj.P.Val")
        ),
        inter_df = dplyr::select(
            interaction_condition_time,
            "feature_nr",
            "feature_names",
            dplyr::any_of("adj.P.Val")
        ),
        max_n = max_hit_number,
        adj_pthresh_avrg_diff_conditions =
            adj_pthresh_avrg_diff_conditions,
        adj_pthresh_interaction =
            adj_pthresh_interaction
    )
}


#' Filter selected features to those present in prediction matrices
#'
#' @noRd
#'
#' @description
#' This internal helper function restricts a set of selected features to
#' those that are present in both prediction matrices used for plotting.
#' This prevents downstream subsetting errors when extracting fitted
#' spline values by feature name.
#'
#' @param features_to_plot A data frame containing at least a
#'   `feature_names` column identifying selected features.
#' @param pred_mat_1 A numeric matrix of predicted values for the first
#'   condition, with row names corresponding to feature identifiers.
#' @param pred_mat_2 A numeric matrix of predicted values for the second
#'   condition, with row names corresponding to feature identifiers.
#'
#' @return
#' A data frame containing only those rows of `features_to_plot` whose
#' `feature_names` occur as row names in both `pred_mat_1` and `pred_mat_2`.
#'
#' @details
#' If `pred_mat_1` has no row names, no filtering is applied and
#' `features_to_plot` is returned unchanged. The function assumes that
#' `pred_mat_1` and `pred_mat_2` use the same feature identifiers.
#'
.filter_features_in_prediction_mats <- function(
        features_to_plot,
        pred_mat_1,
        pred_mat_2
) {
    if (is.null(rownames(pred_mat_1))) {
        return(features_to_plot)
    }
    
    features_to_plot[
        features_to_plot$feature_names %in% rownames(pred_mat_1) &
            features_to_plot$feature_names %in% rownames(pred_mat_2),
        ,
        drop = FALSE
    ]
}


# Level 5 internal functions ---------------------------------------------------


#' Convert an adjusted p-value to an asterisk significance label
#'
#' @noRd
#'
#' @description
#' This internal helper function maps a numeric p-value to a conventional
#' asterisk-based significance label (`*`, `**`, `***`, `****`) or `"ns"`
#' (not significant) based on a supplied adjusted p-value threshold.
#'
#' @param pval A numeric scalar giving the adjusted p-value. If `NA`, an
#'   empty string is returned.
#' @param thresh A numeric scalar specifying the adjusted p-value threshold
#'   used to define significance levels.
#'
#' @return
#' A character scalar containing the significance label corresponding to
#' `pval`, or an empty string if `pval` is `NA`.
#'
#' @details
#' The mapping follows a standard convention:
#' \itemize{
#'   \item `pval < thresh` \eqn{\rightarrow} `*`
#'   \item `pval < thresh / 5` \eqn{\rightarrow} `**`
#'   \item `pval < thresh / 50` \eqn{\rightarrow} `***`
#'   \item `pval < thresh / 500` \eqn{\rightarrow} `****`
#'   \item otherwise \eqn{\rightarrow} `ns`
#' }
#' The function assumes `thresh` is positive and non-zero.
#'
.stars_from <- function(
        pval,
        thresh
) {
    if (is.na(pval)) {
        return("")
    }
    
    if (pval < thresh / 500) {
        "****"
    } else if (pval < thresh / 50) {
        "***"
    } else if (pval < thresh / 5) {
        "**"
    } else if (pval < thresh) {
        "*"
    } else {
        "ns"
    }
}


#' Add dashed lines for treatment timepoints to a plot
#'
#' @noRd
#'
#' @description
#' This internal function adds dashed vertical lines at specified
#' treatment timepoints to a plot, along with text labels that
#' display the corresponding x-axis values.
#'
#' @param p A ggplot object. The plot to which dashed lines and labels
#' will be added.
#' @param treatment_timepoints A numeric vector of timepoints where
#' dashed lines should be drawn.
#' @param treatment_labels A character vector of labels corresponding
#' to each treatment timepoint. These labels are used for coloring
#' the lines, but the x-axis coordinates are displayed as the labels.
#' @param y_pos A numeric value specifying the y-axis position where
#' the text labels should be placed.
#' @param horizontal_labels Boolean flag indicating whether to have a vertical
#' label (default) or horizontal label.
#'
#' @return A ggplot object with added dashed lines and labels.
#'
#' @importFrom ggplot2 geom_vline geom_text aes
#' @importFrom scales hue_pal
#' @importFrom rlang .data
#'
add_dashed_lines <- function(
        p,
        treatment_timepoints,
        treatment_labels,
        y_pos = 1,
        horizontal_labels = FALSE) {
    # Check if treatment labels and timepoints are valid
    if (!is.null(treatment_timepoints) &&
        !is.null(treatment_labels) &&
        all(!is.na(treatment_timepoints)) &&
        all(!is.na(treatment_labels))) {
        # Create a data frame for the treatment lines
        treatment_df <- data.frame(
            Time = treatment_timepoints,
            Label = treatment_labels,
            y_pos = y_pos
        )
        
        # Generate distinct colors for the treatment labels
        treatment_colors <- scales::hue_pal()(length(treatment_labels))
        names(treatment_colors) <- treatment_labels
        
        # Add dashed vertical lines and text labels to the plot
        p <- p +
            ggplot2::geom_vline(
                data = treatment_df,
                ggplot2::aes(
                    xintercept = .data$Time,
                    color = .data$Label
                ),
                linetype = "dashed",
                linewidth = 0.5
            ) +
            ggplot2::geom_text(
                data = treatment_df,
                ggplot2::aes(
                    x = ifelse(
                        horizontal_labels,
                        .data$Time + max(treatment_timepoints) * 0.04,
                        .data$Time - max(treatment_timepoints) * 0.005
                    ),
                    y = .data$y_pos,
                    label = round(.data$Time, 2),
                    color = .data$Label
                ),
                angle = if (horizontal_labels) 0 else 90,
                vjust = if (horizontal_labels) -0.2 else 0,
                hjust = if (horizontal_labels) 0.5 else 1,
                size = 3,
                show.legend = FALSE
            )
    }
    
    return(p) # Return the updated plot object
}


#' Format a p-value for inclusion in plot titles
#'
#' @noRd
#'
#' @description
#' This internal helper function formats a numeric p-value for display in
#' plot titles. Finite values are rounded to two significant digits, while
#' missing values are represented as `"NA"`.
#'
#' @param p A numeric scalar giving a p-value.
#'
#' @return
#' A character scalar suitable for inclusion in plot titles: a formatted
#' numeric string for non-missing values, or `"NA"` if `p` is `NA`.
#'
#' @details
#' The function performs no significance assessment; it only controls
#' textual formatting. Significance annotations (e.g., asterisks or `"ns"`)
#' are handled separately.
#'
.fmt_p_for_title <- function(
        p
) {
    if (is.na(p)) {
        "NA"
    } else {
        as.character(signif(p, 2))
    }
}
