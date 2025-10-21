#' run_ora()
#'
#' @description
#' This function generates a overrepresentation analysis report based
#' on clustered hit levels, gene data, and specified databases. It accomplishes
#' this by using the R package clusterProfiler. As output, you will receive a
#' list of the plot objects it generated, and an HTML report with embedded
#' files containing the enrichment results, and dotplots visualizing the
#' enrichment.
#'
#' @param cluster_table A tibble containing one row per
#'   \code{feature_nr} with metadata and cluster assignments across the
#'   analysis categories. It includes:
#'   \itemize{
#'     \item \code{feature_nr} - Numeric feature identifier.
#'     \item \code{feature_name} - Preferred feature name from the source
#'       data, falling back to the numeric ID if none is available.
#'     \item \code{gene} - Preferred gene symbol from the annotation or
#'       cluster data.
#'     \item \code{cluster_<cond1>} / \code{cluster_<cond2>} - Cluster
#'       assignments for each time-effect condition.
#'     \item \code{cluster_cat2} - (Optional) Combined cluster label for
#'       category 2 hits in the form
#'       \code{"<cluster_<cond1>>_<cluster_<cond2>>"}; \code{NA} if the
#'       feature was not a category 2 hit.
#'     \item \code{cluster_cat3} - (Optional) Combined cluster label for
#'       category 3 hits in the form
#'       \code{"<cluster_<cond1>>_<cluster_<cond2>>"}; \code{NA} if the
#'       feature was not a category 3 hit.
#'   }
#'   For any category-specific cluster column, a value of \code{NA}
#'   indicates that the feature was not significant (not a hit) in that
#'   category.
#'
#' @param databases A \code{data.frame} that defines the gene set collections
#'   to be tested in the overrepresentation analysis. Must contain exactly
#'   three columns:
#'   \describe{
#'     \item{DB}{Character. The database identifier (e.g., KEGG, GO_BP,
#'       Reactome).}
#'     \item{Geneset}{Character. The name of the gene set or pathway within the
#'       database.}
#'     \item{Gene}{Character. A gene identifier belonging to the gene set
#'       (e.g., gene symbol, Ensembl ID).}
#'   }
#'
#'   Each row corresponds to one `(database, geneset, gene)` association. The
#'   same gene may appear in multiple gene sets.
#'
#' @param report_info A list containing information for the report generation,
#' such as omics_data_type and data_description (this is the list used for all
#' report generating functions of this package).
#'
#' @param cluster_hits_report_name Single character string specifying the name
#' of the cluster_hits() function report, that contains the results that were
#' used for the overprepresentation analysis here. Must be specified, because
#' otherwise, the connection is not documented.
#'
#' @param clusterProfiler_params A named list of arguments passed directly to
#'   the corresponding functions in the \strong{clusterProfiler} package.
#'   Typical entries include \code{pvalueCutoff}, \code{pAdjustMethod},
#'   \code{minGSSize}, \code{maxGSSize}, and \code{qvalueCutoff}. The names
#'   must match the argument names in clusterProfiler; see the clusterProfiler
#'   documentation for details. If \code{NULL} (default), the standard
#'   clusterProfiler defaults are used.
#'
#' @param mapping_cfg A named list that controls the optional behavior of
#'        automatically mapping gene symbols across species. This is useful
#'        when your input gene symbols (e.g., from CHO cells) do not match
#'        the species used by the enrichment databases (e.g., human or mouse).
#'        By default, no mapping is performed and gene symbols are used as-is.
#'        If mapping is desired, this list must contain the following three
#'        elements:
#'        \describe{
#'          \item{method}{Mapping method to use. One of \code{none} (default;
#'          no mapping), \code{gprofiler} (online, via the g:Profiler API), or
#'          \code{orthogene} (offline, if installed).}
#'          \item{from_species}{Source species code,
#'          e.g. \code{cgriseus} for CHO. Must match the expected format for
#'          the selected tool.}
#'          \item{to_species}{Target species code,
#'          e.g. \code{hsapiens} for human. This must be the species used in
#'          your ORA database and must also match the expected format for the
#'          selected tool.}
#'        }
#'
#' @param enrichGO_cfg A named list specifying the configuration for
#' running GO enrichment with Bioconductor's
#' \code{\link[clusterProfiler]{enrichGO}}.
#' This is only needed when you want to perform GO Biological Process (BP),
#' Molecular Function (MF), or Cellular Component (CC) enrichment using
#' Bioconductor's organism databases (e.g., \code{org.Mm.eg.db} for mouse).
#'
#' The list must be named according to the GO ontology, e.g., \code{"GO_BP"},
#' \code{"GO_MF"}, \code{"GO_CC"}. Each entry must provide:
#' \itemize{
#'   \item \code{OrgDb}: The organism database, e.g., \code{org.Mm.eg.db}.
#'   \item \code{keyType}: The gene identifier type, e.g., \code{"SYMBOL"}.
#'   \item \code{ontology}: One of \code{"BP"}, \code{"MF"}, or \code{"CC"}.
#' }
#'
#' If \code{enrichGO_cfg} is \code{NULL} (default), no Bioconductor-based GO
#' enrichment is performed. All enrichment runs through
#' \code{\link[clusterProfiler]{enricher}} with the provided TERM2GENE mappings.
#'
#' @param universe Enrichment background data. This is a
#' parameter of clusterProfiler, for the documentation, please check the
#' documentation of the clusterProfiler R package.
#'
#' @param report_dir Character string specifying the directory path where the
#' HTML report and any other output files should be saved. When no path is
#' specified, then the function runs but no HTML report is generated.
#'
#' @param verbose Boolean flag controlling the display of messages.
#'
#' @return
#' A nested, named list whose top-level elements correspond to the
#' limma result categories. The exact set of elements depends on \code{mode}:
#'
#' \describe{
#'   \item{\code{mode == "isolated"}}{
#'     Two elements are returned, one per condition level:
#'     \code{time_effect_condition_<level1>} and
#'     \code{time_effect_condition_<level2>}.
#'   }
#'   \item{\code{mode == "integrated"}}{
#'     The two time-effect elements above, plus (only if there are significant
#'     hits at the chosen thresholds) up to two additional elements:
#'     \code{avrg_diff_conditions} and \code{interaction_condition_time}.
#'     Note that the clusters of \code{interaction_condition_time} are 
#'     "combo-clusters" made of the cluster membership of the feature in
#'     condition 1 and the membership of the same feature in condition 2 (see 
#'     also the respective documentation for the function cluster_hits() about
#'     the cluster_table). For example, if the report generated by this function
#'     has the entries 'time_effect_condition_control' and 
#'     'time_effect_condition_treatment', and the section 
#'     'interaction_condition_time' contains entries such as 'cluster_4_2', then
#'     the first number (4) is the cluster of condition control, and the second
#'     number (2) is the cluster of condition treatment.
#'   }
#' }
#'
#' Each top-level result category element is a list with the fields:
#' \describe{
#'   \item{\code{dotplot}}{A \code{ggplot} object: the dot plot of
#'     over-representation results (clusterProfiler) for that category.}
#'   \item{\code{dotplot_nrows}}{Numeric scalar giving a suggested plot height
#'     (in rows / relative units) that prints nicely for the number of
#'     enriched terms shown.}
#'   \item{\code{ora_results}}{A nested list of the raw enrichment results,
#'     structured as:
#'     \describe{
#'       \item{\emph{cluster} \eqn{\rightarrow} \emph{database}}{
#'         For each cluster in the category, there is a sublist with one entry
#'         per database used in the enrichment. The value of each entry is
#'         either
#'         \code{NA} (no terms enriched for that cluster-database) or a
#'         \code{data.frame} as returned by
#'         \code{clusterProfiler::enricher()} for the enriched terms.}
#'     }}
#' }
#'
#' In summary, the full shape is:
#' \preformatted{
#' list(
#'   time_effect_condition_<level1> = list(
#'     dotplot        = ggplot,
#'     dotplot_nrows  = numeric(1),
#'     ora_results    = list(
#'       <cluster_1> = list(<database_1> = NA|data.frame, ...),
#'       <cluster_2> = list(<database_1> = NA|data.frame, ...),
#'       ...
#'     )
#'   ),
#'   time_effect_condition_<level2> = list(...),
#'   avrg_diff_conditions          = list(...),
#'   interaction_condition_time    = list(...)
#' )
#' }
#'
#' @seealso \code{clusterProfiler::enricher()}
#'
#' @importFrom purrr map2 flatten
#' @importFrom here here
#'
#' @examplesIf requireNamespace("clusterProfiler", quietly = TRUE)
#' {
#'     set.seed(1)
#'
#'     # --- toy cluster table (two "conditions") ------------------------------
#'     toy_genes <- paste0("G", 1:8)
#'     cluster_table <- tibble::tibble(
#'         feature_nr    = 1:8,
#'         feature_name  = paste0("feat_", 1:8),
#'         gene          = toy_genes,
#'         cluster_condA = c(1, 1, 2, 2, NA, NA, 1, 2),
#'         cluster_condB = c(NA, 1, NA, 2, 1, 2, 1, NA)
#'     )
#'
#'     # --- toy TERM2GENE database -------------------------------------------
#'     databases <- data.frame(
#'         DB = rep("ToyDB", 6),
#'         Geneset = c(rep("SetA", 3), rep("SetB", 3)),
#'         Gene = c("G1", "G2", "G7", "G3", "G4", "G6"),
#'         stringsAsFactors = FALSE
#'     )
#'
#'     # --- minimal report info ----------------------------------------------
#'     report_info <- list(
#'         omics_data_type = "TOY",
#'         data_description = "Toy dataset for run_ora() example",
#'         data_collection_date = "2025",
#'         analyst_name = "Example Analyst",
#'         contact_info = "analyst@example.org",
#'         project_name = "ToyProject"
#'     )
#'
#'     # --- output directory (temp) -------------------------------------------
#'     report_dir <- file.path(tempdir(), "run_ora_demo")
#'     dir.create(report_dir, showWarnings = FALSE, recursive = TRUE)
#'
#'     # --- permissive params for tiny example --------------------------------
#'     clusterProfiler_params <- list(
#'         pvalueCutoff = 1,
#'         qvalueCutoff = 1,
#'         minGSSize    = 1,
#'         maxGSSize    = 500
#'     )
#'
#'     # --- run ORA -----------------------------------------------------------
#'     res <- run_ora(
#'         cluster_table            = cluster_table,
#'         databases                = databases,
#'         report_info              = report_info,
#'         cluster_hits_report_name = "cluster_hits_demo",
#'         clusterProfiler_params   = clusterProfiler_params,
#'         report_dir               = report_dir,
#'         verbose                  = TRUE
#'     )
#'
#'     # see sections and files written
#'     names(res)
#'     list.files(report_dir, recursive = TRUE)
#' }
#'
#' @export
#'
run_ora <- function(
    cluster_table,
    databases,
    report_info,
    cluster_hits_report_name,
    clusterProfiler_params = NA,
    mapping_cfg = list(
        method = "none",
        from_species = NULL,
        to_species = NULL
    ),
    enrichGO_cfg = NULL,
    universe = NULL,
    report_dir = NULL,
    verbose = TRUE) {
    args <- lapply(
        as.list(
            match.call()[-1]
        ),
        eval,
        parent.frame()
    )
    args[["verbose"]] <- verbose
    input_control <- InputControl$new(args)
    input_control$auto_validate()

    # Control the test not covered by the InputControl class
    control_inputs_run_ora(
        cluster_table = cluster_table,
        databases = databases,
        params = clusterProfiler_params,
        mapping_cfg = mapping_cfg,
        enrichGO_cfg = enrichGO_cfg,
        background = universe,
        cluster_hits_report_name = cluster_hits_report_name,
        verbose = verbose
    )

    ensure_clusterProfiler() # Deals with clusterProfiler installation.

    cluster_table[["gene"]] <- map_gene_symbols(
        genes = cluster_table[["gene"]],
        mapping_cfg = mapping_cfg
    )
    if (!is.null(universe)) {
        if (!is.character(universe)) {
            stop(
                "'universe' must be a character vector of gene identifiers, ",
                "not an object of class '", class(universe)[1], "'."
            )
        }
        universe <- map_gene_symbols(
            genes = universe,
            mapping_cfg = mapping_cfg
        )
    }

    lvl_cols <- level_columns(cluster_table)
    lvl_names <- level_display_names(cluster_table)

    all_results <- mapply(function(col, nm) {
        manage_ora_result_cat(
            cluster_table = cluster_table,
            level_col = col,
            databases = databases,
            clusterProfiler_params = clusterProfiler_params,
            enrichGO_cfg = enrichGO_cfg,
            universe = universe,
            plot_title = nm, # use display name for the plot title
            verbose = verbose
        )
    }, col = lvl_cols, nm = lvl_names, SIMPLIFY = FALSE)
    names(all_results) <- lvl_names

    processed_results <- Map(
        function(res, nm) process_result(res, nm),
        all_results, names(all_results)
    )

    # Build a structured list of sections for the report
    sections <- lapply(seq_along(processed_results), function(i) {
        pr <- processed_results[[i]]
        raw <- all_results[[i]] # <- has $ora_results from run_ora_level()

        plots <- pr$plot
        if (inherits(plots, "ggplot") || !is.list(plots)) plots <- list(plots)

        sizes <- pr$plot_size
        if (!is.list(sizes)) sizes <- as.list(rep_len(sizes, length(plots)))

        list(
            level = names(processed_results)[i],
            header_name =
                pr$header_info$header_name %||% names(processed_results)[i],
            plots = plots,
            plot_sizes = sizes,
            ora_results = if (is.list(raw) && !is.null(raw$ora_results)) {
                raw$ora_results
            } else {
                NULL
            }
        )
    })
    names(sections) <- names(processed_results)

    # Update report_info payload
    report_info$databases <- unique(databases[["DB"]])
    report_info$clusterProfiler_params <- clusterProfiler_params
    report_info$cluster_hits_report_name <- cluster_hits_report_name
    report_info$cluster_table <- cluster_table
    report_info$background_genes <- universe

    # If *all* sections have zero plots, return NULL
    total_plots <- sum(
        vapply(sections, function(s) length(s$plots), integer(1))
        )
    if (total_plots == 0L && verbose) {
        message("No results --> Not generating a report and returning NULL.")
        return(NULL)
    }

    if (!is.null(report_dir)) {
        report_dir <- normalizePath(
            report_dir,
            mustWork = FALSE
        )
        generate_report_html(
            plots       = sections,
            report_info = report_info,
            report_type = "run_ora",
            filename    = "run_ora_report",
            report_dir  = report_dir
        )

        if (verbose) {
            print_info_message(
                message_prefix = "ORA analysis",
                report_dir = report_dir
            )
        }
    }

    return(all_results)
}



# Level 1 internal functions ---------------------------------------------------


#' Control Inputs for ora Report
#'
#' @noRd
#'
#' @description
#' Validates the inputs for generating a ora report, including clustered
#' hits, genes, databases, parameters, plot titles, and background genes.
#'
#' @param cluster_table A tibble containing one row per
#'   \code{feature_nr} with metadata and cluster assignments across the
#'   analysis categories. It includes:
#'   \itemize{
#'     \item \code{feature_nr} - Numeric feature identifier.
#'     \item \code{feature_name} - Preferred feature name from the source
#'       data, falling back to the numeric ID if none is available.
#'     \item \code{gene} - Preferred gene symbol from the annotation or
#'       cluster data.
#'     \item \code{cluster_<cond1>} / \code{cluster_<cond2>} - Cluster
#'       assignments for each time-effect condition.
#'     \item \code{cluster_cat2} - (Optional) Combined cluster label for
#'       category 2 hits in the form
#'       \code{"<cluster_<cond1>>_<cluster_<cond2>>"}; \code{NA} if the
#'       feature was not a category 2 hit.
#'     \item \code{cluster_cat3} - (Optional) Combined cluster label for
#'       category 3 hits in the form
#'       \code{"<cluster_<cond1>>_<cluster_<cond2>>"}; \code{NA} if the
#'       feature was not a category 3 hit.
#'   }
#'   For any category-specific cluster column, a value of \code{NA}
#'   indicates that the feature was not significant (not a hit) in that
#'   category.
#' @param databases A list of databases to be used in the ora analysis.
#' @param params A list of parameters for the ora analysis.
#' @param mapping_cfg A named list that controls the optional behavior of
#'        automatically mapping gene symbols across species. This is useful
#'        when your input gene symbols (e.g., from CHO cells) do not match
#'        the species used by the enrichment databases (e.g., human or mouse).
#'        By default, no mapping is performed and gene symbols are used as-is.
#'        If mapping is desired, this list must contain the following **three**
#'        elements:
#'        \describe{
#'          \item{method}{Mapping method to use. One of `"none"` (default;
#'          no mapping), `"gprofiler"` (online, via the g:Profiler API), or
#'           `"orthogene"` (offline, if installed).}
#'          \item{from_species}{Source species code
#'          (e.g., `"cgriseus"` for CHO). Must match the expected format for
#'           the selected tool.}
#'          \item{to_species"}{Target species code
#'          (e.g., `"hsapiens"` for human). This must be the species used in
#'           your ORA database.}
#'        }
#'
#' @param enrichGO_cfg A named list specifying the configuration for
#' running GO enrichment with Bioconductor's
#' \code{\link[clusterProfiler]{enrichGO}}.
#' This is only needed when you want to perform GO Biological Process (BP),
#' Molecular Function (MF), or Cellular Component (CC) enrichment using
#' Bioconductor's organism databases (e.g., \code{org.Mm.eg.db} for mouse).
#'
#' The list must be named according to the GO ontology, e.g., \code{"GO_BP"},
#' \code{"GO_MF"}, \code{"GO_CC"}. Each entry must provide:
#' \itemize{
#'   \item \code{OrgDb}: The organism database, e.g., \code{org.Mm.eg.db}.
#'   \item \code{keyType}: The gene identifier type, e.g., \code{"SYMBOL"}.
#'   \item \code{ontology}: One of \code{"BP"}, \code{"MF"}, or \code{"CC"}.
#' }
#'
#' If \code{enrichGO_cfg} is \code{NULL} (default), no Bioconductor-based GO
#' enrichment is performed. All enrichment runs through
#' \code{\link[clusterProfiler]{enricher}} with the provided TERM2GENE mappings.
#'
#' @param plot_titles A character vector of titles for the plots, with length
#' matching `levels_clustered_hits`.
#' @param background A character vector of background genes or NULL.
#' @param cluster_hits_report_name Single character string specifying the name
#' of the cluster_hits() function report, that contains the results that were
#' used for the overprepresentation analysis here. Must be specified, because
#' otherwise, the connection is not documented.
#' @param verbose Boolean flag controlling the display of messages.
#'
control_inputs_run_ora <- function(
    cluster_table,
    databases,
    params,
    mapping_cfg,
    enrichGO_cfg,
    background,
    cluster_hits_report_name,
    verbose) {
    check_cluster_table(cluster_table)

    check_databases(databases)

    check_params(params)

    check_mapping_cfg(
        mapping_cfg = mapping_cfg,
        verbose = verbose
    )

    check_enrichGO_cfg(enrichGO_cfg)

    if (!is.null(background)) {
        if (!is.character(background)) {
            stop_call_false("background must be a character vector or NULL.")
        } else {
            check_genes(background)
        }
    }

    if (!is.character(cluster_hits_report_name) ||
        length(cluster_hits_report_name) != 1) {
        stop_call_false(
            "`cluster_hits_report_name` must be a single character string."
        )
    }
}


#' Ensure 'clusterProfiler' is installed and loaded
#'
#' @noRd
#'
#' @description
#' This function checks if the 'clusterProfiler' package is installed.
#' If not, it prompts the user to choose whether to install it automatically,
#' install it manually, or cancel the operation. Once installed, the package
#' is loaded for use.
ensure_clusterProfiler <- function() {
    ok <- try(requireNamespace("clusterProfiler", quietly = TRUE),
              silent = TRUE)
    if (isTRUE(ok)) return(invisible(TRUE))
    
    is_installed <- "clusterProfiler" %in%
        rownames(utils::installed.packages())
    
    if (!is_installed) {
        stop_call_false(
            "The 'clusterProfiler' package is not installed.\n",
            "Please install it manually using:\n\n",
            "  if (!requireNamespace('BiocManager', quietly = TRUE)) ",
            "install.packages('BiocManager')\n",
            "  BiocManager::install('clusterProfiler')\n\n",
            "This is an optional dependency of the SplineOmics package."
        )
    }
    
    ns_err <- try(loadNamespace("clusterProfiler"), silent = TRUE)
    msg <- if (inherits(ns_err, "try-error")) {
        conditionMessage(attr(ns_err, "condition"))
    } else {
        "Unknown error."
    }
    
    stop_call_false(
        "The 'clusterProfiler' package is installed but failed to load.\n\n",
        "Underlying error:\n  ", msg, "\n\n",
        "Hints:\n",
        "  - Update/install dependencies with BiocManager::install(...)\n",
        "  - Check .libPaths() and R version; rebuild from source if needed\n",
        "  - Run BiocManager::valid() to fix mismatches"
    )
}


#' Manage ORA analysis for a summary level
#'
#' @noRd
#'
#' @description
#' Runs over-representation analysis (ORA) for a single level column from
#' the unified \code{cluster_table} table. Foreground gene sets are
#' built by grouping non-missing \code{gene} values by the labels in the
#' chosen level column (via \code{make_clustered_genes()}). For cat2 this
#' yields at most two clusters (\code{"<cond1>_higher"},
#' \code{"<cond2>_higher"}); other levels use their distinct cluster labels.
#' The resulting foregrounds are passed to \code{run_ora_level()}.
#'
#' @param cluster_table A tibble containing one row per
#'   \code{feature_nr} with metadata and cluster assignments across the
#'   analysis categories. It includes:
#'   \itemize{
#'     \item \code{feature_nr} - Numeric feature identifier.
#'     \item \code{feature_name} - Preferred feature name from the source
#'       data, falling back to the numeric ID if none is available.
#'     \item \code{gene} - Preferred gene symbol from the annotation or
#'       cluster data.
#'     \item \code{cluster_<cond1>} / \code{cluster_<cond2>} - Cluster
#'       assignments for each time-effect condition.
#'     \item \code{cluster_cat2} - (Optional) Combined cluster label for
#'       category 2 hits in the form
#'       \code{"<cluster_<cond1>>_<cluster_<cond2>>"}; \code{NA} if the
#'       feature was not a category 2 hit.
#'     \item \code{cluster_cat3} - (Optional) Combined cluster label for
#'       category 3 hits in the form
#'       \code{"<cluster_<cond1>>_<cluster_<cond2>>"}; \code{NA} if the
#'       feature was not a category 3 hit.
#'   }
#'   For any category-specific cluster column, a value of \code{NA}
#'   indicates that the feature was not significant (not a hit) in that
#'   category.
#' @param level_col Character scalar naming the level column in
#'   \code{cluster_table} to analyze (must exist).
#' @param databases Gene set collections (data frame or list) consumed by
#'   downstream helpers inside \code{run_ora_level()}.
#' @param clusterProfiler_params List of clusterProfiler parameters (e.g.,
#'   \code{pvalueCutoff}, \code{pAdjustMethod}, \code{minGSSize}).
#' @param enrichGO_cfg Optional named list with GO settings (e.g.,
#'   \code{OrgDb}, \code{keyType}, \code{ontology}) used when the selected
#'   collection is GO.
#' @param universe Optional character vector of background genes to be used
#'   as the ORA universe.
#' @param plot_title String specifying the title of the plot.
#' @param verbose Boolean flag controlling the display of messages.
#'
#' @return A list as returned by \code{run_ora_level()} with elements
#'   \code{dotplot}, \code{dotplot_nrows}, and \code{ora_results}. If no
#'   enrichment was obtained for any cluster, \code{NA} may be returned.
#'
manage_ora_result_cat <- function(
    cluster_table,
    level_col,
    databases,
    clusterProfiler_params,
    enrichGO_cfg,
    universe,
    plot_title,
    verbose) {
    cg <- make_clustered_genes(
        cluster_table,
        level_col
    )
    if (verbose) {
        message(paste0("\n\n\n Running clusterProfiler for: ", level_col))
    }

    run_ora_level(
        clustered_genes = cg,
        databases = databases,
        params = clusterProfiler_params,
        enrichGO_cfg = enrichGO_cfg,
        plot_title = plot_title,
        universe = universe,
        verbose = verbose
    )
}


#' Process ORA Result for a Specific Level
#'
#' @noRd
#'
#' @description
#' This function processes the ORA result for a specific level. It handles
#' cases where the result contains `NA` values by adding a section break.
#' Otherwise, it extracts the plot, plot size, and header information from
#' the result.
#'
#' @param level_result A list containing the ORA result for a specific level.
#' @param level_name A character string representing the name of the level.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{plot}{A plot object or "section_break" if the result contains `NA`.}
#'   \item{plot_size}{An integer indicating the size of the plot.}
#'   \item{header_info}{A list with header information, including the level
#'   name, full enrichment results, and raw enrichment results if available.}
#' }
#'
process_result <- function(
    level_result,
    level_name) {
    result <- list()

    if (any(is.na(level_result))) {
        result$plot <- "section_break"
        result$plot_size <- 999
        result$header_info <- list(
            header_name = level_name,
            full_enrich_results = NA
        )
    } else {
        result$plot <- list(level_result$dotplot)
        result$plot_size <- level_result$dotplot_nrows
        result$header_info <- list(
            header_name = level_name,
            ora_results = level_result$ora_results
        )
    }

    return(result)
}


#' Build ORA report
#'
#' @noRd
#'
#' @description
#' Generates an HTML report for over-representation analysis (ORA) using a
#' structured list of sections. Each section may include one or more plots,
#' optional per-section ORA results (to render a summary table and a download
#' link), and a human-readable header. A table of contents is created and the
#' final HTML is written to disk.
#'
#' @param header_section Character scalar with HTML for the report header that
#'   appears before the table of contents.
#' @param sections A named list of sections. Each element is a list with fields:
#'   \itemize{
#'     \item \code{header_name} (character): Section title.
#'     \item \code{plots} (list): One or more plot objects (e.g. ggplot) or
#'       renderable items for \code{process_plots()}.
#'     \item \code{plot_sizes} (list or scalar): Sizes matching \code{plots}.
#'     \item \code{ora_results} (optional, list): Raw ORA output used by
#'       \code{generate_section_content()} to append a filtered table and a
#'       download button under the plots.
#'   }
#' @param level_headers_info Kept for backward compatibility; ignored when
#'   \code{sections} is provided.
#' @param report_info Named list with auxiliary report fields (e.g., databases,
#'   parameters, background genes) used by downstream helpers and the header.
#' @param output_file_path Character path where the HTML report is written.
#'
#' @return This function writes an HTML file to \code{output_file_path} and
#'   returns \code{NULL} invisibly.
#'
#' @details
#' The function initializes the HTML with \code{header_section} and a TOC
#' placeholder, then iterates over \code{sections}. For each section it:
#' \enumerate{
#'   \item Renders an \code{<h2>} header (with an \code{<hr>} after the first).
#'   \item Filters out invalid/empty plot objects and renders remaining plots
#'         via \code{process_plots()}.
#'   \item If no plot is available, inserts a centered note:
#'         "No enrichment plot available for this section."
#'   \item If \code{ora_results} is present, calls
#'         \code{generate_section_content()} to append a filtered ORA table and
#'         a base64 download link for the full results. If no ORA is present
#'         and no plots were rendered, a centered note is shown:
#'         "No enrichment results for this section."
#' }
#' The TOC is updated with one entry per section. Finally,
#' \code{generate_and_write_html()} writes the complete HTML document.
#'
build_run_ora_report <- function(
    header_section,
    sections,
    report_info,
    output_file_path) {
    html_content <- paste(header_section, "<!--TOC-->", sep = "\n")

    toc <- create_toc()
    styles <- define_html_styles()
    section_header_style <- styles$section_header_style
    toc_style <- styles$toc_style

    # Helper to decide if a "plot" is really renderable
    plot_is_valid <- function(p) {
        if (is.null(p)) {
            return(FALSE)
        }
        if (is.character(p)) {
            return(FALSE)
        } # avoid old sentinels/paths
        if (inherits(p, "ggplot")) {
            return(TRUE)
        } # ggplot object OK
        if (is.list(p) && !is.null(p$path)) {
            return(is.character(p$path) && file.exists(p$path))
        }
        FALSE
    }

    # Progress bar over valid plots only
    total_plots <- sum(vapply(sections, function(s) {
        ps <- s$plots
        if (is.null(ps)) {
            return(0L)
        }
        if (!is.list(ps) || inherits(ps, "ggplot")) ps <- list(ps)
        sum(vapply(ps, plot_is_valid, logical(1)))
    }, integer(1)))
    pb <- create_progress_bar(seq_len(max(total_plots, 1)))

    for (s_idx in seq_along(sections)) {
        sec <- sections[[s_idx]]
        header_name <- if (
            !is.null(sec$header_name) && nzchar(sec$header_name)
            ) {
            sec$header_name
        } else {
            paste0("Section ", s_idx)
        }

        # Header + optional HR
        section_header <- sprintf(
            "<h2 style='%s' id='section%d'>%s</h2>",
            section_header_style, s_idx, header_name
        )
        html_content <- paste(
            html_content,
            if (s_idx > 1) "<hr>" else "",
            section_header,
            sep = "\n"
        )

        # TOC entry
        toc_entry <- sprintf(
            "<li style='%s'><a href='#section%d'>%s</a></li>",
            toc_style, s_idx, header_name
        )
        toc <- paste(toc, toc_entry, sep = "\n")

        # Normalize and filter plots
        plots <- sec$plots
        if (!is.list(plots) || inherits(plots, "ggplot")) plots <- list(plots)

        sizes <- sec$plot_sizes
        if (!is.list(sizes)) {
            sizes <- as.list(rep_len(
                if (is.null(sizes)) 999 else sizes,
                length(plots)
            ))
        }

        valid <- vapply(plots, plot_is_valid, logical(1))
        plots <- plots[valid]
        sizes <- sizes[valid]

        if (length(plots) == 0L) {
            # No valid plot -> tidy note instead of broken image
            html_content <- paste(
                html_content,
                "<p style='font-size:14px;color:#666;margin:0.5em 0 1em 0;
             text-align:center;'>",
                "No enrichment plot available for this section.",
                "</p>",
                sep = "\n"
            )
        } else {
            for (i in seq_along(plots)) {
                res <- process_plots(
                    plots_element = plots[[i]],
                    plots_size    = sizes[[i]],
                    html_content  = html_content,
                    toc           = toc,
                    header_index  = s_idx
                )
                html_content <- res$html_content
                toc <- res$toc
                pb$tick()
            }
        }

        # Append filtered table + download (your existing function)
        if (!is.null(sec$ora_results)) {
            add <- generate_section_content(
                section_info = list(
                    header_name = header_name,
                    ora_results = sec$ora_results
                ),
                index = s_idx,
                toc = toc,
                html_content = html_content,
                section_header_style = section_header_style,
                toc_style = toc_style
            )
            html_content <- add$html_content
        } else if (length(plots) == 0L) {
            # No ORA results AND no plots -> show a single concise note
            html_content <- paste(
                html_content,
                "<p style='font-size:14px;color:#666;margin:0.25em 0 1.25em 0;
             text-align:center;'>",
                "No enrichment results for this section.",
                "</p>",
                sep = "\n"
            )
        }
    }

    generate_and_write_html(
        toc = toc,
        html_content = html_content,
        report_info = report_info,
        output_file_path = output_file_path
    )
}


#' Map gene symbols across species
#'
#' @noRd
#'
#' @description
#' This function maps gene symbols from one species to another using either
#' the g:Profiler API (`gprofiler`) or the `orthogene` package (`orthogene`).
#' It preserves the input order and length, returning a vector of mapped gene
#' symbols aligned to the input. If `method = "none"`, the input is returned
#' unchanged.
#'
#' @param genes Character vector of gene symbols.
#' @param mapping_cfg Named list controlling cross-species mapping:
#'   \describe{
#'     \item{method}{One of `"none"` (default), `"gprofiler"`, `"orthogene"`.}
#'     \item{from_species}{Source organism code (e.g., `"cgriseus"`,
#'       `"mmusculus"`).}
#'     \item{to_species}{Target organism code (e.g., `"hsapiens"`,
#'       `"mmusculus"`).}
#'   }
#'
#' @return Character vector of mapped symbols, same length and order as
#'   \code{genes}.
#'
map_gene_symbols <- function(
    genes,
    mapping_cfg = list(method = "none")) {
    stopifnot(is.character(genes))
    if (is.null(mapping_cfg$method)) mapping_cfg$method <- "none"
    method <- tolower(mapping_cfg$method)

    if (identical(method, "none")) {
        return(genes)
    }

    from <- mapping_cfg$from_species
    to <- mapping_cfg$to_species
    if (is.null(from) || is.null(to)) {
        stop("For method='", method,
            "', please provide from_species and to_species.",
            call. = FALSE
        )
    }

    has_suffix <- !is.na(genes) & grepl("_[0-9]+$", genes)
    if (any(has_suffix)) {
        warning("Some gene IDs end in '_<digits>'. No automatic stripping ",
            "performed.",
            call. = FALSE
        )
    }

    unique_genes <- unique(genes[!is.na(genes)])

    pick_first_per_input <- function(df, in_col, out_col) {
        df <- df[!is.na(df[[out_col]]), , drop = FALSE]
        df <- df[!duplicated(df[[in_col]]), , drop = FALSE]
        stats::setNames(df[[out_col]], df[[in_col]])
    }

    map_vec <- switch(method,
        "gprofiler" = {
            if (!requireNamespace("gprofiler2", quietly = TRUE)) {
                stop_call_false(
                    "'gprofiler2' not installed; install it or set",
                    "method='none'."
                    )
            }
            gp <- gprofiler2::gorth(
                unique_genes,
                source_organism = from,
                target_organism = to,
                mthreshold = Inf,
                filter_na = TRUE
            )
            if (!all(c("input", "ortholog_name") %in% names(gp))) {
                stop(
                    "gprofiler2::gorth returned unexpected columns: ",
                    paste(names(gp), collapse = ", ")
                )
            }
            pick_first_per_input(gp,
                in_col = "input",
                out_col = "ortholog_name"
            )
        },
        "orthogene" = {
            if (!requireNamespace("orthogene", quietly = TRUE)) {
                stop_call_false(
                    "'orthogene' not installed; install it or set",
                    "method='none'."
                    )
            }
            or <- orthogene::convert_orthologs(
                gene_df        = unique_genes,
                input_species  = from,
                output_species = to,
                method         = "gprofiler"
            )
            in_col <- intersect(
                names(or),
                c("input_gene", "input", "gene")
            )
            out_col <- intersect(
                names(or),
                c(
                    "ortholog_gene", "ortholog_symbol",
                    "output_gene", "ortholog_name"
                )
            )
            if (length(in_col) == 0 || length(out_col) == 0) {
                stop(
                    "orthogene::convert_orthologs returned unexpected",
                    "columns: ",
                    paste(names(or), collapse = ", ")
                )
            }
            pick_first_per_input(or,
                in_col = in_col[1],
                out_col = out_col[1]
            )
        },
        stop("Unknown mapping method: ", method, call. = FALSE)
    )

    id_map <- stats::setNames(unique_genes, unique_genes)
    id_map[names(map_vec)] <- map_vec
    idx <- match(genes, names(id_map))
    mapped <- genes
    mapped[!is.na(idx)] <- unname(id_map[idx[!is.na(idx)]])

    return(mapped)
}


# Level 2 internal functions ---------------------------------------------------


#' Validate clustering results table
#'
#' @noRd
#'
#' @description
#' Validates a unified clustering results tibble produced by
#' \code{construct_cluster_summary()}. Ensures required metadata columns
#' exist (\code{feature_nr}, \code{feature_name}, \code{gene}), verifies
#' that there are at least two time-effect cluster columns named like
#' \code{cluster_<cond1>} and \code{cluster_<cond2>}, checks that
#' \code{feature_nr} is integer-like and unique, and that
#' \code{feature_name} is non-missing and non-empty. Also validates the
#' optional category columns \code{cluster_cat2} and \code{cluster_cat3}
#' (when present) to be atomic vectors with values either \code{NA} or
#' strings of the form \code{"<cluster_<cond1>>_<cluster_<cond2>>"}.
#'
#' @param cluster_table A tibble containing one row per
#'   \code{feature_nr} with metadata and cluster assignments across the
#'   analysis categories. It includes:
#'   \itemize{
#'     \item \code{feature_nr} - Numeric feature identifier.
#'     \item \code{feature_name} - Preferred feature name from the source
#'       data, falling back to the numeric ID if none is available.
#'     \item \code{gene} - Preferred gene symbol from the annotation or
#'       cluster data.
#'     \item \code{cluster_<cond1>} / \code{cluster_<cond2>} - Cluster
#'       assignments for each time-effect condition.
#'     \item \code{cluster_cat2} - (Optional) Combined cluster label for
#'       category 2 hits in the form
#'       \code{"<cluster_<cond1>>_<cluster_<cond2>>"}; \code{NA} if the
#'       feature was not a category 2 hit.
#'     \item \code{cluster_cat3} - (Optional) Combined cluster label for
#'       category 3 hits in the form
#'       \code{"<cluster_<cond1>>_<cluster_<cond2>>"}; \code{NA} if the
#'       feature was not a category 3 hit.
#'   }
#'   For any category-specific cluster column, a value of \code{NA}
#'   indicates that the feature was not significant (not a hit) in that
#'   category.
#'
#' @return Invisibly returns \code{TRUE} on success. The function
#'   \emph{stops with an informative error} if any validation fails.
#'
check_cluster_table <- function(cluster_table) {
    if (!is.data.frame(cluster_table)) {
        stop("cluster_table must be a data.frame/tibble.", call. = FALSE)
    }

    req <- c("feature_nr", "feature_name", "gene")
    missing <- setdiff(req, names(cluster_table))
    if (length(missing) > 0) {
        stop(paste0(
            "Missing required columns: ",
            paste(missing, collapse = ", ")
        ), call. = FALSE)
    }

    clust_cols <- grep("^cluster_", names(cluster_table), value = TRUE)
    cond_cols <- setdiff(clust_cols, c("cluster_cat2", "cluster_cat3"))
    if (length(cond_cols) < 2) {
        stop("Expected at least two condition cluster columns named like ",
            "'cluster_<cond1>' and 'cluster_<cond2>'.",
            call. = FALSE
        )
    }

    fn <- cluster_table$feature_nr
    int_like <- is.integer(fn) ||
        (is.numeric(fn) && all(fn == as.integer(fn), na.rm = TRUE))
    if (!int_like) {
        stop("'feature_nr' must be integer-like.", call. = FALSE)
    }

    if (any(duplicated(cluster_table$feature_nr))) {
        dups <- unique(cluster_table$feature_nr[
            duplicated(cluster_table$feature_nr)
        ])
        stop(
            paste0(
                "Duplicate feature_nr values: ",
                paste(head(dups, 10), collapse = ", "),
                if (length(dups) > 10) " ..." else ""
            ),
            call. = FALSE
        )
    }

    fnm <- cluster_table$feature_name
    if (any(is.na(fnm) | fnm == "")) {
        stop("'feature_name' must be non-missing and non-empty.", call. = FALSE)
    }

    for (cc in cond_cols) {
        v <- cluster_table[[cc]]
        if (!is.atomic(v)) {
            stop(paste0("Condition cluster column '", cc, "' must be atomic."),
                call. = FALSE
            )
        }
    }

    for (cc in c("cluster_cat2", "cluster_cat3")) {
        if (cc %in% names(cluster_table)) {
            v <- cluster_table[[cc]]
            if (!is.atomic(v)) {
                stop(paste0("'", cc, "' must be atomic."), call. = FALSE)
            }

            if (cc == "cluster_cat2") {
                bad <- !is.na(v) & !grepl("^[^_]+_[^_]+$", v)
                if (any(bad)) {
                    stop(paste0("'", cc, "' must be of the form 'a_b' or NA."),
                        call. = FALSE
                    )
                }
            }

            if (cc == "cluster_cat3") {
                bad <- !is.na(v) & !is.character(v)
                if (any(bad)) {
                    stop(paste0("'", cc, "' must be character or NA."),
                        call. = FALSE
                    )
                }
            }
        }
    }

    invisible(TRUE)
}


#' Check Valid Gene IDs (permissive, no mutation)
#'
#' @noRd
#'
#' @description
#' Minimal validation that \code{genes} is a character vector and, if
#' \code{max_index_overall} is provided, that \code{length(genes)} is at
#' least that value. No content validation is performed to avoid false
#' negatives. This function does not modify or return \code{genes}.
#'
#' @param genes Character vector of gene IDs/symbols.
#' @param max_index_overall Optional single integer: required minimum length
#'   of \code{genes} (e.g., highest feature index expected).
#'
#' @return Invisibly returns \code{TRUE}. Stops only on type/length errors.
#'
check_genes <- function(
    genes,
    max_index_overall = NA_integer_) {
    if (!is.character(genes)) {
        stop("'genes' must be a character vector.", call. = FALSE)
    }

    if (!is.na(max_index_overall)) {
        if (!is.numeric(max_index_overall) ||
            length(max_index_overall) != 1L ||
            is.na(max_index_overall)) {
            stop("'max_index_overall' must be a single numeric value.",
                call. = FALSE
            )
        }
        if (length(genes) < as.integer(max_index_overall)) {
            stop(sprintf(
                "`genes` must have length >= %d (got %d).",
                as.integer(max_index_overall), length(genes)
            ), call. = FALSE)
        }
    }

    invisible(TRUE)
}


#' Check Valid Databases Dataframe
#'
#' @noRd
#'
#' @description
#' This function checks if the dataframe has exactly three columns
#' named DB, Geneset, and Gene, and all columns must be of type
#' character.
#'
#' @param databases A dataframe to check.
#'
#' @return None. This function stops execution and provides an
#'         error message if the dataframe is not valid.
#'
check_databases <- function(databases) {
    if (!is.data.frame(databases)) {
        stop("The input must be a dataframe.", call. = FALSE)
    }

    if (ncol(databases) != 3) {
        stop("The dataframe must have exactly three columns.", call. = FALSE)
    }

    expected_colnames <- c("DB", "Geneset", "Gene")
    if (!all(colnames(databases) == expected_colnames)) {
        stop("The dataframe must have columns named DB, Geneset, and Gene.",
            call. = FALSE
        )
    }

    if (!all(vapply(databases, is.character, logical(1)))) {
        stop_call_false(
            "All columns in the dataframe must be of type character."
            )
    }
}


#' Check Params List for Required Conditions
#'
#' @noRd
#'
#' @description
#' This function checks if a given list `params` contains only the allowed
#' named elements. The elements do not have to be present, but if they are,
#' they must be
#' named exactly as specified and must contain the correct data types: float,
#' character, int,
#' int, and float. If any condition is not met, the function stops the script
#' and produces an
#' informative error message. `params` can also be `NA`.
#'
#' @param params A list to be checked for the required conditions, or `NA`.
#'
#' @return This function does not return a value. It stops with an error message
#'         if the conditions are not met.
#'
#' @importFrom stats p.adjust
#'
check_params <- function(params) {
    required_params <- list(
        pvalueCutoff = "numeric",
        pAdjustMethod = "character",
        minGSSize = "numeric",
        maxGSSize = "numeric",
        qvalueCutoff = "numeric"
    )

    # Check if params is a list
    if (!is.list(params)) {
        stop("The input must be a list.", call. = FALSE)
    }

    # Check for extra elements
    extra_params <- setdiff(names(params), names(required_params))
    if (length(extra_params) > 0) {
        stop_call_false(
            paste(
                "The list contains extra elements besides the allowed elements",
                "pvalueCutoff, pAdjustMethod, minGSSize, maxGSSize and",
                "qvalueCutoff:",
                paste(extra_params, collapse = ", ")
            )
        )
    }

    valid_adj_p_value_methods <- stats::p.adjust.methods

    # Check for required elements and their data types
    for (param in names(params)) {
        if (param %in% names(required_params)) {
            actual_value <- params[[param]]
            expected_type <- required_params[[param]]
            actual_type <- class(params[[param]])

            if (is.null(actual_value)) {
                next
            }

            if (expected_type == "integer" && !is.integer(params[[param]])) {
                # Check if the value can be coerced to integer
                if (is.numeric(params[[param]]) &&
                    all(params[[param]] == as.integer(params[[param]]))) {
                    actual_type <- "integer"
                }
            }
            if (expected_type != actual_type) {
                stop("The element '", param, "' must be of type ",
                    expected_type, ".",
                    call. = FALSE
                )
            }
        } else {
            stop(
                "Unexpected element '",
                param,
                "' in the list.",
                call. = FALSE
            )
        }

        if (param == "pAdjustMethod") {
            if (!(actual_value %in% valid_adj_p_value_methods)) {
                stop(
                    paste("pAdjustMethod must be one of",
                        valid_adj_p_value_methods,
                        collapse = ", "
                    ),
                    call. = FALSE
                )
            }
        }
    }
}


#' Check the user-supplied mapping configuration
#'
#' @noRd
#'
#' @description
#'
#' This routine validates a `mapping_cfg` list **without altering it** and
#' raises a clear, actionable `stop()` message on the first problem detected.
#' If everything is acceptable it returns `invisible(TRUE)`.
#'
#' @param mapping_cfg list. Must contain **exactly** the elements
#'   `method`, `from_species`, and `to_species`.
#' @param verbose Boolean flag controlling the display of messages.
#'
#' @details
#' * **method** - character scalar, one of `"none"`, `"gprofiler"`,
#'   or `"orthogene"` (case-insensitive).
#' * **from_species**, **to_species** - character scalars giving the species
#'   identifiers expected by the chosen tool.  They are only checked (not
#'   auto-corrected).
#'
#' The function deliberately performs *no* coercion or guessing: any deviation
#' from the expected input is considered an error so that the user must supply
#' an unambiguous, reproducible configuration.
#'
check_mapping_cfg <- function(
    mapping_cfg,
    verbose) {
    # structural checks
    if (!is.list(mapping_cfg) || inherits(mapping_cfg, "data.frame")) {
        stop(
            "`mapping_cfg` must be a *list*, e.g.\n",
            "  list(method = 'gprofiler', from_species = 'cgchok1gshd', ",
            "to_species = 'hsapiens')."
        )
    }

    expected <- c("method", "from_species", "to_species")
    missing <- setdiff(expected, names(mapping_cfg))
    if (length(missing)) {
        stop(
            "`mapping_cfg` is missing field(s): ",
            paste(missing, collapse = ", "), "."
        )
    }

    extra <- setdiff(names(mapping_cfg), expected)
    if (length(extra)) {
        stop(
            "`mapping_cfg` contains unknown field(s): ",
            paste(extra, collapse = ", "), ".\n",
            "Only the fields ", paste(expected, collapse = ", "),
            " are allowed."
        )
    }

    # method checks
    method <- tolower(mapping_cfg$method[[1]])
    if (!is.character(mapping_cfg$method) || length(mapping_cfg$method) != 1) {
        stop("`method` must be a single character string.")
    }
    ok_methods <- c("none", "gprofiler", "orthogene")
    if (!method %in% ok_methods) {
        stop(
            "`method` must be one of ",
            paste(shQuote(ok_methods), collapse = ", "),
            ". You supplied: ", shQuote(mapping_cfg$method), "."
        )
    }

    # species checks
    need_species <- method != "none"
    if (need_species) {
        for (sp in c("from_species", "to_species")) {
            val <- mapping_cfg[[sp]]
            if (!is.character(val) || length(val) != 1 || !nzchar(val)) {
                stop("`", sp, "` must be a non-empty character scalar.")
            }
        }

        if (identical(mapping_cfg$from_species, mapping_cfg$to_species)) {
            stop(
                "`from_species` and `to_species` must not be the same: ",
                shQuote(mapping_cfg$from_species), "."
            )
        }
    }

    # tool-specific checks -
    if (method == "gprofiler") {
        if (!requireNamespace("gprofiler2", quietly = TRUE)) {
            stop_call_false(
                "`gprofiler2` is not installed. Install with:
        install.packages('gprofiler2')"
            )
        }

        # Skip organism code checking entirely
        if (isTRUE(verbose)) {
            message(
                "Note: `from_species` and `to_species` are not checked for",
                "validity.\n",
                "Make sure they match g:Profiler's expected organism codes.",
                "See:\n",
                "https://biit.cs.ut.ee/gprofiler/page/organism-list"
            )
        }
    }

    if (method == "orthogene") {
        if (!requireNamespace("orthogene", quietly = TRUE)) {
            stop_call_false(
                "`orthogene` is not installed. Install with:\n",
                "  BiocManager::install('orthogene')"
            )
        }
        # orthogene's own mapper returns NA for unknown species
        bad_from <- is.na(orthogene::map_species(
            mapping_cfg$from_species,
            quiet = TRUE,
            verbose = FALSE
        ))
        if (bad_from) {
            stop_call_false(
                "`orthogene` cannot recognise `from_species` = '",
                mapping_cfg$from_species, "'.\n",
                "Use a common/scientific name (e.g. 'Chinese hamster') ",
                "or a valid NCBI tax-ID (10029)."
            )
        }
        bad_to <- is.na(orthogene::map_species(
            mapping_cfg$to_species,
            quiet = TRUE,
            verbose = FALSE
        ))
        if (bad_to) {
            stop_call_false(
                "`orthogene` cannot recognise `to_species` = '",
                mapping_cfg$to_species, "'.\n",
                "Use a common/scientific name (e.g. 'human') ",
                "or a valid NCBI tax-ID (9606)."
            )
        }
    }

    invisible(TRUE)
}


#' Validate the `enrichGO_cfg` Configuration
#'
#' @noRd
#'
#' @description
#' This function validates the user-supplied \code{enrichGO_cfg} argument to
#' ensure it is properly specified for running GO enrichment through
#' \code{\link[clusterProfiler]{enrichGO}}. It checks that the input is a named
#' list with appropriate entries and that each entry provides the necessary
#' fields in the correct format.
#'
#' Each list element must be named according to the ontology type
#' (e.g., \code{"GO_BP"}, \code{"GO_MF"}, \code{"GO_CC"}) and must contain:
#' \itemize{
#'   \item \code{OrgDb}: The organism annotation database, e.g.,
#'    \code{org.Mm.eg.db}.
#'   \item \code{keyType}: The gene identifier type as a string, e.g.,
#'    \code{"SYMBOL"}.
#'   \item \code{ontology}: One of \code{"BP"}, \code{"MF"}, or \code{"CC"}.
#' }
#'
#' If any required element is missing or improperly specified, the function
#' raises an informative error via \code{stop_call_false()}.
#'
#' @param enrichGO_cfg A named list specifying the configuration for
#' running GO enrichment with Bioconductor's
#' \code{\link[clusterProfiler]{enrichGO}}.
#' This is only needed when you want to perform GO Biological Process (BP),
#' Molecular Function (MF), or Cellular Component (CC) enrichment using
#' Bioconductor's organism databases (e.g., \code{org.Mm.eg.db} for mouse).
#'
#' The list must be named according to the GO ontology, e.g., \code{"GO_BP"},
#' \code{"GO_MF"}, \code{"GO_CC"}. Each entry must provide:
#' \itemize{
#'   \item \code{OrgDb}: The organism database, e.g., \code{org.Mm.eg.db}.
#'   \item \code{keyType}: The gene identifier type, e.g., \code{"SYMBOL"}.
#'   \item \code{ontology}: One of \code{"BP"}, \code{"MF"}, or \code{"CC"}.
#' }
#'
#' If \code{enrichGO_cfg} is \code{NULL} (default), no Bioconductor-based GO
#' enrichment is performed. All enrichment runs through
#' \code{\link[clusterProfiler]{enricher}} with the provided TERM2GENE mappings.
#'
#' @return Invisibly returns \code{NULL} if all checks pass; otherwise,
#' stops with an informative error message.
#'
check_enrichGO_cfg <- function(enrichGO_cfg) {
    if (is.null(enrichGO_cfg)) {
        return(invisible(NULL))
    }

    if (!is.list(enrichGO_cfg)) {
        stop_call_false("`enrichGO_cfg` must be a named list.")
    }

    cfg_names <- names(enrichGO_cfg)

    if (is.null(cfg_names) || any(cfg_names == "")) {
        stop_call_false(
            "`enrichGO_cfg` must be a named list with non-empty names."
            )
    }

    if (!all(startsWith(cfg_names, "GO_"))) {
        stop_call_false(
            "All names in `enrichGO_cfg` must start with 'GO_', e.g., 'GO_BP'."
        )
    }

    required_fields <- c("OrgDb", "keyType", "ontology")

    for (cfg_name in cfg_names) {
        cfg <- enrichGO_cfg[[cfg_name]]

        if (!is.list(cfg)) {
            stop_call_false("Each entry in `enrichGO_cfg` must be a list.")
        }

        missing_fields <- setdiff(required_fields, names(cfg))
        if (length(missing_fields) > 0) {
            stop_call_false(
                paste0(
                    "Entry '", cfg_name, "' is missing required fields: ",
                    paste(missing_fields, collapse = ", ")
                )
            )
        }

        if (!inherits(cfg$OrgDb, "OrgDb")) {
            stop_call_false(
                paste0(
                    "In entry '",
                    cfg_name,
                    "': 'OrgDb' must be a valid Bioconductor OrgDb object",
                    "(e.g., org.Mm.eg.db)."
                )
            )
        }

        if (!is.character(cfg$keyType) || length(cfg$keyType) != 1) {
            stop_call_false(
                paste0(
                    "In entry '",
                    cfg_name,
                    "': 'keyType' must be a single string (e.g., 'SYMBOL')."
                )
            )
        }

        valid_ontologies <- c("BP", "MF", "CC")
        if (!cfg$ontology %in% valid_ontologies) {
            stop_call_false(
                paste0(
                    "In entry '",
                    cfg_name,
                    "': 'ontology' must be one of 'BP', 'MF', or 'CC'."
                )
            )
        }
    }

    invisible(NULL)
}


#' Perform ORA for clustered genes and plot results
#'
#' @noRd
#'
#' @description
#' Runs over-representation analysis (ORA) with \pkg{clusterProfiler}
#' for each distinct value in a \code{cluster} column, using the genes
#' in \code{gene} as foreground. Gene-set collections are taken from
#' \code{databases} via \code{dbs_to_term2genes()}. For GO collections
#' provided in \code{enrichGO_cfg}, \code{enrichGO()} (optionally
#' followed by \code{simplify()}) is used; all other collections are
#' analyzed with \code{enricher()} and their \code{TERM2GENE} maps.
#' Invalid/empty gene IDs are removed and duplicates are deduplicated
#' before enrichment. If any cluster yields results, a dot plot is
#' generated; otherwise \code{NA} is returned.
#'
#' @param clustered_genes A data.frame/tibble with two columns:
#'   \code{gene} (character gene identifiers compatible with the chosen
#'   collections) and \code{cluster} (scalar label; may be numeric or
#'   character such as \code{"cond1_higher"}). One row per
#'   feature/gene-instance; foregrounds are built by grouping on
#'   \code{cluster}.
#' @param databases A data.frame or list describing the gene-set
#'   collections to analyze. It is converted to \code{TERM2GENE} maps by
#'   \code{dbs_to_term2genes()}.
#' @param params A list of clusterProfiler settings (e.g.,
#'   \code{pvalueCutoff}, \code{pAdjustMethod}, \code{minGSSize},
#'   \code{maxGSSize}, \code{qvalueCutoff}). Passed through to the
#'   enrichment functions.
#' @param enrichGO_cfg Optional named list configuring GO runs for those
#'   collections that are GO. Each entry should provide:
#'   \code{OrgDb}, \code{keyType}, and \code{ontology} (one of
#'   \code{"BP"}, \code{"MF"}, \code{"CC"}). When present for a given
#'   collection name, \code{enrichGO()} is used instead of
#'   \code{enricher()}.
#' @param plot_title Optional character string used as a title/label for
#'   the generated dot plot(s).
#' @param universe Optional character vector of background genes used as
#'   the universe for ORA (passed to \code{enrichGO()}/\code{enricher()}).
#' @param verbose Boolean flag controlling the display of messages.
#'
#' @return A list with:
#'   \itemize{
#'     \item \code{dotplot}: ggplot object (or list of ggplots) showing
#'       ORA results across clusters/collections.
#'     \item \code{dotplot_nrows}: suggested plot height/rows.
#'     \item \code{ora_results}: nested list
#'       \code{cluster_<label> -> collection -> data.frame} with raw
#'       enrichment tables (or \code{NA} when empty).
#'   }
#'   If no enrichment is found for any cluster, returns \code{NA}.
#'
run_ora_level <- function(
    clustered_genes,
    databases,
    params = NA,
    enrichGO_cfg = NULL,
    plot_title = "",
    universe = NULL,
    verbose = TRUE) {
    params <- set_default_params(params)
    gene_set_collections <- dbs_to_term2genes(databases)
    unique_clusters <- sort(unique(clustered_genes$cluster))
    ora_results <- list()

    for (cluster in unique_clusters) {
        cluster_label <- paste0("cluster_", cluster)
        if (isTRUE(verbose)) message(paste("\nCluster:", cluster_label))
        ora_results[[cluster_label]] <- list()

        fg <- as.character(clustered_genes$gene[
            clustered_genes$cluster == cluster
        ])
        fg <- unique(fg[!is.na(fg) & fg != ""])
        if (length(fg) == 0L) next

        for (gene_set_name in names(gene_set_collections)) {
            if (isTRUE(verbose)) message(paste("Database:", gene_set_name))
            if (
                !is.null(enrichGO_cfg) && gene_set_name %in% names(enrichGO_cfg)
                ) {
                cfg <- enrichGO_cfg[[gene_set_name]]
                ora_result <- clusterProfiler::enrichGO(
                    gene = fg,
                    OrgDb = cfg$OrgDb,
                    keyType = cfg$keyType,
                    ont = cfg$ontology,
                    pvalueCutoff = params$pvalueCutoff,
                    pAdjustMethod = params$pAdjustMethod,
                    universe = universe,
                    minGSSize = params$minGSSize,
                    maxGSSize = params$maxGSSize,
                    qvalueCutoff = params$qvalueCutoff
                )
                # Only simplify if there are results
                if (
                    !is.null(ora_result) && nrow(as.data.frame(ora_result)) > 0
                    ) {
                    ora_result <- clusterProfiler::simplify(ora_result)
                }
            } else {
                gene_set_map <- gene_set_collections[[gene_set_name]]
                check_gene_overlap(
                    fg = fg,
                    gene_set_map = gene_set_map,
                    verbose = verbose
                )
                ora_result <- clusterProfiler::enricher(
                    gene = fg,
                    pvalueCutoff = params$pvalueCutoff,
                    pAdjustMethod = params$pAdjustMethod,
                    universe = universe,
                    minGSSize = params$minGSSize,
                    maxGSSize = params$maxGSSize,
                    qvalueCutoff = params$qvalueCutoff,
                    gson = NULL,
                    TERM2GENE = gene_set_map,
                    TERM2NAME = NA
                )
            }

            df <- as.data.frame(ora_result)
            if (nrow(df) == 0) df <- NA
            ora_results[[cluster_label]][[gene_set_name]] <- df
        }
    }

    ora_results <- add_odds_ratios_to_ora(ora_results)

    any_result <- any(vapply(ora_results, function(cluster_entry) {
        any(vapply(cluster_entry, function(res) {
            is.data.frame(res) && nrow(res) > 0L
        }, logical(1)))
    }, logical(1)))

    if (any_result) {
        result <- make_enrich_dotplot(
            ora_results = ora_results,
            title = plot_title
            )
    } else {
        if (isTRUE(verbose)) message("No cluster led to an enrichment result!")
        return(NA)
    }

    list(
        dotplot = result[["dotplot"]],
        dotplot_nrows = result[["dotplot_height"]],
        ora_results = ora_results
    )
}


#' Generate Section Content
#'
#' @noRd
#'
#' @description
#' Generates the HTML content for a section, including headers and enrichment
#' results.
#'
#' @param section_info A list containing the information for the section.
#' @param index The index of the current section.
#' @param toc The current state of the Table of Contents.
#' @param html_content The current state of the HTML content.
#' @param section_header_style The CSS style for the section headers.
#' @param toc_style The CSS style for the TOC entries.
#'
#' @return A list with updated HTML content and TOC.
#'
generate_section_content <- function(
    section_info,
    index,
    toc,
    html_content,
    section_header_style,
    toc_style) {
    ora_results <- section_info$ora_results

    # Filtered results (only Count >= 2 and adjusted p < 0.05)
    top_df <- prepare_plot_data(ora_results)

    if (nrow(top_df) == 0) {
        no_results_message <- paste0(
            "<p style='font-size: 40px; color: #FF0000;'>",
            "No gene set showed statistically significant ",
            "overrepresentation in any cluster.",
            "</p>"
        )

        html_content <- paste(
            html_content,
            no_results_message,
            sep = "\n"
        )

        return(list(html_content = html_content))
    }

    # Full unfiltered results
    full_df <- flatten_ora_results(ora_results)

    # Header for filtered ORA results shown in HTML
    ora_results_header <- paste0(
        "<h3 style='font-size: 30px; font-weight: bold; color: #333;'>",
        "Filtered Overrepresentation Analysis (ORA) Results</h3>"
    )

    # Create HTML table for filtered results
    html_table <- "<table style='width:100%;border-collapse:collapse;'>"
    html_table <- paste0(html_table, "<thead><tr>")
    for (header in colnames(top_df)) {
        html_table <- paste0(html_table, "<th>", header, "</th>")
    }
    html_table <- paste0(html_table, "</tr></thead><tbody>")

    for (i in seq_len(nrow(top_df))) {
        html_table <- paste0(html_table, "<tr>")
        for (j in seq_len(ncol(top_df))) {
            html_table <- paste0(html_table, "<td>", top_df[i, j], "</td>")
        }
        html_table <- paste0(html_table, "</tr>")
    }

    ora_results_html <- paste0(html_table, "</tbody></table>")

    # Header for Excel download of full ORA results
    full_ora_header <- paste0(
        "<h3 style='font-size: 30px; font-weight: bold; color: #333;'>",
        "Full ORA Results (including terms supported by < 2 genes)</h3>"
    )

    # Generate base64 download for full ORA result
    base64_df <- sprintf(
        '<a href="%s" download="full_ora_results.xlsx">
  <button>Download full_ora_results.xlsx</button></a>',
        encode_df_to_base64(
            full_df,
            "run_ora"
        )
    )

    html_content <- paste(
        html_content,
        ora_results_header,
        ora_results_html,
        full_ora_header,
        base64_df,
        sep = "\n"
    )

    list(html_content = html_content)
}


#' Determine level columns in clustering summary table
#'
#' @noRd
#'
#' @description
#' Returns the names of columns in a \code{clustering_results} data frame
#' that represent enrichment levels. The order is:
#' (1) time-effect condition columns (\code{cluster_*} except
#' \code{cluster_cat2} and \code{cluster_cat3}), followed by
#' (2) \code{cluster_cat2} if present, and
#' (3) \code{cluster_cat3} if present.
#'
#' @param cr A tibble/data frame like that produced by
#'   \code{construct_cluster_summary()}, containing \code{cluster_*}
#'   columns.
#'
#' @return A character vector of level column names in the preferred
#'   processing order.
#'
#' @details
#' Columns are identified via \code{grep("^cluster_", names(cr))}. The
#' time-effect condition columns retain their original ordering from
#' \code{cr}; category columns are appended if present.
#'
level_columns <- function(cr) {
    allc <- grep("^cluster_", names(cr), value = TRUE)
    cond <- setdiff(allc, c("cluster_cat2", "cluster_cat3"))
    c(
        cond,
        intersect("cluster_cat2", allc),
        intersect("cluster_cat3", allc)
    )
}


#' Derive display names for enrichment levels
#'
#' @noRd
#'
#' @description
#' Generates human-readable labels for the enrichment levels present in a
#' \code{clustering_results} table. Time-effect condition columns
#' (\code{cluster_*} except \code{cluster_cat2}/\code{cluster_cat3}) are
#' mapped to \code{"time_effect_condition_<cond>"}, while category columns
#' (if present) are named \code{"avrg_diff_conditions"} (cat2) and
#' \code{"interaction_condition_time"} (cat3).
#'
#' @param cr A tibble/data frame like that produced by
#'   \code{construct_cluster_summary()} containing \code{cluster_*}
#'   columns.
#'
#' @return A character vector of display names ordered to match the level
#'   columns returned by \code{level_columns(cr)}.
#'
#' @details
#' Condition names are derived by stripping the \code{"cluster_"} prefix
#' from the corresponding column names. Category labels are included only
#' if the respective columns exist in \code{cr}.
#'
level_display_names <- function(cr) {
    allc <- grep("^cluster_", names(cr), value = TRUE)
    cond <- setdiff(allc, c("cluster_cat2", "cluster_cat3"))
    c(
        paste0("time_effect_condition_", sub("^cluster_", "", cond)),
        if ("cluster_cat2" %in% allc) "avrg_diff_conditions",
        if ("cluster_cat3" %in% allc) "interaction_condition_time"
    )
}


#' Build foreground gene sets for a level
#'
#' @noRd
#'
#' @description
#' Extracts non-missing genes and their cluster labels from a unified
#' \code{clustering_results} table for a specific level column, producing
#' a two-column data frame suitable for ORA (\code{gene}, \code{cluster}).
#' Rows with missing/empty genes or missing level labels are dropped, and
#' duplicate gene-cluster pairs are removed.
#'
#' @param clustering_results A tibble/data frame produced by
#'   \code{construct_cluster_summary()} that includes a \code{gene} column
#'   and one or more level columns named \code{cluster_*}.
#' @param level_col Character scalar naming the level column in
#'   \code{clustering_results} to use (must exist in \code{names()}).
#'
#' @return A data frame with columns:
#'   \itemize{
#'     \item \code{gene} (character): foreground gene identifiers.
#'     \item \code{cluster} (character): labels from \code{level_col}.
#'   }
#'   May be empty (0 rows) if no valid gene-cluster pairs are found.
#'
#' @details
#' Rows are kept when \code{!is.na(level_col) & !is.na(gene) & gene != ""}.
#' Both columns are coerced to \code{character}, and duplicates are removed
#' via \code{unique()}.
#'
make_clustered_genes <- function(
    clustering_results,
    level_col) {
    stopifnot(level_col %in% names(clustering_results))
    x <- clustering_results
    keep <- !is.na(x[[level_col]]) & !is.na(x$gene) & x$gene != ""
    if (!any(keep)) {
        return(data.frame(
            gene = character(), cluster = character(),
            stringsAsFactors = FALSE
        ))
    }
    data.frame(
        gene = as.character(x$gene[keep]),
        cluster = as.character(x[[level_col]][keep]),
        stringsAsFactors = FALSE
    ) |>
        unique()
}


# Level 3 internal functions ---------------------------------------------------


#' Set Default Parameters
#'
#' @noRd
#'
#' @description
#' This function checks if the provided `params` list is `NA` or missing any
#' elements. If `params` is `NA`, it assigns a list of default parameters.
#' If any element is missing from `params`, it adds the missing element with
#' its respective default value.
#'
#' @param params A list of parameters to be checked and updated with default
#'               values if necessary.
#'
#' @return A list of parameters with all required elements, either from the
#'         input `params` or with added default values for any missing elements.
#'
set_default_params <- function(params) {
    default_params <- list(
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        minGSSize = 10,
        maxGSSize = 500,
        qvalueCutoff = 0.2
    )

    if (any(is.na(params))) {
        params <- default_params
    } else {
        # Check for missing elements and add default values if not present
        missing_params <- setdiff(names(default_params), names(params))
        for (param in missing_params) {
            params[[param]] <- default_params[[param]]
        }
    }

    return(params)
}


#' Convert DB/Geneset/Gene table to TERM2GENE list
#'
#' @noRd
#'
#' @description
#' Takes a data.frame with columns DB, Geneset, Gene and returns a named list,
#' one element per DB. Each element is a data.frame with two columns named
#' `term` (the Geneset) and `gene` (the Gene), suitable for
#' clusterProfiler::enricher(TERM2GENE = ...).
#'
#' @param databases data.frame with columns DB, Geneset, Gene
#' @return named list of data.frames (columns: term, gene), one per DB
#'
dbs_to_term2genes <- function(databases) {
    db_split <- split(databases, databases$DB)

    # Transform into long format
    all_term2genes <- lapply(db_split, function(db_df) {
        df_with_renamed_columns <- db_df[, c("Geneset", "Gene")]
        # Remove duplicate terms
        df_with_renamed_columns <- unique(df_with_renamed_columns)
        colnames(df_with_renamed_columns) <- c("term", "gene")
        return(df_with_renamed_columns)
    })

    names(all_term2genes) <- names(db_split)

    return(all_term2genes)
}


#' Add Odds Ratios to ORA Results
#'
#' @noRd
#'
#' @description
#' Computes and adds odds ratio values to each enrichment result data
#' frame within a nested ORA result structure. This function parses the
#' \code{GeneRatio} and \code{BgRatio} columns returned by
#' \code{clusterProfiler::enricher()}, calculates the odds ratio for each
#' term, and stores the result in a new column \code{odds_ratio}.
#'
#' The odds ratio is computed as:
#' \deqn{ (a / b) / (c / d) }
#' where \code{a/b} is the gene ratio (number of foreground genes in term
#' over total foreground genes), and \code{c/d} is the background ratio
#' (number of background genes in term over total background genes).
#'
#' @param ora_results A nested list of ORA results. Each top-level element
#'   represents a gene cluster. Each sub-list contains enrichment result
#'   data frames for specific gene sets. Data frames must include
#'   \code{GeneRatio} and \code{BgRatio} columns in "x/y" string format.
#'
#' @return A nested list with the same structure as \code{ora_results}, but
#'   with an additional numeric column \code{odds_ratio} in each result data
#'   frame, if applicable.
#'
#' @seealso \code{\link{make_enrich_dotplot}} for visualizing the results,
#'   and \code{\link{prepare_plot_data}} to extract filtered plot-ready data.
#'
#' @importFrom stats setNames
#'
add_odds_ratios_to_ora <- function(ora_results) {
    parse_ratio_vec <- function(x) {
        parts <- strsplit(x, "/", fixed = TRUE)
        numer_chr <- vapply(
            parts,
            function(p) {
                if (length(p) >= 2L) p[1L] else NA_character_
            },
            character(1)
        )
        
        denom_chr <- vapply(
            parts,
            function(p) {
                if (length(p) >= 2L) p[2L] else NA_character_
            },
            character(1)
        )
        numer <- as.numeric(numer_chr)
        denom <- as.numeric(denom_chr)
        numer / denom
    }

    for (cluster in names(ora_results)) {
        for (gene_set_name in names(ora_results[[cluster]])) {
            df <- ora_results[[cluster]][[gene_set_name]]

            if (is.data.frame(df) &&
                "GeneRatio" %in% names(df) &&
                "BgRatio" %in% names(df)) {
                gene_ratios <- parse_ratio_vec(df$GeneRatio)
                bg_ratios <- parse_ratio_vec(df$BgRatio)

                # Compute odds ratio (handles NA/Inf if bg_ratios == 0)
                df$odds_ratio <- gene_ratios / bg_ratios

                ora_results[[cluster]][[gene_set_name]] <- df
            }
        }
    }

    ora_results
}


#' Create Dotplot of ORA Results by Cluster and Gene Set
#'
#' @noRd
#'
#' @description
#' Generates a dotplot visualizing overrepresentation analysis (ORA)
#' results across multiple gene clusters and gene set databases. This
#' function consumes a nested ORA result structure and uses significance
#' and odds ratio information to represent enriched terms per cluster.
#'
#' The function relies on \code{prepare_plot_data()} to extract and filter
#' top enrichment results. Only terms with adjusted p-values below 0.05
#' and supported by at least two genes are shown. Term names are truncated
#' for readability if longer than 100 characters.
#'
#' @param ora_results A nested list of ORA results, where each top-level
#'   element corresponds to a gene cluster, and each sub-list contains
#'   enrichment result data frames for specific gene set databases. These
#'   data frames must include columns \code{p.adjust}, \code{Description},
#'   \code{odds_ratio}, and \code{Count}.
#' @param title A character string used as the plot title. Defaults to
#'   \code{"Title"}.
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{dotplot}{A \code{ggplot2} object representing the ORA dotplot.}
#'   \item{dotplot_height}{A numeric value for the optimal plot height, based
#'   on the number of unique terms.}
#' }
#'
#' @seealso \code{\link{prepare_plot_data}} for preparing the filtered input,
#'   and \code{\link{add_odds_ratios_to_ora}} to precompute odds ratios.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_blank ylab
#'   scale_color_gradient scale_size_area theme_bw theme element_blank
#'   element_text labs coord_cartesian scale_y_discrete guide_colorbar
#' @importFrom scales oob_squish
#' @importFrom grid unit
#' @importFrom rlang .data
#'
make_enrich_dotplot <- function(
    ora_results,
    title = "Title") {
    top_plot_data <- prepare_plot_data(ora_results)
    top_plot_data$cluster <- gsub(
        "^cluster_",
        "",
        top_plot_data$cluster
        )
    if (!is.data.frame(top_plot_data) || nrow(top_plot_data) == 0L) {
        # Optional: return a placeholder plot instead of NULL
        p_empty <- ggplot2::ggplot() +
            ggplot2::theme_void() +
            ggplot2::geom_text(
                ggplot2::aes(
                    x = 0,
                    y = 0,
                    label = "No enriched terms to display"
                    )
            )
        return(list(dotplot = p_empty, dotplot_height = 0.7))
    }

    height_per_label <- 0.1
    num_labels <- length(unique(top_plot_data$term))
    plot_height <- max(num_labels * height_per_label, 0.70)

    top_plot_data$term <- as.character(top_plot_data$term)
    top_plot_data$term <- ifelse(
        nchar(top_plot_data$term) > 100,
        paste0(substr(top_plot_data$term, 1, 97), "..."),
        top_plot_data$term
    )

    p <- ggplot2::ggplot(
        top_plot_data,
        ggplot2::aes(
            x = .data$cluster,
            y = .data$term,
            size = -log10(.data$adj.p_value)
        )
    ) +
        ggplot2::geom_point(aes(color = .data$odds_ratio), na.rm = TRUE) +
        ggplot2::geom_blank(aes(.data$cluster, .data$term)) +
        ggplot2::ylab("database: term") +
        ggplot2::scale_color_gradient(
            "odds\nratio",
            low = "blue",
            high = "red",
            labels = function(x) round(x, 2),
            guide = ggplot2::guide_colorbar(
                barheight = unit(12, "mm"),
                barwidth = unit(2, "mm"),
                ticks = FALSE
            )
        ) +
        ggplot2::scale_size_area(
            max_size = 3,
            limits = c(0, 2),
            breaks = c(0, 1, 2),
            labels = c("0", "1", "2 or higher"),
            oob = scales::oob_squish
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            panel.grid.major.x = ggplot2::element_blank(),
            panel.grid.minor.x = ggplot2::element_blank(),
            panel.grid.major.y = ggplot2::element_line(
                linewidth = 0.2,
                color = "gray80"
            ),
            panel.grid.minor.y = ggplot2::element_line(
                linewidth = 0.1,
                color = "gray90"
            ),
            plot.title = ggplot2::element_text(size = 12, hjust = 0.5),
            axis.text.y = ggplot2::element_text(size = 5),
            axis.text.x = ggplot2::element_text(size = 5),
            legend.text = ggplot2::element_text(size = 5),
            legend.title = ggplot2::element_text(size = 6),
            legend.key.height = unit(4, "mm"),
            legend.key.width = unit(3, "mm"),
            legend.spacing = ggplot2::unit(0, "mm"),
            legend.spacing.y = ggplot2::unit(0.5, "mm"),
            plot.margin = ggplot2::unit(c(2, 0, 2, 2), "mm"),
            legend.position = "right",
            legend.box.margin = ggplot2::margin(0, 0, 0, 0)
        ) +
        ggplot2::labs(title = title) +
        ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = c(1, 1))) +
        ggplot2::coord_cartesian(clip = "off")

    list(
        dotplot = p,
        dotplot_height = plot_height
    )
}


#' Flatten Nested ORA Results
#'
#' @noRd
#'
#' @description
#' Converts a nested ORA results list into a tidy long-format data frame. This
#' function extracts all result data frames from a structure like
#' \code{ora_results[[cluster]][[gene_set]]}, appends the cluster and gene set
#' as metadata columns, and returns a combined data frame.
#'
#' @param ora_results A nested list of ORA results. The first level of keys
#'   represents cluster names, and each sub-list contains result data frames
#'   from ORA against specific gene set databases.
#'
#' @return A data frame combining all ORA results, with additional columns:
#' \describe{
#'   \item{cluster}{The name of the cluster (character).}
#'   \item{gene_set}{The name of the gene set database.}
#'   \item{term}{A combined label of \code{gene_set: Description}.}
#' }
#'
#' Data frames that are \code{NULL}, \code{NA}, or have 0 rows are skipped.
#'
#' @importFrom dplyr bind_rows
#'
flatten_ora_results <- function(ora_results) {
    rows <- list()

    for (cluster in names(ora_results)) {
        for (gene_set in names(ora_results[[cluster]])) {
            df <- ora_results[[cluster]][[gene_set]]

            if (is.data.frame(df) && nrow(df) > 0) {
                df$cluster <- cluster
                df$gene_set <- gene_set
                df$term <- paste(gene_set, df$Description, sep = ": ")
                rows[[length(rows) + 1]] <- df
            }
        }
    }

    if (length(rows) == 0) {
        return(data.frame())
    }

    dplyr::bind_rows(rows)
}


#' Check overlap between foreground genes and a gene set map
#'
#' @noRd
#'
#' @description
#' This function compares a character vector of foreground genes (`fg`)
#' against the `gene` column of a gene set map (`gene_set_map`).
#' It reports the number and proportion of overlapping genes.
#'
#' @param fg Character vector of foreground genes.
#' @param gene_set_map Data frame with a column named `gene`.
#' @param verbose Boolean flag controlling the display of messages.
#'
#' @return Invisibly returns a list with counts and percentage overlap.
#'
check_gene_overlap <- function(
    fg,
    gene_set_map,
    verbose) {
    stopifnot(is.character(fg))
    stopifnot("gene" %in% colnames(gene_set_map))

    fg <- unique(fg[!is.na(fg) & nzchar(fg)])
    genes_avail <- unique(
        gene_set_map$gene[!is.na(gene_set_map$gene) & nzchar(gene_set_map$gene)]
    )

    overlap <- intersect(fg, genes_avail)
    n_fg <- length(fg)
    n_overlap <- length(overlap)
    pct_overlap <- if (n_fg == 0) 0 else round(100 * n_overlap / n_fg, 1)

    if (isTRUE(verbose)) message("Foreground genes:", n_fg)
    if (isTRUE(verbose)) {
        message(
            "Foreground genes overlapping with database: ",
            n_overlap,
            " (", pct_overlap, "%)\n"
        )
    }

    invisible(list(
        n_fg = n_fg,
        n_overlap = n_overlap,
        pct_overlap = pct_overlap
    ))
}


# Level 4 internal functions ---------------------------------------------------


#' Prepare Plot Data from Nested ORA Results
#'
#' @noRd
#'
#' @description
#' Flattens a nested ORA result list into a tidy data frame for plotting. The
#' input should be a list of clusters, each containing a list of gene set
#' names mapped to enrichment result data frames. This function extracts all
#' results, adds cluster and gene set metadata, and filters for terms that are
#' statistically significant and supported by at least two genes.
#'
#' @param ora_results A nested list of ORA results. Each top-level entry
#'   corresponds to a gene cluster (e.g., "cluster_1") and contains a named
#'   list of enrichment result data frames (one per gene set). Data frames
#'   must include the columns \code{p.adjust}, \code{Description}, \code{Count},
#'   and \code{odds_ratio}.
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{cluster}{Cluster identifier (character).}
#'   \item{term}{Combined gene set and term description label.}
#'   \item{adj.p_value}{Adjusted p-value for the enrichment term.}
#'   \item{odds_ratio}{Odds ratio of the enrichment.}
#' }
#'
#' Only rows with \code{p.adjust < 0.05} and \code{Count >= 2} are included in
#' the result. This output is intended for dotplot visualizations.
#'
#' @seealso \code{\link{add_odds_ratios_to_ora}},
#'   \code{\link{make_enrich_dotplot}}
#'
prepare_plot_data <- function(ora_results) {
    all_rows <- list()

    for (cluster in names(ora_results)) {
        for (gene_set in names(ora_results[[cluster]])) {
            df <- ora_results[[cluster]][[gene_set]]
            if (is.data.frame(df) && nrow(df) > 0) {
                df$cluster <- cluster
                df$gene_set <- gene_set
                df$term <- paste(
                    gene_set,
                    df$Description,
                    sep = ": "
                )
                all_rows[[length(all_rows) + 1]] <- df
            }
        }
    }

    full_df <- do.call(rbind, all_rows)

    # Filter for top plot: only supported by >= 2 genes
    top_df <- full_df[full_df$Count >= 2, ]
    top_df <- top_df[, c("cluster", "term", "p.adjust", "odds_ratio")]
    colnames(top_df) <- c(
        "cluster",
        "term",
        "adj.p_value",
        "odds_ratio"
    )

    top_df
}
