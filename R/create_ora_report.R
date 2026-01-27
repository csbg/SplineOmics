#' Generate an HTML report for over-representation analysis results
#'
#' @description
#' This function generates an HTML report for over-representation analysis
#' (ORA) results produced by \code{run_ora()}. It consumes a structured
#' \code{report_payload} object and creates dotplots summarizing enrichment
#' results for each grouping column in the cluster table.
#'
#' The function is responsible for all visualization and reporting side
#' effects. Statistical analysis is not recomputed. Calling this function
#' implies that an HTML report should be generated; therefore, a valid output
#' directory must be provided.
#'
#' @param report_payload `list`:
#'   A structured ORA report payload as returned by \code{run_ora()}. It must
#'   contain the raw enrichment results and all metadata required for report
#'   generation.
#'
#' @param report_info `list` | `NULL`:
#'   Optional list with experiment metadata used to annotate the HTML report
#'   (e.g., omics data type, project name, analyst, and data description). If
#'   \code{NULL}, a minimal report is generated.
#'
#' @param cluster_hits_report_name `character(1)` | `NULL`:
#'   Optional name of the \code{cluster_hits()} report that produced the
#'   clustering results used for ORA. When provided, it is displayed in the
#'   report to document provenance.
#'
#' @param report_dir `character(1)`:
#'   Directory where the HTML report and all associated output files are
#'   written. The directory is created if it does not already exist.
#'
#' @param verbose `logical(1)`:
#'   Logical flag controlling the display of progress and status messages.
#'
#' @return
#' A named list of \code{ggplot} objects containing the ORA dotplots generated
#' for the report. Each element is named after the corresponding grouping
#' column in the cluster table, allowing easy identification and reuse of the
#' plots outside the HTML report.
#'
#' @seealso
#' \code{run_ora()}, \code{cluster_hits()}, \code{clusterProfiler::enricher()}
#'
#' @export
#' 
create_ora_report <- function(
        report_payload,
        report_info = NULL,
        cluster_hits_report_name = NULL,
        verbose = TRUE,
        report_dir = here::here()
) {
    check_inputs_create_ora_report(
        report_payload = report_payload,
        report_info = report_info,
        cluster_hits_report_name = cluster_hits_report_name,
        verbose = verbose,
        report_dir = report_dir
    )
    
    all_results <- report_payload[["all_results"]]
    databases <- report_payload[["databases"]]
    clusterProfiler_params <- report_payload[["clusterProfiler_params"]]
    cluster_table <- report_payload[["cluster_table"]]
    universe <- report_payload[["universe"]]
    
    if (is.null(all_results) || length(all_results) == 0L) {
        stop_call_false(
            "No ORA results found in report_payload[['all_results']]."
        )
    }

    all_results <- setNames(lapply(names(all_results), function(nm) {

        ora_results <- all_results[[nm]]
        
        if (!is.list(ora_results)) {
            return(list(
                ora_results = NULL,
                dotplot = NULL,
                dotplot_nrows = NULL
            ))
        }
        
        any_result <- any(vapply(ora_results, function(cluster_entry) {
            any(vapply(cluster_entry, function(res) {
                is.data.frame(res) && nrow(res) > 0L
            }, logical(1)))
        }, logical(1)))
        
        if (!any_result) {
            if (isTRUE(verbose)) {
                message(
                    "No enrichment results for section: ", nm
                )
            }
            return(list(
                ora_results = ora_results,
                dotplot = NULL,
                dotplot_nrows = NULL
            ))
        }
        
        dp <- make_enrich_dotplot(
            ora_results = ora_results,
            title = nm
        )
        
        list(
            ora_results = ora_results,
            dotplot = dp[["dotplot"]],
            dotplot_nrows = dp[["dotplot_height"]]
        )
    }), names(all_results))
    
    processed_results <- Map(
        function(res, nm) process_result(res, nm),
        all_results,
        names(all_results)
    )
    
    sections <- lapply(seq_along(processed_results), function(i) {
        pr <- processed_results[[i]]
        raw <- all_results[[i]]
        
        plots <- pr$plot
        if (is.null(plots)) {
            plots <- list()
        } else if (inherits(plots, "ggplot") || !is.list(plots)) {
            plots <- list(plots)
        }
        
        sizes <- pr$plot_size
        if (is.null(sizes)) {
            sizes <- as.list(rep_len(NA_real_, length(plots)))
        } else if (!is.list(sizes)) {
            sizes <- as.list(rep_len(sizes, length(plots)))
        }
        
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
    
    report_info$databases <- unique(databases[["DB"]])
    report_info$clusterProfiler_params <- clusterProfiler_params
    report_info$cluster_hits_report_name <- cluster_hits_report_name
    report_info$cluster_table <- cluster_table
    report_info$background_genes <- universe
    
    total_plots <- sum(vapply(
        sections,
        function(s) length(s$plots),
        integer(1)
    ))
    
    if (total_plots == 0L) {
        stop_call_false(
            "No enrichment results to plot; report not generated."
        )
    }
    
    report_dir <- normalizePath(
        report_dir,
        mustWork = FALSE
    )
    
    generate_report_html(
        plots = sections,
        report_info = report_info,
        report_type = "run_ora",
        filename = "run_ora_report",
        report_dir = report_dir
    )
    
    if (isTRUE(verbose)) {
        print_info_message(
            message_prefix = "ORA analysis",
            report_dir = report_dir
        )
    }
    
    setNames(lapply(sections, function(s) {
        if (length(s$plots) == 0L) {
            return(NULL)
        }
        s$plots[[1]]
    }), names(sections))
}


# Level 1 internal functions ---------------------------------------------------


#' Validate inputs for create_ora_report()
#' 
#' @noRd
#'
#' @description
#' Internal helper function that validates the inputs passed to
#' \code{create_ora_report()}. It performs structural checks on the ORA report
#' payload, ensures that required arguments are present and correctly typed,
#' and enforces constraints implied by report generation (e.g., a mandatory
#' output directory).
#'
#' The function performs shallow validation only: it checks the presence and
#' names of required payload fields but does not inspect their internal
#' contents.
#'
#' @param report_payload `SplineOmicsOraReportPayload`:
#'   Report payload object as returned by \code{run_ora()}. It must inherit from
#'   class \code{"SplineOmicsOraReportPayload"} and contain exactly the expected
#'   fields.
#'
#' @param report_info `list` | `NULL`:
#'   Optional named list with experiment metadata used for annotating the ORA
#'   report.
#'
#' @param cluster_hits_report_name `character(1)` | `NULL`:
#'   Optional name of the \code{cluster_hits()} report whose results were used
#'   for the ORA.
#'
#' @param verbose `logical(1)`:
#'   Logical flag controlling the display of validation messages.
#'
#' @param report_dir `character(1)`:
#'   Directory path where the HTML report will be written. Must be a non-empty
#'   character string.
#'
#' @return
#' Invisibly returns \code{TRUE} if all input checks pass. The function throws
#' an error with an informative message if any validation fails.
#'
check_inputs_create_ora_report <- function(
        report_payload,
        report_info,
        cluster_hits_report_name,
        verbose,
        report_dir
) {
    # report_payload: class + expected fields only (no deep validation)
    if (!inherits(report_payload, "SplineOmicsOraReportPayload")) {
        stop_call_false(
            paste0(
                "`report_payload` must inherit from class ",
                "'SplineOmicsOraReportPayload'."
            )
        )
    }
    
    expected_fields <- c(
        "all_results",
        "cluster_table",
        "databases",
        "clusterProfiler_params",
        "enrichGO_cfg",
        "mapping_cfg",
        "universe"
    )
    missing_fields <- setdiff(expected_fields, names(report_payload))
    extra_fields <- setdiff(names(report_payload), expected_fields)
    
    if (length(missing_fields) > 0L) {
        stop_call_false(
            paste(
                "`report_payload` is missing required field(s):",
                paste(missing_fields, collapse = ", ")
            )
        )
    }
    if (length(extra_fields) > 0L) {
        stop_call_false(
            paste(
                "`report_payload` contains unexpected field(s):",
                paste(extra_fields, collapse = ", ")
            )
        )
    }
    
    # report_dir: required and must be a single, non-empty string
    if (is.null(report_dir) || !is.character(report_dir) ||
        length(report_dir) != 1L || !nzchar(report_dir)) {
        stop_call_false(
            "`report_dir` must be a non-empty character(1) path."
        )
    }
    
    # cluster_hits_report_name: NULL or character(1)
    if (!is.null(cluster_hits_report_name) &&
        (!is.character(cluster_hits_report_name) ||
         length(cluster_hits_report_name) != 1L ||
         !nzchar(cluster_hits_report_name))) {
        stop_call_false(
            paste(
                "`cluster_hits_report_name` must be a non-empty",
                "character(1) string or NULL."
            )
        )
    }
    
    # verbose: logical(1)
    if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
        stop_call_false("`verbose` must be a logical(1) value.")
    }
    
    # report_info: NULL or named list
    if (!is.null(report_info)) {
        if (!is.list(report_info) || is.null(names(report_info))) {
            stop_call_false("`report_info` must be a named list or NULL.")
        }
        if (any(!nzchar(names(report_info)))) {
            stop_call_false("`report_info` must not contain empty names.")
        }
    }
    
    invisible(TRUE)
}


# Level 2 internal functions ---------------------------------------------------


# Level 3 internal functions ---------------------------------------------------

