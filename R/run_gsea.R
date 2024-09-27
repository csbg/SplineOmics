# Export functions definitions -------------------------------------------------


#' Generate a GSEA Report
#'
#' @description
#' This function generates a Gene Set Enrichment Analysis (GSEA) report based
#' on clustered hit levels, gene data, and specified databases. It processes
#' the input data, manages GSEA levels, and produces an HTML report with plots.
#'
#' @param levels_clustered_hits A list of clustered hits at different levels.
#' @param databases A list of databases for the gene set enrichment analysis.
#' @param report_info A list containing information for the report generation.
#' @param clusterProfiler_params Additional parameters for the GSEA analysis,
#'                               default is NA. Those include adj_p_value,
#'                               pAdjustMethod, etc (see clusterProfiler
#'                               documentation).
#' @param plot_titles Titles for the plots, default is NA.
#' @param universe Enrichment background data, default is NULL.
#' @param report_dir Directory where the report will be saved, default is
#' `here::here()`.
#'
#' @return A list of plots generated for the GSEA report.
#'
#' @importFrom purrr map2 flatten
#' @importFrom here here
#'
#' @export
#'
run_gsea <- function(
    levels_clustered_hits,
    databases,
    report_info,
    clusterProfiler_params = NA,
    plot_titles = NA,
    universe = NULL,
    report_dir = here::here()
    ) {

  report_dir <- normalizePath(
    report_dir,
    mustWork = FALSE
  )

  # Check report_info and report_dir
  args <- lapply(
    as.list(
      match.call()[-1]),
    eval,
    parent.frame()
    )

  input_control <- InputControl$new(args)
  input_control$auto_validate()

  # Remove levels that
  levels_clustered_hits <- levels_clustered_hits[
    !sapply(
      levels_clustered_hits,
      is.character
      )
    ]

  # Control the test not covered by the InputControl class
  control_inputs_create_gsea_report(
    levels_clustered_hits = levels_clustered_hits,
    databases = databases,
    params = clusterProfiler_params,
    plot_titles = plot_titles,
    background = universe
    )

  ensure_clusterProfiler()  # Deals with clusterProfiler installation.

  all_results <- map2(
    levels_clustered_hits,
    names(levels_clustered_hits),
    ~ manage_gsea_level(
      .x,
      .y,
      databases,
      clusterProfiler_params,
      universe
      )
  )

  names(all_results) <- names(levels_clustered_hits)

  processed_results <- map2(
    all_results,
    names(all_results),
    process_result
  )

  # Extract the plots, plot sizes, and header info from the processed results
  plots <- purrr::flatten(map(processed_results, "plot"))
  plots_sizes <- unlist(map(processed_results, "plot_size"))
  
  insert_after_each <- function(lst, value) {
    result <- vector("list", 2 * length(lst))  # Create a list twice the size
    result[seq(1, length(result), by = 2)] <- lst  # Insert original elements at odd positions
    result[seq(2, length(result), by = 2)] <- value  # Insert value at even positions
    return(result)
  }
  
  plots <- insert_after_each(plots, "section_break")
  plots_sizes <- insert_after_each(plots_sizes, 999)

  level_headers_info <- map(processed_results, "header_info")

  names(level_headers_info) <- names(all_results)


  report_info$databases <- unique(databases[["DB"]])

  generate_report_html(
    plots = plots,
    plots_sizes = plots_sizes,
    report_info = report_info,
    level_headers_info = level_headers_info,
    report_type = "create_gsea_report",
    filename = "create_gsea_report",
    report_dir = report_dir
    )

  print_info_message(
    message_prefix = "Gene set enrichment analysis",
    report_dir = report_dir
  )

  return(Filter(function(x) !is.character(x), plots))
}



# Level 1 internal functions ---------------------------------------------------


#' Control Inputs for GSEA Report
#'
#' @description
#' Validates the inputs for generating a GSEA report, including clustered
#' hits, genes, databases, parameters, plot titles, and background genes.
#'
#' @param levels_clustered_hits A list containing clustered hits at various
#' levels.
#' @param databases A list of databases to be used in the GSEA analysis.
#' @param params A list of parameters for the GSEA analysis.
#' @param plot_titles A character vector of titles for the plots, with length
#' matching `levels_clustered_hits`.
#' @param background A character vector of background genes or NULL.
#'
control_inputs_create_gsea_report <- function(
    levels_clustered_hits,
    databases,
    params,
    plot_titles,
    background
    ) {

  check_clustered_hits(levels_clustered_hits)

  check_databases(databases)

  check_params(params)


  if (!is.na(plot_titles)) {

    if (!is.character(plot_titles) ||
        length(plot_titles) != length(levels_clustered_hits)) {
      stop(paste("plot_titles must be a character vector with length == length",
                 "length levels_clustered_hits"), call. = FALSE)
    }
  }


  if (!is.null(background)) {
    if (!is.character(background)) {
      stop("background must be a character vector or NULL.", call. = FALSE)
    } else {
      check_genes(background)
    }
  }
}


#' Ensure 'clusterProfiler' is installed and loaded
#'
#' @description
#' This function checks if the 'clusterProfiler' package is installed.
#' If not, it prompts the user to choose whether to install it automatically,
#' install it manually, or cancel the operation. Once installed, the package
#' is loaded for use.
#'
ensure_clusterProfiler <- function() {

  # Check if clusterProfiler is installed; if not, prompt the user
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    message("The 'clusterProfiler' package is not installed.")

    # Prompt user for action
    repeat {
      user_input <- readline(
        prompt = paste0(
          "What would you like to do?\n",
          "1: Automatically install clusterProfiler\n",
          "2: Manually install clusterProfiler\n",
          "3: Cancel\n",
          "Please enter 1, 2, or 3: "
        )
      )

      if (user_input == "1") {
        # Try to install clusterProfiler automatically from Bioconductor
        message(
          "Attempting to install 'clusterProfiler' automatically
          from Bioconductor..."
          )

        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          utils::install.packages("BiocManager")
        }

        tryCatch(
          {
            BiocManager::install("clusterProfiler", update = FALSE, ask = FALSE)
          },
          error = function(e) {
            stop(
              "Automatic installation of 'clusterProfiler' failed.
              Please install it manually and try again.",
              call. = FALSE
              )
          }
        )
        break  # Exit the loop if installation is successful
      } else if (user_input == "2") {
        stop(
          "Please install 'clusterProfiler' manually
          using BiocManager::install('clusterProfiler') and
          then re-run the function.",
          call. = FALSE
          )
      } else if (user_input == "3") {
        stop("Operation canceled by the user.", call. = FALSE)
      } else {
        message("Invalid input. Please enter 1, 2, or 3.")
      }
    }
  }
}


#' Manage GSEA Analysis for a Specific Level
#'
#' @description
#' This function manages the GSEA analysis for a specific level. It extracts
#' genes associated with the clustered hits, removes rows with `NA` values,
#' and runs the GSEA analysis using the `create_gsea_report` function.
#'
#' @param clustered_hits A dataframe containing the clustered hits for a
#' specific level. It must include a column named `feature` to extract genes.
#' @param level_name A character string representing the name of the level.
#' @param databases A list of databases for the gene set enrichment analysis.
#' @param clusterProfiler_params Additional parameters for the GSEA analysis,
#'                               default is NA. Those include adj_p_value,
#'                               pAdjustMethod, etc (see clusterProfiler
#'                               documentation).
#' @param universe Enrichment background data, default is NULL.
#'
#' @return The result of the `create_gsea_report` function, which typically
#'         includes various plots and enrichment results.
#'
manage_gsea_level <- function(
    clustered_hits,
    level_name,
    databases,
    clusterProfiler_params,
    universe
    ) {

  clustered_hits <- na.omit(clustered_hits)

  message(paste(
    "\n\n Running clusterProfiler for the level:",
    level_name
    ))

  result <- create_gsea_report_level(
    clustered_genes = clustered_hits,
    databases = databases,
    params = clusterProfiler_params,
    plot_title = level_name,
    universe = universe
    )
}


#' Process GSEA Result for a Specific Level
#'
#' @description
#' This function processes the GSEA result for a specific level. It handles
#' cases where the result contains `NA` values by adding a section break.
#' Otherwise, it extracts the plot, plot size, and header information from
#' the result.
#'
#' @param level_result A list containing the GSEA result for a specific level.
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
    level_name
    ) {

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
      full_enrich_results = level_result$full_enrich_results,
      raw_enrich_results = level_result$raw_results
    )
  }

  return(result)
}


#' Build GSEA Report
#'
#' @description
#' Generates an HTML report for Gene Set Enrichment Analysis (GSEA) based on
#' provided plot data, header information, and other content. The report
#' includes sections for each level of clustered hits, along with a table of
#' contents and various plots.
#'
#' @param header_section A string containing the HTML content for the header
#' section of the report.
#' @param plots A list of plots to be included in the report.
#' @param plots_sizes A list of sizes for the plots.
#' @param level_headers_info A list containing header information for each
#' level of clustered hits.
#' @param report_info A named list containg the report info fields. Here used
#'                    for the email hotkey functionality.
#' @param output_file_path A string specifying the file path where the report
#' will be saved.
#'
#' @return None. The function generates and writes an HTML report to the
#' specified output file path.
#'
#' @details
#' The function first initializes the HTML content with the provided header
#' section and a placeholder for the table of contents (TOC). It then iterates
#' through the plots, generating sections for each level of clustered hits and
#' processing individual plots. The TOC is inserted into the HTML content,
#' which is then finalized and written to the specified output file.
#'
build_create_gsea_report <- function(
    header_section,
    plots,
    plots_sizes,
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
  level_headers_info <- Filter(
    Negate(is.null),
    level_headers_info
    )

  pb <- create_progress_bar(plots)
  # Generate the sections and plots
  for (index in seq_along(plots)) {

    # means jump to next level
    if (any(class(plots[[index]]) == "character")) {

      section_info <- level_headers_info[[current_header_index]]

      section_content <- generate_section_content(
        section_info,
        index,
        toc,
        html_content,
        section_header_style,
        toc_style
      )

      html_content <- section_content$html_content
      # toc <- section_content$toc

      current_header_index <- current_header_index + 1

      pb$tick()
      next
    }
    
    # Add the section header and horizontal line just before the plot
    section_info <- level_headers_info[[current_header_index]]
    section_header <- sprintf(
      "<h2 style='%s' id='section%d'>%s</h2>",
      section_header_style,
      index,
      section_info$header_name
    )
    
    horizontal_line <- ""
    
    if (current_header_index > 1) {
      horizontal_line <- "<hr>"
    }
    
    # Update the HTML content with the section header and horizontal line
    html_content <- paste(
      html_content,
      horizontal_line,
      section_header,
      sep = "\n"
    )
    
    toc_entry <- sprintf(
      "<li style='%s'><a href='#section%d'>%s</a></li>",
      toc_style,
      index,
      section_info$header_name
    )
    
    toc <- paste(
      toc,
      toc_entry,
      sep = "\n"
      )

    result <- process_plots(
      plots_element = plots[[index]],
      plots_size = plots_sizes[[index]],
      html_content = html_content,
      toc = toc,
      header_index = index
    )
    html_content <- result$html_content
    toc <- result$toc

    pb$tick()
  }

  generate_and_write_html(
    toc = toc,
    html_content = html_content,
    report_info,
    output_file_path = output_file_path
  )
}



# Level 2 internal functions ---------------------------------------------------


#' Check Clustered Genes Dataframe for Required Conditions
#'
#' @description
#' This function checks if a given dataframe `clustered_genes` contains the
#' required columns `gene` and `cluster`. The `gene` column must contain only
#' character strings of length 1, and the `cluster` column must contain only
#' integers. If any condition is not met, the function stops the script and
#' produces an informative error message.
#'
#' @param levels_clustered_hits A list of dataframes to be checked for the
#'                              required format.
#'
#' @return This function does not return a value. It stops with an error message
#'         if the conditions are not met.
#'
check_clustered_hits <- function(levels_clustered_hits) {

  if (!is.list(levels_clustered_hits)) {
    stop(paste("levels_clustered_hits must be a list"), call. = FALSE)
  }

  for (i in seq_along(levels_clustered_hits)) {

    clustered_hits <- levels_clustered_hits[[i]]

    # Check if clustered_hits is a dataframe
    if (!is.data.frame(clustered_hits)) {
      stop(paste("Element", i ,"of levels_clustered_hits is not a dataframe",
                 "but must be one."), call. = FALSE)
    }

    # Check if the dataframe contains the columns 'gene' and 'cluster'
    required_columns <- c("feature", "cluster")
    if (!all(required_columns %in% colnames(clustered_hits))) {
      stop(paste("The dataframe must contain the columns 'feature' and",
                 "'cluster'."),
           call. = FALSE)
    }

    # Check if the 'feature' column contains only integers
    if (!is.integer(clustered_hits$feature) &&
        !all(clustered_hits$feature == as.integer(clustered_hits$feature))) {
      stop("The 'feature' column must contain only integers.", call. = FALSE)
    }

    # Check if the 'cluster' column contains only integers
    if (!is.integer(clustered_hits$cluster) &&
        !all(clustered_hits$cluster == as.integer(clustered_hits$cluster))) {
      stop("The 'cluster' column must contain only integers.", call. = FALSE)
    }
  }
}


#' Check Valid Gene IDs
#'
#' @description
#' This function checks whether a character vector `genes`
#' contains only valid gene IDs. Each gene ID must consist
#' solely of alphabetic letters and numbers.
#'
#' @param genes A character vector containing gene IDs.
#' @param max_index_overall An integer, specifying the highest index of all
#'                          features across all levels.
#'
#' @return None. This function stops execution and provides
#' an error message if the vector does not meet the criteria,
#' including the first offending element and its index.
#'
check_genes <- function(
    genes,
    max_index_overall = NA
    ) {

  if (!is.na(max_index_overall)) {
    if (length(genes) < max_index_overall) {
      stop(paste("genes must at least have over", max_index_overall,
                 "elements"),
           call. = FALSE)
    }
  }


  if (!is.character(genes)) {
    stop("genes must be a character vector", call. = FALSE)
  }

  valid <- grepl("^[a-zA-Z0-9]+$", genes) | is.na(genes)
  if (any(!valid)) {
    first_invalid_index <- which(!valid)[1]
    first_invalid_value <- genes[first_invalid_index]
    num_invalid <- sum(!valid) - 1
    stop(sprintf(paste("Invalid gene found at index %d: '%s'.",
                 "There are %d more invalid elements."),
                 first_invalid_index, first_invalid_value, num_invalid),
         call. = FALSE)
  }
}


#' Check Valid Databases Dataframe
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
         call. = FALSE)
  }

  if (!all(sapply(databases, is.character))) {
    stop("All columns in the dataframe must be of type character.",
         call. = FALSE)
  }
}


#' Check Params List for Required Conditions
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
             expected_type, ".", call. = FALSE)
      }
    } else {
      stop(
        "Unexpected element '",
        param,
        "' in the list.",
        call. = FALSE
        )
    }

    if (param == 'pAdjustMethod') {
      if (!(actual_value %in% valid_adj_p_value_methods)) {
        stop(paste("pAdjustMethod must be one of",
                   valid_adj_p_value_methods, collapse = ", "),
             call. = FALSE)

      }
    }

  }
}


#' Perform Gene Set Enrichment Analysis and plot it.
#'
#' This function conducts a Gene Set Enrichment Analysis (GSEA) using either the
#' clusterProfiler package. Afterwards, it plots the results.
#' It allows for customization of enrichment parameters, selection of databases,
#' and optionally specifying a custom plot title and background gene list.
#'
#' @param clustered_genes A list of dataframes with two columns: the first
#' column contains the standard gene symbol, and the second column contains an
#' integer specifying the cluster.
#' @param databases A dataframe containing the data of the downloaded Enrichr
#'                  databases
#' @param params A list specifying the clusterProfiler parameters for the
#'               enrichment analysis.
#' @param plot_title An optional string specifying the title of the plot. If not
#' provided, a default title based on the analysis will be used.
#' @param universe An optional list of standard gene symbols to be used as the
#' background for the enrichment analysis instead of the background chosen by
#' the `enricher`. The default is an empty list, which implies the use of the
#' default background set by the enrichment tool.
#'
#' @return An object containing the results of the Gene Set Enrichment Analysis,
#' including any plots generated during the analysis.
#'
create_gsea_report_level <- function(
    clustered_genes,
    databases,
    params = NA,
    plot_title = "",
    universe = NULL
    ) {

  set_default_params(params)

  all_term2genes <- dbs_to_term2genes(databases)


  ## Prepare objects
  unique_clusters <- sort(unique(clustered_genes$cluster))

  all_db_results <- vector("list", length(all_term2genes))

  for (i in seq_along(all_db_results)) {
    all_db_results[[i]] <- data.frame(BioProcess = character())
  }

  raw_results <- list()

  for (cluster in unique_clusters) {
    enrichment_results <- list()

    for (i in seq_along(all_db_results)) {
      column_name <- paste("Cluster", cluster, sep = "_")
      count_column_name <- paste0(column_name, "_odds_ratios")

      # Initialize the column as an empty vector for the df if it has no rows
      if (nrow(all_db_results[[i]]) == 0) {
        all_db_results[[i]][[column_name]] <- vector("logical", length = 0)
        all_db_results[[i]][[count_column_name]] <-
          vector("logical", length = 0)
      } else {
        all_db_results[[i]][[column_name]] <- NA
        all_db_results[[i]][[count_column_name]] <- NA
      }
    }

    # Select rows where 'Cluster' equals the current cluster number
    selected_rows <- clustered_genes[clustered_genes$cluster == cluster, ]

    cluster_genes <- selected_rows$gene
    gene_list <- as.list(cluster_genes)
    gene_list <- unlist(gene_list)

    at_least_one_result = FALSE

    # Run clusterProfiler and process output
    for (database_name in names(all_term2genes)) {
      term2gene <- all_term2genes[[database_name]]

      message(paste(
        "\nDatabase:",
        database_name
        ))

      enrichment <-
        clusterProfiler::enricher(
          gene = gene_list,
          pvalueCutoff  = params$pvalueCutoff,
          pAdjustMethod = params$pAdjustMethod,
          universe      = universe,
          minGSSize     = params$minGSSize,
          maxGSSize     = params$maxGSSize,
          qvalueCutoff  = params$qvalueCutoff,
          gson          = NULL,
          TERM2GENE     = term2gene,
          TERM2NAME     = NA
          )

      enrichment <- as.data.frame(enrichment)

      if (is.null(enrichment) || (nrow(enrichment) == 0 &
                                  ncol(enrichment) == 0)) {
        next
      }

      at_least_one_result = TRUE

      enrichment_results[[length(enrichment_results) + 1]] <- enrichment

      # Store all for returning in the end.
      name <- sprintf(
        "cluster: %s, database: %s",
        cluster,
        database_name
        )

      raw_results[[name]] <- enrichment
    }

    if (is.null(universe)) {
      use_background = TRUE
    } else {
      use_background = FALSE
    }

    if (at_least_one_result) {
      all_db_results <- process_enrichment_results(
        all_db_results,
        enrichment_results,
        params$pvalueCutoff,
        column_name,
        count_column_name,
        background = use_background
      )
    } else {
      print(paste0("No enrichment results for cluster", cluster))
    }
  }

  # Make dotplot
  any_result <- sapply(all_db_results, function(df) nrow(df) > 0)
  has_true <- any(any_result)

  if (has_true) {
    result <- make_enrich_dotplot(
      all_db_results,
      names(all_term2genes),
      plot_title
      )
  } else {
    message("No database led to an enrichment result!")
    return(NA)
  }

  list(
    dotplot = result[["dotplot"]],
    dotplot_nrows = result[["dotplot_height"]],
    full_enrich_results = result[["full_enrich_results"]],
    raw_results = raw_results
    )
}


#' Generate Section Content
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
generate_section_content <- function(
    section_info,
    index,
    toc,
    html_content,
    section_header_style,
    toc_style
) {

  if (any(is.na(section_info$full_enrich_results))) {

    no_results_message <- paste0(
      "<p style='font-size: 40px; color: #FF0000;'>",
      "No database specified led to identification of any enrichment ",
      "terms.",
      "</p>"
    )

    html_content <- paste(
      html_content,
      no_results_message,
      sep = "\n"
    )

    return(list(
      html_content = html_content
    ))
  }

  full_enrich_results_header <- paste0(
    "<h3 style='font-size: 30px; font-weight: bold; color: #333;",
    "'>Enrichment Results</h3>"
    )

  df <- section_info$full_enrich_results

  # Start the HTML table with the specified attributes
  html_table <- "<table style='width:100%;border-collapse:collapse;'>"

  # Add the table header
  html_table <- paste0(html_table, "<thead><tr>")
  for (header in colnames(df)) {
    html_table <- paste0(html_table, "<th>", header, "</th>")
  }
  html_table <- paste0(html_table, "</tr></thead><tbody>")

  # Add the table rows
  for (i in 1:nrow(df)) {
    html_table <- paste0(html_table, "<tr>")
    for (j in 1:ncol(df)) {
      html_table <- paste0(html_table, "<td>", df[i, j], "</td>")
    }
    html_table <- paste0(html_table, "</tr>")
  }

  # Close the table body and the table tag
  full_enrich_results_html <- paste0(html_table, "</tbody></table>")

  raw_enrich_results_header <- paste0(
    "<h3 style='font-size: 30px; font-weight: bold; color: #333;",
    "'>Count smaller 2 Enrichment Results</h3>"
    )

  base64_df <- sprintf(
    '<a href="%s" download="count2small_results.xlsx">
  <button>Download count2small_results.xlsx</button></a>',
    encode_df_to_base64(
      section_info$raw_enrich_results,
      "create_gsea_report"
    )
  )

  html_content <- paste(
    html_content,
    full_enrich_results_header,
    full_enrich_results_html,
    raw_enrich_results_header,
    base64_df,
    sep = "\n"
  )

  list(
    html_content = html_content
  )
}


# Level 3 internal functions ---------------------------------------------------


#' Set Default Parameters
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


#' Convert Database File to TERM2GENE List
#'
#' Reads a specified .tsv file containing information about databases,
#' gene sets, and genes. The file should have three columns: 'DB' for database
#' names, Geneset' for gene set identifiers, and 'Gene' for gene names. This
#' function organizes this information into a nested list. Each top-level
#' element corresponds to a unique database, and within each, gene sets map to
#' lists of associated genes.
#'
#' @param databases A dataframe, containing the three columns DB, Geneset, and
#'                  gene. This dataframe contains the databases downloaded from
#'                  Enrichr with the SplineOmics package function:
#'                  download_enrichr_databases.
#' @return A nested list where the first level of names corresponds to database
#' names ('DB'),
#'         the second level to gene sets ('Geneset'), and the innermost lists
#'         contain gene names ('Gene') associated with each gene set.
#'
dbs_to_term2genes <- function(databases) {

  db_split <- split(databases, databases$DB)

  # Transform into long format
  all_term2genes <- lapply(db_split, function(db_df) {

    df_with_renamed_columns <- db_df[, c("Geneset", "Gene")]
    colnames(df_with_renamed_columns) <- c("term", "gene")
    return(df_with_renamed_columns)
  })

  names(all_term2genes) <- names(db_split)

  return(all_term2genes)
}


#' Process Enrichment Results
#'
#' Process enrichment results for visualization.
#'
#' @param all_db_results A list of data frames containing enrichment results
#'                       for all databases.
#' @param enrichment_results A list of data frames containing enrichment
#'                           results for individual databases.
#' @param adjP_threshold The threshold for adjusted p-values.
#' @param column_name The name of the column to store adjusted p-values.
#' @param count_column_name The name of the column to store gene counts.
#' @param background Logical indicating whether background ratios are included.
#'
#' @return A list of data frames containing processed enrichment results.
#'
process_enrichment_results <- function(
    all_db_results,
    enrichment_results,
    adjP_threshold,
    column_name,
    count_column_name,
    background = FALSE
    ) {

  column_indices <- list(2, 6, 3, 4, 9)

  # Process results for all databases.
  for (i in seq_along(enrichment_results)) {
    df <- enrichment_results[[i]]
    df <- subset(df, df[[column_indices[[2]]]] < adjP_threshold)

    if (nrow(df) == 0) {
      next
    }

    term_list <- as.list(df[[column_indices[[1]]]])
    adjP_list <- as.list(df[[column_indices[[2]]]])
    odds_ratio <- as.list(df[[column_indices[[3]]]])

    odds_ratio <- sapply(odds_ratio, function(x) eval(parse(text = x)))
    bg_ratio <- as.list(df[[column_indices[[4]]]])
    bg_ratio <- sapply(bg_ratio, function(x) eval(parse(text = x)))
    odds_ratio <- mapply("/", odds_ratio, bg_ratio)

    gene_count_list <- as.list(df[[column_indices[[5]]]])

    named_list <- list()

    # Loop through the terms
    for(j in seq_along(term_list)) {
      # Create a sublist for each term with its adjP value and gene count

      # Skip terms that are just supported by one gene.
      if (gene_count_list[j] < 2) {
        next
      }

      sublist <- list(adjP_list[j], odds_ratio[j])

      # Assign this sublist to named_list with the term as its name
      term <- as.character(term_list[j])
      named_list[[term]] <- sublist
    }

    for (name in names(named_list)) {

      if (!name %in% all_db_results[[i]]$BioProcess) {
        # Create a list with the same structure as your DataFrame
        row_index <- nrow(all_db_results[[i]]) + 1
        all_db_results[[i]][row_index, "BioProcess"] <- name
        all_db_results[[i]][row_index, column_name] <- named_list[[name]][[1]]
        all_db_results[[i]][row_index, count_column_name] <-
          named_list[[name]][[2]]

      } else {
        row_index <- which(all_db_results[[i]]$BioProcess == name)
        all_db_results[[i]][row_index, column_name] <- named_list[[name]][[1]]
        all_db_results[[i]][row_index, count_column_name] <-
          named_list[[name]][[2]]
      }
    }
  }
  return(all_db_results)
}


#' Make Enrich Dotplot
#'
#' Make an enriched dotplot for visualization.
#'
#' @param enrichments_list A list of enrichments containing data frames for
#'                         different databases.
#' @param databases A character vector specifying the databases to be included.
#' @param title A character string specifying the title of the dotplot.
#'
#' @return A list containing:
#' \describe{
#'   \item{p}{The ggplot object representing the dotplot.}
#'   \item{dotplot_nrows}{An integer specifying the number of rows in the
#'                        dotplot.}
#'   \item{full_enrich_results}{A data frame containing the full enrichments
#'                              results.}
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point ylab coord_fixed guide_colorbar
#'                     theme_bw theme element_blank element_text labs
#' @importFrom scales oob_squish
#' @importFrom rlang .data
#'
# make_enrich_dotplot <- function(
#     enrichments_list,
#     databases,
#     title = "Title"
# ) {
#
#   results <- prepare_plot_data(
#     enrichments_list,
#     databases
#   )
#
#   top_plot_data <- results$top_plot_data
#   full_enrich_results <- results$full_enrich_results
#
#   # Calculate plot height based on the number of y-axis labels
#   height_per_label <- 0.13
#   num_labels <- length(unique(top_plot_data$term))
#   plot_height <- num_labels * height_per_label
#
#   if (plot_height > 0.9) {  # to always have a minimum size.
#     plot_height <- 0.9
#   }
#
#   # Ensure term labels are truncated to a maximum of 100 characters
#   top_plot_data$term <- as.character(top_plot_data$term)
#   top_plot_data$term <- ifelse(
#     nchar(top_plot_data$term) > 100,
#     paste0(substr(top_plot_data$term, 1, 97), "..."),
#     top_plot_data$term
#   )
#
#   p <- ggplot2::ggplot(
#     top_plot_data,
#     ggplot2::aes(
#       .data$cluster,
#       .data$term,
#       size = -log10(.data$adj.p_value)
#     )
#   ) +
#     ggplot2::geom_point(
#       aes(color = .data$odds_ratios),
#       na.rm = TRUE
#     ) +
#     ggplot2::geom_blank(aes(
#       .data$cluster,
#       .data$term
#     )) + # Ensure all columns are shown
#     ggplot2::ylab("database: term") +
#     ggplot2::scale_color_gradient(
#       "odds\nratio",
#       low = "blue",
#       high = "red",
#       labels = function(x) round(x, 2),
#       guide = ggplot2::guide_colorbar(
#         barheight = unit(12, "mm"),
#         barwidth = unit(2, "mm"),
#         ticks = FALSE
#       )
#     ) +
#     ggplot2::scale_size_area(
#       max_size = 3,
#       limits = c(0, 2),
#       breaks = c(0, 1, 2),
#       labels = c("0", "1", "2 or higher"),
#       oob = scales::oob_squish
#     ) +
#     ggplot2::theme_bw() +
#     ggplot2::theme(
#       panel.grid = ggplot2::element_blank(),
#       plot.title = ggplot2::element_text(size = 12, hjust = 0.5),
#       axis.text.y = ggplot2::element_text(size = 5),
#       axis.text.x = ggplot2::element_text(size = 5),
#       legend.text = ggplot2::element_text(size = 5),
#       legend.title = ggplot2::element_text(size = 6),
#       legend.key.height = unit(4, "mm"),
#       legend.key.width = unit(3, "mm"),
#       legend.spacing = ggplot2::unit(4, "mm"),
#       plot.margin = ggplot2::unit(c(2, 2, 2, 2), "mm")
#     ) +
#     ggplot2::labs(title = title) +
#     ggplot2::scale_y_discrete(expand = expansion(mult = c(0.01, 0.01))) +
#     ggplot2::coord_cartesian(clip = "off")
#
#   list(
#     dotplot = p,
#     dotplot_height = plot_height,
#     full_enrich_results = full_enrich_results
#   )
# }
make_enrich_dotplot <- function(
    enrichments_list,
    databases,
    title = "Title"
) {

  results <- prepare_plot_data(
    enrichments_list,
    databases
  )

  top_plot_data <- results$top_plot_data
  full_enrich_results <- results$full_enrich_results

  # Calculate plot height based on the number of y-axis labels
  height_per_label <- 0.1
  num_labels <- length(unique(top_plot_data$term))
  plot_height <- num_labels * height_per_label

  if (plot_height < 0.70) {  # to always have a minimum size.
    plot_height <- 0.70
  }

  # Ensure term labels are truncated to a maximum of 100 characters
  top_plot_data$term <- as.character(top_plot_data$term)
  top_plot_data$term <- ifelse(
    nchar(top_plot_data$term) > 100,
    paste0(substr(top_plot_data$term, 1, 97), "..."),
    top_plot_data$term
  )

  p <- ggplot2::ggplot(
    top_plot_data,
    ggplot2::aes(
      .data$cluster,
      .data$term,
      size = -log10(.data$adj.p_value)
    )
  ) +
    ggplot2::geom_point(
      aes(color = .data$odds_ratios),
      na.rm = TRUE
    ) +
    ggplot2::geom_blank(aes(
      .data$cluster,
      .data$term
    )) + # Ensure all columns are shown
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
      panel.grid.major.y = ggplot2::element_line(size = 0.2, color = "gray80"),
      panel.grid.minor.y = ggplot2::element_line(size = 0.1, color = "gray90"),
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
    ggplot2::scale_y_discrete(expand = expansion(mult = c(0.1, 0.1))) +
    ggplot2::coord_cartesian(clip = "off")

  list(
    dotplot = p,
    dotplot_height = plot_height,
    full_enrich_results = full_enrich_results
  )
}



# Level 4 internal functions ---------------------------------------------------


#' Prepare Plot Data
#'
#' This function prepares plot data for visualization based on enrichments
#' lists and specified databases.
#'
#' @param enrichments_list A list of enrichments containing data frames for
#'        different databases.
#' @param databases A character vector specifying the databases to be included.
#'
#' @return A list containing two data frames:
#' \describe{
#'   \item{top_plot_data}{A data frame containing the prepared plot data for
#'                          visualization of top combinations.}
#'   \item{full_enrich_results}{A data frame containing the full enrichments
#'                               results.}
#' }
#'
#' @importFrom purrr set_names keep
#' @importFrom dplyr bind_rows starts_with mutate group_by ungroup arrange
#'                   semi_join filter
#' @importFrom tidyr pivot_longer separate_wider_regex replace_na pivot_wider
#'                   unite
#' @importFrom stats na.omit
#' @importFrom rlang .data
#'
prepare_plot_data <- function(
    enrichments_list,
    databases
    ) {

  plot_data <-
    enrichments_list |>
    purrr::set_names(databases) |>
    purrr::keep(is.data.frame) |>
    dplyr::bind_rows(.id = "db") |>
    tidyr::pivot_longer(dplyr::starts_with("Cluster"), names_to = "p_odd") |>
    tidyr::separate_wider_regex(
      .data$p_odd,
      c("Cluster_", cluster = "\\d+", "_", type = ".*"),
      too_few = "align_start"
    ) |>
    tidyr::replace_na(list(type = "adj.p_value")) |>
    tidyr::pivot_wider(names_from = .data$type, values_from = .data$value)

  plot_data <- plot_data |>
    dplyr::group_by(.data$db, .data$BioProcess) |>
    dplyr::mutate(avg_odds_ratio = mean(.data$odds_ratios, na.rm = TRUE)) |>
    dplyr::ungroup() |>
    dplyr::arrange(
      dplyr::desc(.data$avg_odds_ratio),
      .data$db, .data$BioProcess
      )

  # Initialize the cluster counts
  cluster_counts <-
    vector("integer", length = max(plot_data$cluster, na.rm = TRUE))
  names(cluster_counts) <- as.character(seq_along(cluster_counts))

  selected_combos <- list()

  min_threshold <- 5

  # Iterate through combinations
  i <- 1
  while (i <= nrow(plot_data)) {
    combo <- plot_data[i, ]

    # Get the relevant rows from the original data
    # relevant_rows <-
    #   dplyr::filter(
    #     plot_data,
    #     db == combo$db,
    #     BioProcess == combo$BioProcess
    #     )
    relevant_rows <-
      dplyr::filter(
        plot_data,
        plot_data[["db"]] == combo[["db"]],
        plot_data[["BioProcess"]] == combo[["BioProcess"]]
      )

    # Update counts for non-NA odds_ratios clusters
    update_counts <-
      table(relevant_rows$cluster[!is.na(relevant_rows$odds_ratios)])

    # Identify clusters in the current combination
    combo_clusters <- as.numeric(names(update_counts))

    # Check if any cluster in the combination is below the threshold
    include_combo <- any(cluster_counts[combo_clusters] < min_threshold)

    if (include_combo) {
      # Temporarily update cluster counts for evaluation
      temp_cluster_counts <- cluster_counts
      temp_cluster_counts[combo_clusters] <-
        temp_cluster_counts[combo_clusters] + update_counts

      # Final check to ensure we don't exclude necessary clusters under the
      # threshold. This step might seem redundant given the current logic but
      # could be adjusted for more complex conditions
      if (any(temp_cluster_counts <= min_threshold |
              temp_cluster_counts > min_threshold)) {
        # Commit the update if the combination is still eligible
        cluster_counts <- temp_cluster_counts
        selected_combos[[length(selected_combos) + 1]] <- combo
      }
    }

    # Check stopping conditions
    # Stop if all combos have been evaluated or the sum of cluster_counts
    # exceeds 5 times the number of clusters
    if (i == nrow(plot_data) || sum(cluster_counts) > 5 *
        length(cluster_counts)) {
      break
    }

    i <- i + length(cluster_counts)   # Jump to next combo
  }

  # Combine selected combos into a dataframe
  top_combos <- do.call(
    rbind,
    selected_combos
    )

  # Filter the original data to keep only rows matching the top 5 combinations
  top_plot_data <- plot_data |>
    dplyr::semi_join(
      top_combos,
      by = c(
        "db",
        "BioProcess"
        )
      )

  full_enrich_results <- stats::na.omit(plot_data)

  top_plot_data <- top_plot_data |>
    tidyr::unite(
      .data$db,
      .data$BioProcess,
      col = "term",
      sep = ": "
      )

  top_plot_data$term <- factor(
    top_plot_data$term,
    levels = rev(unique(top_plot_data$term))
    )

  list(
    top_plot_data = top_plot_data,
    full_enrich_results = full_enrich_results
    )
}