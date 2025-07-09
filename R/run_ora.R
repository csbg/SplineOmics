#' run_ora()
#'
#' @description
#' This function generates a overrepresentation analysis report based
#' on clustered hit levels, gene data, and specified databases. It accomplishes 
#' this by using the R package clusterProfiler. As output, you will receive a
#' list of the plot objects it generated, and an HTML report with embedded files
#' containing the enrichment results, and dotplots visualizing the enrichment.
#'
#' @param levels_clustered_hits A list of dataframes that contain the clustered
#' hits of the different levels. When clustering_results is the variable that 
#' collects the output of the SplineOmics::cluster_hits() function, then an 
#' easy way to get this is clustering_results$clustered_hits_levels. Every 
#' element of that list is a dataframe, with the three columns feature, cluster,
#' gene. feature contains the index number of the feature (for example a protein
#' ), cluster is an integer specifying in which cluster this feature was placed,
#' and gene contains the gene name. It is essential that the gene name matches
#' the gene names used in the databases that are used for this enrichment here.
#' @param databases A dataframe with the three columns: DB containing the 
#' database name, Geneset containng the name of the geneset, and Gene, 
#' containing the name of the gene. This dataframe can be obtained by specifying
#' the desired Enrichr databases and downloading them to a for example .tsv file
#' with the help of the SplineOmics::download_enrichr_databases function, and 
#' then loading this .tsv file as a dataframe. In essence, this dataframe then
#' contains all the database info used for the gene set enrichment analysis with
#' clusterProfiler in this function. 
#' @param report_info A list containing information for the report generation,
#' such as omics_data_type and data_description (this is the list used for all
#' report generating functions of this package).
#' @param cluster_hits_report_name Single character string specifying the name 
#' of the cluster_hits() function report, that contains the results that were
#' used for the overprepresentation analysis here. Must be specified, because
#' otherwise, the connection is not documented.
#' @param clusterProfiler_params A list that specifies the parameters for the 
#' clusterProfiler, such as for example: clusterProfiler_params <- list(
#'   pvalueCutoff = 0.05,
#'   pAdjustMethod = "BH",
#'   minGSSize = 10,
#'   maxGSSize = 500,
#'   qvalueCutoff = 0.2
#' )
#' (Those are all the parameters that can be controlled here). The names are 
#' equivalent to the argument names of clusterProfiler, therefore, check out the 
#' documentation of clusterProfiler for their description. When this argument is
#' not specified, it is per default NULL, in which case default parameters for
#' those are selected, which are equivalent to the parameter values shown in the 
#' example definition above.
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
#' @param plot_titles Titles for the enrichment dotplots generated in the HTML
#' report, default is NA.
#' @param universe Enrichment background data, default is NULL. This is a
#' parameter of clusterProfiler, for the documentation, please check the 
#' documentation of the clusterProfiler R package.
#' @param report_dir Directory where the report will be saved, default is
#' `here::here()`.
#'
#' @return A list of all plot objects, generated for the ora report.
#'
#' @importFrom purrr map2 flatten
#' @importFrom here here
#'
#' @export
#'
run_ora <- function(
    levels_clustered_hits,
    databases,
    report_info,
    cluster_hits_report_name,
    clusterProfiler_params = NA,
    mapping_cfg = list(
      method = "none",
      from_species = NULL,
      to_species   = NULL
      ),
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
      match.call()[-1]
    ),
    eval,
    parent.frame()
  )
  
  input_control <- InputControl$new(args)
  input_control$auto_validate()
  
  levels_clustered_hits <- levels_clustered_hits[
    !vapply(
      levels_clustered_hits,
      is.character,
      logical(1)
    )
  ]
  
  # Control the test not covered by the InputControl class
  control_inputs_run_ora(
    levels_clustered_hits = levels_clustered_hits,
    databases = databases,
    params = clusterProfiler_params,
    mapping_cfg = mapping_cfg,
    plot_titles = plot_titles,
    background = universe,
    cluster_hits_report_name = cluster_hits_report_name
  )
  
  ensure_clusterProfiler() # Deals with clusterProfiler installation.

  levels_clustered_hits <- map_gene_symbols(
    levels_clustered_hits = levels_clustered_hits,
    mapping_cfg = mapping_cfg
  )
  
  all_results <- map2(
    levels_clustered_hits,
    names(levels_clustered_hits),
    ~ manage_ora_level(
      .x,
      .y,
      databases,
      clusterProfiler_params,
      universe
    )
  )
  
  processed_results <- map2(
    all_results,
    names(all_results),
    process_result
  )
  
  # Extract the plots, plot sizes, and header info from the processed results
  plots <- purrr::flatten(map(
    processed_results,
    "plot"
  ))
  plots_sizes <- unlist(map(
    processed_results,
    "plot_size"
  ))
  
  insert_after_each <- function(lst, value) {
    result <- vector("list", 2 * length(lst))
    result[seq(1, length(result), by = 2)] <- lst
    result[seq(2, length(result), by = 2)] <- value
    return(result)
  }
  
  plots <- insert_after_each(
    plots,
    "section_break"
  )
  plots_sizes <- insert_after_each(plots_sizes, 999)
  
  level_headers_info <- map(
    processed_results,
    "header_info"
  )
  
  names(level_headers_info) <- names(all_results)
  
  
  # Move them to the report by shipping them inside report_info: Convenience.
  report_info$databases <- unique(databases[["DB"]])
  report_info$clusterProfiler_params <- clusterProfiler_params
  report_info$cluster_hits_report_name <- cluster_hits_report_name
  report_info$levels_clustered_hits <- levels_clustered_hits
  report_info$background_genes <- universe
  
  if (all(vapply(plots, is.character, logical(1)))) {
    message("No results --> Not generating a report and returning NULL.")
    return(NULL)
  }
  
  print_info_message(
    message_prefix = "ORA analysis",
    report_dir = report_dir
  )
  
  generate_report_html(
    plots = plots,
    plots_sizes = plots_sizes,
    report_info = report_info,
    level_headers_info = level_headers_info,
    report_type = "run_ora",
    filename = "run_ora_report",
    report_dir = report_dir
  )
  
  return(Filter(function(x) !is.character(x), plots))
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
#' @param levels_clustered_hits A list containing clustered hits at various
#' levels.
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
#' @param plot_titles A character vector of titles for the plots, with length
#' matching `levels_clustered_hits`.
#' @param background A character vector of background genes or NULL.
#' @param cluster_hits_report_name Single character string specifying the name 
#' of the cluster_hits() function report, that contains the results that were
#' used for the overprepresentation analysis here. Must be specified, because
#' otherwise, the connection is not documented.
#'
control_inputs_run_ora <- function(
    levels_clustered_hits,
    databases,
    params,
    mapping_cfg,
    plot_titles,
    background,
    cluster_hits_report_name
) {
  
  check_clustered_hits(levels_clustered_hits)
  
  check_databases(databases)
  
  check_params(params)
  
  check_mapping_cfg(mapping_cfg)
  
  if (!is.na(plot_titles)) {
    if (!is.character(plot_titles) ||
        length(plot_titles) != length(levels_clustered_hits)) {
      stop(paste(
        "plot_titles must be a character vector with length == length",
        "length levels_clustered_hits"
      ), call. = FALSE)
    }
  }
  
  
  if (!is.null(background)) {
    if (!is.character(background)) {
      stop_call_false("background must be a character vector or NULL.")
    } else {
      check_genes(background)
    }
  }
  
  if (!is.character(cluster_hits_report_name) 
      || length(cluster_hits_report_name) != 1) {
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
#'
ensure_clusterProfiler <- function() {
  
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop_call_false(
      "The 'clusterProfiler' package is not installed.\n",
      "Please install it manually using the", 
      "command below and re-run the function:\n\n",
      "  BiocManager::install('clusterProfiler')\n\n",
      "This is an optional dependency of the SplineOmics package, ",
      "only needed for optional functionality and not part of the core", 
      "package, ",
      "which is why it must be installed manually if this function is used."
    )
  }
}


#' Manage ORA Analysis for a Specific Level
#' 
#' @noRd
#'
#' @description
#' This function manages the ora analysis for a specific level. It extracts
#' genes associated with the clustered hits, removes rows with `NA` values,
#' and runs the ora analysis using the `run_ora` function.
#'
#' @param clustered_hits A dataframe containing the clustered hits for a
#' specific level. It must include a column named `feature` to extract genes.
#' @param level_name A character string representing the name of the level.
#' @param databases A list of databases for the gene set enrichment analysis.
#' @param clusterProfiler_params Additional parameters for the ora analysis,
#'                               default is NA. Those include adj_p_value,
#'                               pAdjustMethod, etc (see clusterProfiler
#'                               documentation).
#' @param universe Enrichment background data, default is NULL.
#'
#' @return The result of the `run_ora` function, which typically
#'         includes various plots and enrichment results.
#'
manage_ora_level <- function(
    clustered_hits,
    level_name,
    databases,
    clusterProfiler_params,
    universe
) {
  
  clustered_hits <- na.omit(clustered_hits)
  
  message(paste(
    "\n\n\n Running clusterProfiler for the level:",
    level_name
  ))
  
  result <- run_ora_level(
    clustered_genes = clustered_hits,
    databases = databases,
    params = clusterProfiler_params,
    plot_title = level_name,
    universe = universe
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
      ora_results = level_result$ora_results
    )
  }
  
  return(result)
}


#' Build ORA Report
#' 
#' @noRd
#'
#' @description
#' Generates an HTML report for Gene Set Enrichment Analysis (ora) based on
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
build_run_ora_report <- function(
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
    if (any(is(plots[[index]], "character"))) {
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


#' Map gene symbols across species
#' 
#' @noRd
#' 
#' @description
#' This function maps gene symbols from one species to another using optional 
#' external tools. 
#' It is primarily intended to harmonize gene symbols across species in 
#' preparation for downstream analyses 
#' such as overrepresentation analysis (ORA), where gene identifiers must match
#'  those used by the reference database.
#'
#' The function takes a list of data frames (`levels_clustered_hits`), each of
#'  which must contain a column named `gene`.
#' It replaces the gene symbols in this column based on orthology mappings 
#' between the specified source and target species.
#' Mappings are performed via either the `gprofiler2` or `orthogene` package,
#'  depending on user configuration. 
#' If no mapping is requested (`method = "none"`), the gene symbols are returned
#'  unchanged.
#'
#' Only 1:1 orthologs are retained during mapping; genes with no matching
#'  orthologs or ambiguous mappings are left unchanged.
#' The structure and order of the input list and data frames are preserved.
#'
#' @param levels_clustered_hits  List of data frames. Each must have a `gene`
#'        column. Order is preserved.
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
#' @return Same `levels_clustered_hits` list, same order, `gene` column
#'  replaced by mapped symbols.
#' 
map_gene_symbols <- function(
    levels_clustered_hits,
    mapping_cfg
    ) {
  browser()
  stopifnot(
    is.list(levels_clustered_hits),
    all(vapply(levels_clustered_hits,
               \(x) "gene" %in% names(x), logical(1))),
    is.list(mapping_cfg),
    all(c("method","from_species","to_species") %in%
          names(mapping_cfg))
    )
  
  method <- tolower(mapping_cfg$method)
  from   <- mapping_cfg$from_species
  to     <- mapping_cfg$to_species
  
  all_genes <- unlist(lapply(levels_clustered_hits, `[[`, "gene"),
                      use.names = FALSE)
  
  if (any(grepl("_[0-9]+$", all_genes))) {
    warning(
      "Some gene IDs end in '_<digits>'. No automatic stripping ",
      "will be performed.  Clean them yourself if they fail to map.",
      call. = FALSE
      )
  }
  
  unique_genes <- unique(all_genes)
  
  map_vec <- switch(
    method,
    "none" = setNames(
      unique_genes,
      unique_genes
      ),
                    
    "gprofiler" = {
      if (!requireNamespace("gprofiler2", quietly = TRUE))
        stop_call_false(
          "`gprofiler2` not installed; install it or set method = 'none'."
          )

      gp <- gprofiler2::gorth(
        unique_genes,
        source_organism  = from,
        target_organism  = to,
        mthreshold       = Inf,
        filter_na        = TRUE
        )
      
      # gorth names: 'input' = original, 'ortholog_name' = mapped
      setNames(
        gp$ortholog_name,
        gp$input
        )
    },
    
    "orthogene" = {
      if (!requireNamespace("orthogene", quietly = TRUE))
        stop_call_false(
          "`orthogene` not installed; install it or set method = 'none'."
          )
      
      or <- orthogene::convert_orthologs(
        gene_df        = unique_genes,
        input_species  = from,                       # correct arg names
        output_species = to,
        method         = "gprofiler",
        gene_output    = "dict")                     # returns a named list
      
      map_chr <- unlist(    # make it a character vec
        or,
        use.names = TRUE
        )              
      setNames(
        map_chr,
        names(map_chr)
        )
    },
    
    stop_call_false(
      "Unknown mapping method: ",
      method
      )
  )
  
  ## fallback for unmapped -> keep original symbol
  fallback_map <- setNames(   # start with identity map
    unique_genes,
    unique_genes
    )  
  fallback_map[names(map_vec)] <- map_vec       # overwrite mapped ones
  map_vec <- fallback_map                       # use this from here on
  
  # Write back in the same order
  rebuild <- \(df, idx) { df$gene <- map_vec[idx]; df }
  row_counts <- vapply(
    levels_clustered_hits,
    nrow,
    integer(1)
    )
  gene_split <- split(
    all_genes,
    rep(
      seq_along(row_counts),
      row_counts
      )
    )
  Map(
    rebuild,
    levels_clustered_hits,
    gene_split
    )
}


# Level 2 internal functions ---------------------------------------------------


#' Check Clustered Genes Dataframe for Required Conditions
#' 
#' @noRd
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
      stop(paste(
        "Element", i, "of levels_clustered_hits is not a dataframe",
        "but must be one."
      ), call. = FALSE)
    }
    
    # Check if the dataframe contains the columns 'gene' and 'cluster'
    required_columns <- c("feature", "cluster")
    if (!all(required_columns %in% colnames(clustered_hits))) {
      stop(
        paste(
          "The dataframe must contain the columns 'feature' and",
          "'cluster'."
        ),
        call. = FALSE
      )
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
#' @noRd
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
    max_index_overall = NA) {
  if (!is.na(max_index_overall)) {
    if (length(genes) < max_index_overall) {
      stop(
        paste(
          "genes must at least have over", max_index_overall,
          "elements"
        ),
        call. = FALSE
      )
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
    stop(
      sprintf(
        paste(
          "Invalid gene found at index %d: '%s'.",
          "There are %d more invalid elements."
        ),
        first_invalid_index, first_invalid_value, num_invalid
      ),
      call. = FALSE
    )
  }
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
    stop_call_false("All columns in the dataframe must be of type character.")
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
check_mapping_cfg <- function(mapping_cfg) {
  
  # structural checks
  if (!is.list(mapping_cfg) || inherits(mapping_cfg, "data.frame")) {
    stop("`mapping_cfg` must be a *list*, e.g.\n",
         "  list(method = 'gprofiler', from_species = 'cgchok1gshd', ",
         "to_species = 'hsapiens').")
  }
  
  expected <- c("method", "from_species", "to_species")
  missing  <- setdiff(expected, names(mapping_cfg))
  if (length(missing)) {
    stop("`mapping_cfg` is missing field(s): ",
         paste(missing, collapse = ", "), ".")
  }
  
  extra <- setdiff(names(mapping_cfg), expected)
  if (length(extra)) {
    stop("`mapping_cfg` contains unknown field(s): ",
         paste(extra, collapse = ", "), ".\n",
         "Only the fields ", paste(expected, collapse = ", "),
         " are allowed.")
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
      if (!is.character(val) || length(val) != 1 || !nzchar(val))
        stop("`", sp, "` must be a non-empty character scalar.")
    }
  }
  
  # tool-specific checks -
  if (method == "gprofiler") {
    if (!requireNamespace("gprofiler2", quietly = TRUE))
      stop_call_false(
        "`gprofiler2` is not installed. Install with:  
        install.packages('gprofiler2')"
        )
    
    # Skip organism code checking entirely
    message(
      "Note: `from_species` and `to_species` are not checked for validity.\n",
      "Make sure they match g:Profiler's expected organism codes. See:\n",
      "https://biit.cs.ut.ee/gprofiler/page/organism-list"
      )
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


#' Perform Gene Set Enrichment Analysis and plot it.
#' 
#' @noRd
#'
#' @description
#' This function conducts a Gene Set Enrichment Analysis (ora) using either the
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
run_ora_level <- function(
    clustered_genes,
    databases,
    params = NA,
    plot_title = "",
    universe = NULL
) {
  
  params <- set_default_params(params)
  gene_set_collections <- dbs_to_term2genes(databases)
  
  unique_clusters <- sort(unique(clustered_genes$cluster))
  
  ora_results <- list()   # stores all the results from clusterProfiler
  
  for (cluster in unique_clusters) {
    cluster_label <- paste0("cluster_", cluster)
    ora_results[[cluster_label]] <- list()
    
    foreground_genes <- 
      as.character(clustered_genes$gene[clustered_genes$cluster == cluster])
    
    at_least_one_result <- FALSE
    
    # Run clusterProfiler and process output
    for (gene_set_name in names(gene_set_collections)) {
      gene_set_map <- gene_set_collections[[gene_set_name]]
      
      # There can be a problem with a database, and that potential message is 
      # not coming from the SplineOmics code. That is why here always the 
      # database is send as a message, so that the problems can be connected
      # to a specific database.
      message(paste(    
        "\nDatabase:",
        gene_set_name
      ))
      
      ora_result <-
        clusterProfiler::enricher(
          gene = foreground_genes,
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
      
      ora_result_df <- as.data.frame(ora_result)
      if (nrow(ora_result_df) == 0) {
        ora_result_df <- NA  # Mark explicitly
      }

      ora_results[[cluster_label]][[gene_set_name]] <- ora_result_df
    }
  }
  
  ora_results <- add_odds_ratios_to_ora(ora_results)
  
  # Check if any enrichment results exist
  any_result <- any(sapply(ora_results, function(cluster_entry) {
    any(sapply(cluster_entry, function(res) {
      is.data.frame(res) && nrow(res) > 0
    }))
  }))
  
  if (any_result) {
    result <- make_enrich_dotplot(
      ora_results,
      plot_title
    )
  } else {
    message("No cluster led to an enrichment result!")
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
    toc_style
) {
  
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


#' Convert Database File to TERM2GENE List
#' 
#' @noRd
#'
#' @description
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
  
  for (cluster in names(ora_results)) {
    for (gene_set_name in names(ora_results[[cluster]])) {
      df <- ora_results[[cluster]][[gene_set_name]]
      
      if (
        is.data.frame(df) &&
        "GeneRatio" %in% colnames(df) &&
        "BgRatio" %in% colnames(df)
      ) {
        # Parse GeneRatio
        gene_ratio_parts <- strsplit(df$GeneRatio, "/")
        gene_numer <- as.numeric(sapply(gene_ratio_parts, function(x) x[1]))
        gene_denom <- as.numeric(sapply(gene_ratio_parts, function(x) x[2]))
        gene_ratios <- gene_numer / gene_denom
        
        # Parse BgRatio
        bg_ratio_parts <- strsplit(df$BgRatio, "/")
        bg_numer <- as.numeric(sapply(bg_ratio_parts, function(x) x[1]))
        bg_denom <- as.numeric(sapply(bg_ratio_parts, function(x) x[2]))
        bg_ratios <- bg_numer / bg_denom
        
        # Compute odds ratio
        odds_ratios <- gene_ratios / bg_ratios
        df$odds_ratio <- odds_ratios
        
        ora_results[[cluster]][[gene_set_name]] <- df
      }
    }
  }
  
  return(ora_results)
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
    title = "Title"
) {
  
  top_plot_data <- prepare_plot_data(ora_results)
  
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
  colnames(top_df) <- c("cluster", "term", "adj.p_value", "odds_ratio")
  
  top_df
}