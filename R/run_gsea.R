# Library import ---------------------------------------------------------------
library(here)
library(clusterProfiler)
library(AnnotationDbi)
library(biomaRt)
library(knitr)
library(tidyverse)
library(viridis)
library(data.table)
library(fs)



# Export functions definitions -------------------------------------------------


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
#' @param params A list specifying the parameters for the enrichment analysis.
#' @param databases A character vector containing strings of the database names
#' to be used for enrichment analysis.
#' @param plot_title An optional string specifying the title of the plot. If not
#' provided, a default title based on the analysis will be used.
#' @param background An optional list of standard gene symbols to be used as the
#' background for the enrichment analysis instead of the background chosen by
#' the `enricher`. The default is an empty list, which implies the use of the
#' default background set by the enrichment tool.
#'
#' @return An object containing the results of the Gene Set Enrichment Analysis,
#' including any plots generated during the analysis.
#'
#' @importFrom clusterProfiler enricher
#'
#' @examples
#' # Example usage
#' pcgsea(your_df, path_to_downloaded_dbs, your_params, "My GSEA Plot")
#'
#' @export
#'
run_gsea <- function(clustered_genes,
                     downloaded_dbs_filepath,
                     params = NA,
                     plot_title = "",
                     background = NULL) {

  control_inputs_run_gsea(clustered_genes,
                          downloaded_dbs_filepath,
                          params,
                          plot_title,
                          background)

  set_default_params(params)

  all_term2genes <- dbs_to_term2genes(downloaded_dbs_filepath)



  ## Prepare objects
  unique_clusters <- sort(unique(clustered_genes[, 2]))

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
    selected_rows <- clustered_genes[clustered_genes$Cluster == cluster, ]

    cluster_genes <- selected_rows$Gene
    gene_list <- as.list(cluster_genes)
    gene_list <- unlist(gene_list)
    gene_list <- toupper(gene_list)

    at_least_one_result = FALSE

    ## Run clusterProfiler and process output
    for (database_name in names(all_term2genes)) {
      term2gene <- all_term2genes[[database_name]]

      enrichment <- 
        clusterProfiler::enricher(gene = gene_list,
                                  pvalueCutoff  = params$adj_p_value,
                                  pAdjustMethod = params$pAdjustMethod,
                                  universe      = background,
                                  minGSSize     = params$minGSSize,
                                  maxGSSize     = params$maxGSSize,
                                  qvalueCutoff  = params$qvalueCutoff,
                                  gson          = NULL,
                                  TERM2GENE     = term2gene,
                                  TERM2NAME     = NA)

      enrichment <- as.data.frame(enrichment)

      if (is.null(enrichment) || (nrow(enrichment) == 0 &
                                  ncol(enrichment) == 0)) {
        next
      }

      at_least_one_result = TRUE

      enrichment_results[[length(enrichment_results) + 1]] <- enrichment

      # Store all for returning in the end.
      name <- sprintf("cluster: %s, database: %s", cluster, database_name)
      raw_results[[name]] <- enrichment
    }

    if (at_least_one_result) {
      all_db_results <- process_enrichment_results(all_db_results,
                                                   enrichment_results,
                                                   params$adj_p_value,
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
    result <- make_enrich_dotplot(all_db_results,
                                  names(all_term2genes),
                                  plot_title)
  } else {
    stop("No database led to an enrichment result!")
  }

  list(dotplot = result[[1]],
       dotplot_nrows = result[[2]],
       full_enrich_results = result[[3]],
       raw_results = raw_results)
}



# Level 1 internal functions ---------------------------------------------------


#' Control Input Parameters for Analysis Functions
#'
#' @description
#' This function validates the input parameters for analysis functions.
#' It checks the
#' structure and data types of `clustered_genes`, `downloaded_dbs_filepath`,
#'  `params`,
#' `plot_title`, and `background`. If any condition is not met, the function
#' stops the
#' script and produces an informative error message.
#'
#' @param clustered_genes A data structure representing clustered genes to be
#' checked by `check_clustered_genes`.
#' @param downloaded_dbs_filepath A character string specifying the path to the
#'  TSV file to be checked by `check_tsv_file`.
#' @param params A list of parameters to be checked by `check_params`.
#' @param plot_title A character string specifying the title for the plot. Must
#'  be of length 1.
#' @param background A character vector specifying background genes or NULL. If
#'  provided, must be a character vector.
#'
#' @return This function does not return a value. It stops with an error message
#'  if any
#'         of the input parameters do not meet the required conditions.
#'
control_inputs_run_gsea <- function(clustered_genes,
                                    downloaded_dbs_filepath,
                                    params,
                                    plot_title,
                                    background) {

  check_clustered_genes(clustered_genes)

  check_tsv_file(downloaded_dbs_filepath)

  check_params(params)

  if (!is.character(plot_title) || length(plot_title) != 1) {
    stop("plot_title must be a character with length 1")
  }

  if (!is.null(background)) {
    if (!is.character(background)) {
      stop("background must be a character vector or NULL.")
    }
  }
}


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
    adj_p_value = 0.05,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    qvalueCutoff = 0.2
  )

  if (is.na(params)) {
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
#' @param file_path The file path to the .tsv file containing the database
#' information.
#' @return A nested list where the first level of names corresponds to database
#' names ('DB'),
#'         the second level to gene sets ('Geneset'), and the innermost lists
#'         contain gene names ('Gene') associated with each gene set.
#'
#' @importFrom readr read_tsv
#'
dbs_to_term2genes <- function(downloaded_dbs_filepath) {

  data <- readr::read_tsv(downloaded_dbs_filepath, col_types = cols())

  db_split <- split(data, data$DB)

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
process_enrichment_results <- function(all_db_results,
                                       enrichment_results,
                                       adjP_threshold,
                                       column_name,
                                       count_column_name,
                                       background = FALSE) {

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
#' 
prepare_plot_data <- function(enrichments_list, databases) {
  plot_data <-
    enrichments_list %>%
    purrr::set_names(databases) %>%
    purrr::keep(is.data.frame) %>%
    dplyr::bind_rows(.id = "db") %>%
    tidyr::pivot_longer(dplyr::starts_with("Cluster"), names_to = "p_odd") %>%
    tidyr::separate_wider_regex(
      p_odd,
      c("Cluster_", cluster = "\\d+", "_", type = ".*"),
      too_few = "align_start"
    ) %>%
    tidyr::replace_na(list(type = "adj.p_value")) %>%
    tidyr::pivot_wider(names_from = type, values_from = value) %>%
    # na.omit() %>%
    {.}

  plot_data <- plot_data %>%
    dplyr::group_by(db, BioProcess) %>%
    dplyr::mutate(avg_odds_ratio = mean(odds_ratios, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(desc(avg_odds_ratio), db, BioProcess)

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
    relevant_rows <-
      dyplr::filter(plot_data, db == combo$db, BioProcess == combo$BioProcess)

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
  top_combos <- do.call(rbind, selected_combos)

  # Filter the original data to keep only rows matching the top 5 combinations
  top_plot_data <- plot_data %>%
    dplyr::semi_join(top_combos, by = c("db", "BioProcess"))

  full_enrich_results <- stats::na.omit(plot_data)

  top_plot_data <- top_plot_data %>%
    tidyr::unite(db, BioProcess, col = "term", sep = ": ")

  top_plot_data$term <- factor(top_plot_data$term,
                               levels = rev(unique(top_plot_data$term)))

  list(top_plot_data = top_plot_data,
       full_enrich_results = full_enrich_results)
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
#' @importFrom viridis scale_color_viridis
#' @importFrom scales oob_squish
#' 
make_enrich_dotplot <- function(enrichments_list,
                                databases,
                                title = "Title") {

  results <- prepare_plot_data(enrichments_list, databases)
  top_plot_data <- results$top_plot_data
  full_enrich_results <- results$full_enrich_results

  # dotplot_nrows <- as.integer(nrow(enrichr_results)/6) + 1
  dotplot_nrows <- 2

  p <- ggplot2::ggplot(top_plot_data, 
                       ggplot2::aes(cluster,
                                    term,
                                    size = -log10(adj.p_value))) +
    ggplot2::geom_point(aes(color = odds_ratios)) +
    ggplot2::ylab("database: term") +
    viridis::scale_color_viridis(
      "odds\nratio",
      option = "inferno",
      direction = -1,
      end = 0.9,
      labels = function(x) round(x, 2),
      guide = ggplot2::guide_colorbar(
        barheight = unit(15, "mm"),
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
    )  +
    ggplot2::coord_fixed() +
    # facet_wrap(vars(db), ncol = 1, scales = "free_y") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = 12)
      #   legend.key.height = unit(3, "mm"),
      #   legend.key.width = unit(3, "mm"),
      #   legend.margin = margin(),
      #   panel.spacing = unit(-.5, "pt"),
    ) +
    ggplot2::labs(title = title)

  list(p, 
       dotplot_nrows, 
       full_enrich_results)
}


#' Get Enrichr Gene Sets
#'
#' Get gene sets from specified Enrichr databases.
#'
#' @param databases A character vector of database names to download from 
#'                  Enrichr.
#'
#' @return A named list of gene sets from the specified Enrichr databases. Each
#'         database is represented as a list, with gene set names as list 
#'         names and vectors of human gene symbols as list elements.
#'
#' @importFrom stats tryCatch message
#'
enrichr_get_genesets <- function(databases) {
  
  setNames(lapply(databases, function(dbx){
    fpath <-
      paste0("http://amp.pharm.mssm.edu/Enrichr/",
             "geneSetLibrary?mode=text&libraryName=",dbx)
    fhandle <- file(fpath)
    dblines <- tryCatch({
      readLines(con = fhandle)
    }, error = function(e){
      message(e, "\nFailed reading database: ", dbx)
      NULL
    })
    close(fhandle)
    if (is.null(dblines)) {
      return(list())
    }else {
      res <- strsplit(dblines, "\t")
      names(res) <- sapply(res, function(x) x[1])
      res <- lapply(res, function(x) x[3:length(x)])
      return(res)
    }
  }), databases)
}



# Level 2 internal functions ---------------------------------------------------


#' Check Clustered Genes List for Required Conditions
#'
#' @description
#' This function checks if a given list `clustered_genes` contains dataframes
#' that meet the specified conditions. Each dataframe in the list must have
#' exactly
#' 2 columns. The first column must contain strings, and the second column must
#' contain integers. If any condition is not met, the function stops the script
#' and produces an informative error message.
#'
#' @param clustered_genes A list of dataframes to be checked for the required
#' conditions.
#'
#' @return This function does not return a value. It stops with an error message
#'         if the conditions are not met.
#'
check_clustered_genes <- function(clustered_genes) {

  # Check if clustered_genes is a list
  if (!is.list(clustered_genes)) {
    stop("clustered_genes must be a list.", call. = FALSE)
  }

  # Iterate over each dataframe in the list
  for (df in clustered_genes) {

    # Check if each element is a dataframe with exactly 2 columns
    if (!is.data.frame(df) || ncol(df) != 2) {
      stop(paste("Each element of clustered_genes must be a dataframe with",
                 "exactly 2 columns."), call. = FALSE)
    }

    # Check if the first column of each dataframe contains strings
    if (!is.character(df[[1]])) {
      stop(paste("The first column of each dataframe in clustered_genes must",
                 "contain strings."), call. = FALSE)
    }

    # Check if the second column of each dataframe contains integers
    if (!is.integer(df[[2]]) && !all(df[[2]] == as.integer(df[[2]]))) {
      stop(paste("The second column of each dataframe in clustered_genes must",
                 "contain integers."), call. = FALSE)
    }
  }
}



#' Check TSV File for Required Conditions
#'
#' @description
#' This function checks if a given TSV file has the required columns: DB,
#' Geneset,
#' Gene, and that all elements in all columns are strings. If any condition
#' is not met,
#' it stops the script and produces an informative error message.
#'
#' @param filepath A character string specifying the path to the TSV file.
#'
#' @return This function does not return a value. It stops with an error message
#'         if the conditions are not met.
#'
#' @importFrom readr read_tsv cols col_character
#'
check_tsv_file <- function(filepath) {

  required_cols <- c("DB", "Geneset", "Gene")

  data <- tryCatch({
    readr::read_tsv(filepath,
                    col_types = readr::cols(.default = readr::col_character()))
  }, error = function(e) {
    stop("Error reading the file: ", filepath, "\n", e$message, call. = FALSE)
  })

  # Check for required columns
  if (!all(required_cols %in% names(data))) {
    stop("The file must contain the following columns: ",
         paste(required_cols, collapse = ", "), call. = FALSE)
  }

  # Check that all elements are strings
  if (!all(sapply(data, is.character))) {
    stop("All elements in all columns must be strings.", call. = FALSE)
  }

  # Check for extra columns
  extra_cols <- setdiff(names(data), required_cols)
  if (length(extra_cols) > 0) {
    stop("The file contains extra columns: ",
         paste(extra_cols, collapse = ", "),
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
check_params <- function(params) {

  required_params <- list(
    float_param1 = "numeric",
    char_param = "character",
    int_param1 = "integer",
    int_param2 = "integer",
    float_param2 = "numeric"
  )

  # Allow params to be NA
  if (is.na(params)) {
    return(invisible(NULL))
  }

  # Check if params is a list
  if (!is.list(params)) {
    stop("The input must be a list.", call. = FALSE)
  }

  # Check for extra elements
  extra_params <- setdiff(names(params), names(required_params))
  if (length(extra_params) > 0) {
    stop("The list contains extra elements: ",
         paste(extra_params, collapse = ", "), call. = FALSE)
  }

  # Check for required elements and their data types
  for (param in names(params)) {
    if (param %in% names(required_params)) {
      expected_type <- required_params[[param]]
      actual_type <- class(params[[param]])
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
      stop("Unexpected element '", param, "' in the list.", call. = FALSE)
    }
  }
}


