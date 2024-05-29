#' # Library import ---------------------------------------------------------------
#' library(here)
#' library(clusterProfiler)
#' library(AnnotationDbi)
#' library(biomaRt)
#' library(knitr)
#' library(tidyverse)
#' library(viridis)
#' library(data.table)
#' library(fs)
#' 
#' 
#' 
#' # Internal functions definitions -----------------------------------------------
#' 
#' 
#' process_enrichment_results <- function(all_db_results, enrichment_results, 
#'                                        adjP_threshold, column_name, 
#'                                        count_column_name, background = FALSE) {
#'   
#'   column_indices <- list(2, 6, 3, 4, 9)
#'   
#'   # Process results for all databases.
#'   for (i in seq_along(enrichment_results)) {
#'     df <- enrichment_results[[i]]
#'     df <- subset(df, df[[column_indices[[2]]]] < adjP_threshold)
#'     
#'     if (nrow(df) == 0) {
#'       next
#'     }
#'     
#'     term_list <- as.list(df[[column_indices[[1]]]])
#'     adjP_list <- as.list(df[[column_indices[[2]]]])
#'     odds_ratio <- as.list(df[[column_indices[[3]]]])
#'     
#'     odds_ratio <- sapply(odds_ratio, function(x) eval(parse(text = x)))
#'     bg_ratio <- as.list(df[[column_indices[[4]]]])
#'     bg_ratio <- sapply(bg_ratio, function(x) eval(parse(text = x)))
#'     odds_ratio <- mapply("/", odds_ratio, bg_ratio)
#'     
#'     gene_count_list <- as.list(df[[column_indices[[5]]]])
#'     
#'     named_list <- list()
#'     
#'     # Loop through the terms
#'     for(j in seq_along(term_list)) {
#'       # Create a sublist for each term with its adjP value and gene count
#'       
#'       # Skip terms that are just supported by one gene.
#'       if (gene_count_list[j] < 2) {
#'         next
#'       }
#'       
#'       sublist <- list(adjP_list[j], odds_ratio[j])
#'       
#'       # Assign this sublist to named_list with the term as its name
#'       term <- as.character(term_list[j])
#'       named_list[[term]] <- sublist
#'     }
#'     
#'     for (name in names(named_list)) {
#'       if (!name %in% all_db_results[[i]]$BioProcess) {
#'         # Create a list with the same structure as your DataFrame
#'         row_index <- nrow(all_db_results[[i]]) + 1
#'         all_db_results[[i]][row_index, "BioProcess"] <- name
#'         all_db_results[[i]][row_index, column_name] <- named_list[[name]][[1]]
#'         all_db_results[[i]][row_index, count_column_name] <- 
#'           named_list[[name]][[2]]
#'       } else {
#'         row_index <- which(all_db_results[[i]]$BioProcess == name)
#'         all_db_results[[i]][row_index, column_name] <- named_list[[name]][[1]]
#'         all_db_results[[i]][row_index, count_column_name] <- 
#'           named_list[[name]][[2]]
#'       }
#'     }
#'   }
#'   return(all_db_results)
#' }
#' 
#' 
#' prepare_plot_data <- function(enrichments_list, databases) {
#'   plot_data <- 
#'     enrichments_list %>%
#'     set_names(databases) %>% 
#'     keep(is.data.frame) %>%
#'     bind_rows(.id = "db") %>%
#'     pivot_longer(starts_with("Cluster"), names_to = "p_odd") %>%
#'     separate_wider_regex(
#'       p_odd,
#'       c("Cluster_", cluster = "\\d+", "_", type = ".*"),
#'       too_few = "align_start"
#'     ) %>% 
#'     replace_na(list(type = "adj.p_value")) %>% 
#'     pivot_wider(names_from = type, values_from = value) %>% 
#'     # na.omit() %>% 
#'     {.}
#' 
#'   plot_data <- plot_data %>%
#'     group_by(db, BioProcess) %>%
#'     mutate(avg_odds_ratio = mean(odds_ratios, na.rm = TRUE)) %>%
#'     ungroup() %>%
#'     arrange(desc(avg_odds_ratio), db, BioProcess)
#'   
#'   # Initialize the cluster counts
#'   cluster_counts <- 
#'     vector("integer", length = max(plot_data$cluster, na.rm = TRUE))
#'   names(cluster_counts) <- as.character(seq_along(cluster_counts))
#'   
#'   selected_combos <- list()
#'   
#'   min_threshold <- 5
#'   
#'   # Iterate through combinations
#'   i <- 1
#'   while (i <= nrow(plot_data)) {
#'     combo <- plot_data[i, ]
#'     
#'     # Get the relevant rows from the original data
#'     relevant_rows <- 
#'       filter(plot_data, db == combo$db, BioProcess == combo$BioProcess)
#'     
#'     # Update counts for non-NA odds_ratios clusters
#'     update_counts <- 
#'       table(relevant_rows$cluster[!is.na(relevant_rows$odds_ratios)])
#'     
#'     # Identify clusters in the current combination
#'     combo_clusters <- as.numeric(names(update_counts))
#'     
#'     # Check if any cluster in the combination is below the threshold
#'     include_combo <- any(cluster_counts[combo_clusters] < min_threshold)
#'     
#'     if (include_combo) {
#'       # Temporarily update cluster counts for evaluation
#'       temp_cluster_counts <- cluster_counts
#'       temp_cluster_counts[combo_clusters] <- 
#'         temp_cluster_counts[combo_clusters] + update_counts
#'       
#'       # Final check to ensure we don't exclude necessary clusters under the 
#'       # threshold. This step might seem redundant given the current logic but 
#'       # could be adjusted for more complex conditions
#'       if (any(temp_cluster_counts <= min_threshold | 
#'               temp_cluster_counts > min_threshold)) {
#'         # Commit the update if the combination is still eligible
#'         cluster_counts <- temp_cluster_counts
#'         selected_combos[[length(selected_combos) + 1]] <- combo
#'       }
#'     }
#'     
#'     # Check stopping conditions
#'     # Stop if all combos have been evaluated or the sum of cluster_counts 
#'     # exceeds 5 times the number of clusters
#'     if (i == nrow(plot_data) || sum(cluster_counts) > 5 * 
#'         length(cluster_counts)) {
#'       break
#'     }
#'     
#'     i <- i + length(cluster_counts)   # Jump to next combo
#'   }
#'   
#'   # Combine selected combos into a dataframe
#'   top_combos <- do.call(rbind, selected_combos)
#'   
#'   # Filter the original data to keep only rows matching the top 5 combinations
#'   top_plot_data <- plot_data %>%
#'     dplyr::semi_join(top_combos, by = c("db", "BioProcess"))
#'   
#'   full_enrich_results <- na.omit(plot_data)
#' 
#'   top_plot_data <- top_plot_data %>% 
#'     unite(db, BioProcess, col = "term", sep = ": ") 
#'   
#'   top_plot_data$term <- factor(top_plot_data$term, levels = 
#'                                  rev(unique(top_plot_data$term)))
#'   
#'   return(list(top_plot_data = top_plot_data, 
#'          full_enrich_results = full_enrich_results))
#' }
#' 
#' 
#' make_enrich_dotplot <- function(enrichments_list, databases, title = "Title") {
#'   
#'   results <- prepare_plot_data(enrichments_list, databases)
#'   top_plot_data <- results$top_plot_data
#'   full_enrich_results <- results$full_enrich_results
#'   
#'   # dotplot_nrows <- as.integer(nrow(enrichr_results)/6) + 1
#'   dotplot_nrows <- 2
#'   
#'   p <- ggplot(top_plot_data, aes(cluster, term, size = -log10(adj.p_value))) +
#'     geom_point(aes(color = odds_ratios)) + 
#'     ylab("database: term") +
#'     scale_color_viridis(
#'       "odds\nratio",
#'       option = "inferno",
#'       direction = -1,
#'       end = 0.9,
#'       labels = function(x) round(x, 2),
#'       guide = guide_colorbar(
#'         barheight = unit(15, "mm"),
#'         barwidth = unit(2, "mm"),
#'         ticks = FALSE
#'       )
#'     ) +
#'     scale_size_area(
#'       max_size = 3,
#'       limits = c(0, 2),
#'       breaks = c(0, 1, 2),
#'       labels = c("0", "1", "2 or higher"),
#'       oob = scales::oob_squish
#'     )  +
#'     coord_fixed() +
#'     # facet_wrap(vars(db), ncol = 1, scales = "free_y") +
#'     theme_bw() +
#'     theme(
#'       panel.grid = element_blank(),
#'       plot.title = element_text(size = 12)
#'       #   legend.key.height = unit(3, "mm"),
#'       #   legend.key.width = unit(3, "mm"),
#'       #   legend.margin = margin(),
#'       #   panel.spacing = unit(-.5, "pt"),
#'     ) +
#'     labs(title = title)
#'   
#'   result <- list(p, dotplot_nrows, full_enrich_results)
#' }
#' 
#' 
#' # Control if the function input is valid.
#' control_input <- function(clustered_genes, all_term2genes, params,  
#'                           plot_title, background) {
#'   # Check if `clustered_genes` is a dataframe with exactly 2 columns
#'   if (!is.data.frame(clustered_genes) || ncol(clustered_genes) != 2) {
#'     stop("clustered_genes must be a dataframe with exactly 2 columns.")
#'   }
#'   
#'   # Check if the first column of `clustered_genes` contains strings
#'   if (!is.character(clustered_genes[[1]])) {
#'     stop("The first column of clustered_genes must contain strings.")
#'   }
#'   
#'   # Check if the second column of `clustered_genes` contains integers
#'   if (!is.integer(clustered_genes[[2]]) && 
#'       !all(clustered_genes[[2]] == as.integer(clustered_genes[[2]]))) {
#'     stop("The second column of clustered_genes must contain integers.")
#'   }
#' 
#'   if (!is.list(params)) {
#'     stop("params must be a list.")
#'   }
#'   
#'   # Check if `all_term2genes` is a named list of dataframes
#'   if (!is.list(all_term2genes) || is.null(names(all_term2genes))) {
#'     stop("The input must be a named list.")
#'   } else {
#'     for (name in names(all_term2genes)) {
#'       if (!is.data.frame(all_term2genes[[name]])) {
#'         stop(paste("The list element", sQuote(name), "is not a dataframe."))
#'       }
#'     }
#'   }
#'   
#'   if (!is.character(plot_title) || length(plot_title) != 1) {
#'     stop("plot_title must be a single string.")
#'   }
#'   
#'   if (!is.null(background)) {
#'     if (!is.character(background)) {
#'       stop("background must be a list of strings or NULL.")
#'     }
#'   }
#' }
#' 
#' 
#' enrichr_get_genesets <- function(databases){
#'   setNames(lapply(databases, function(dbx){
#'     fpath <- 
#'       paste0("http://amp.pharm.mssm.edu/Enrichr/",
#'              "geneSetLibrary?mode=text&libraryName=",dbx)
#'     fhandle <- file(fpath)
#'     dblines <- tryCatch({
#'       readLines(con = fhandle)
#'     }, error = function(e){
#'       message(e, "\nFailed reading database: ", dbx)
#'       NULL
#'     })
#'     close(fhandle)
#'     if (is.null(dblines)) {
#'       return(list())
#'     }else {
#'       res <- strsplit(dblines, "\t")
#'       names(res) <- sapply(res, function(x) x[1])
#'       res <- lapply(res, function(x) x[3:length(x)])
#'       return(res)
#'     }
#'   }), databases)
#' }
#' 
#' 
#' 
#' # Export functions definitions -------------------------------------------------
#' 
#' 
#' #' Download Gene Sets from Enrichr and Save as TSV
#' #'
#' #' Retrieves gene sets from specified Enrichr libraries, formats the data,
#' #' and saves it to a TSV file in a specified directory. Each row in the output
#' #' file represents a gene in a gene set from a specific database.
#' #'
#' #' @param gene_set_lib A character vector of Enrichr libraries to download gene
#' #'        sets from, e.g., c("KEGG_2021_Human").
#' #' @param dirname Optional. The directory name where the gene sets will be 
#' #'        saved. Defaults to "databases".
#' #'
#' #' @return Invisibly returns the path to the saved TSV file containing all
#' #'         downloaded gene sets. The file is named "all_databases_" followed by
#' #'         a timestamp, saved within the specified directory. Each row contains
#' #'         a gene set name, the associated gene, and the source database, with
#' #'         columns "DB", "Geneset", and "Gene".
#' #'
#' #' @examples
#' #' download_databases(gene_set_lib = c("KEGG_2021_Human"), dirname = 
#' #'                    "my_gene_sets")
#' #'
#' #' @importFrom data.table data.table
#' #' @importFrom fs dir_create
#' #' @importFrom readr write_tsv
#' #' @importFrom here here
#' #' @export
#' #'
#' #' @details Leverages `enrichr_get_genesets` to fetch gene set data from 
#' #'          Enrichr. Processes the data into a long-format data.table, each row 
#' #'          corresponding to a gene in a gene set, along with metadata 
#' #'          indicating the gene set and source database. The data is then saved 
#' #'          as a TSV file for downstream analyses. An internet connection is 
#' #'          required to fetch data from Enrichr.
#' #'          
#' download_databases <- function(gene_set_lib, 
#'                                dirname = "databases") {
#'   
#'   genesets <- enrichr_get_genesets(databases = gene_set_lib)
#'   
#'   genesets <- do.call(rbind, lapply(names(genesets), function(db.nam){
#'     do.call(rbind,lapply(names(genesets[[db.nam]]), function(set.nam){
#'       data.table(DB = db.nam, Geneset = set.nam, Gene = 
#'                    genesets[[db.nam]][[set.nam]])
#'     }))
#'   }))
#'   
#'   genesets[,Gene := gsub(",.+$", "", Gene)]
#'   
#'   fs::dir_create(dirname)
#'   
#'   timestamp <- format(Sys.time(), "%d_%m_%Y-%H_%M_%S")
#'   filename <- paste0("all_databases_", timestamp, ".tsv")
#'   write_tsv(x = genesets, here::here(dirname, filename))
#' }
#' 
#' 
#' #' Convert Database File to TERM2GENE List
#' #'
#' #' Reads a specified .tsv file containing information about databases, 
#' #' gene sets, and genes. The file should have three columns: 'DB' for database 
#' #' names, Geneset' for gene set identifiers, and 'Gene' for gene names. This 
#' #' function organizes this information into a nested list. Each top-level 
#' #' element corresponds to a unique database, and within each, gene sets map to 
#' #' lists of associated genes.
#' #'
#' #' @param file_path The file path to the .tsv file containing the database 
#' #' information.
#' #' @return A nested list where the first level of names corresponds to database 
#' #' names ('DB'), 
#' #'         the second level to gene sets ('Geneset'), and the innermost lists 
#' #'         contain gene names ('Gene') associated with each gene set.
#' #' @examples
#' #' # Assuming 'example.tsv' is in the working directory and formatted correctly:
#' #' list_structure <- dbs_to_term2genes("example.tsv")
#' #' @importFrom readr read_tsv
#' #' @importFrom stats setNames
#' #' @export
#' #'
#' #' @details The function assumes the input file is a tab-separated values (TSV) 
#' #'          file with at least three columns: 'DB', 'Geneset', and 'Gene'. It is 
#' #'          important that the file follows this format for the function to work 
#' #'          as expected. The resulting nested list structure allows for easy 
#' #'          access to genes associated with each gene set within each database, 
#' #'          facilitating further analysis or processing.
#' #'          
#' dbs_to_term2genes <- function(file_path) {
#'   
#'   library(readr) # Ensure readr is loaded for read_tsv
#'   
#'   # Read the TSV file into a dataframe
#'   data <- read_tsv(file_path, col_types = cols())
#'   
#'   # Split the data by the 'DB' column to create a list of dataframes, one for 
#'   # each database
#'   db_split <- split(data, data$DB)
#'   
#'   # Transform each database's data into the desired long dataframe format with 
#'   # renamed columns
#'   all_TERM2GENEs <- lapply(db_split, function(db_df) {
#'     # Directly return the dataframe format with renamed columns
#'     df_with_renamed_columns <- db_df[, c("Geneset", "Gene")]
#'     colnames(df_with_renamed_columns) <- c("term", "gene")
#'     return(df_with_renamed_columns)
#'   })
#'   
#'   # Set the names of the list to be the database names
#'   names(all_TERM2GENEs) <- names(db_split)
#'   
#'   return(all_TERM2GENEs)
#' }
#' 
#' 
#' 
#' #' Perform Gene Set Enrichment Analysis and plot it.
#' #'
#' #' This function conducts a Gene Set Enrichment Analysis (GSEA) using either the 
#' #' clusterProfiler package. Afterwards, it plots the results.
#' #' It allows for customization of enrichment parameters, selection of databases, 
#' #' and optionally specifying a custom plot title and background gene list.
#' #'
#' #' @param clustered_genes A dataframe with two columns: the first column 
#' #' contains the standard gene symbol, and the second column contains an integer 
#' #' specifying the cluster.
#' #' @param params A list specifying the parameters for the enrichment analysis. 
#' #' @param databases A character vector containing strings of the database names 
#' #' to be used for enrichment analysis.
#' #' @param plot_title An optional string specifying the title of the plot. If not 
#' #' provided, a default title based on the analysis will be used.
#' #' @param background An optional list of standard gene symbols to be used as the 
#' #' background for the enrichment analysis instead of the background chosen by 
#' #' the `enricher`. The default is an empty list, which implies the use of the 
#' #' default background set by the enrichment tool.
#' #'
#' #' @return An object containing the results of the Gene Set Enrichment Analysis, 
#' #' including any plots generated during the analysis. 
#' #'
#' #' @examples
#' #' # Example usage
#' #' pcgsea(your_df, your_params, your_databases, "My GSEA Plot")
#' #'
#' #' @export
#' pcgsea <- function(clustered_genes, all_term2genes, params = NULL,
#'                    plot_title = "", background = NULL) {
#'   
#'   if (is.null(params)) {
#'     params <- list(adj_p_value = 0.05, pAdjustMethod = "BH",
#'                    minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2)
#'   }
#'   
#'   control_input(clustered_genes, all_term2genes, params, plot_title, background)
#'   
#'   ## Prepare objects -----------------------------------------------------------
#'   unique_clusters <- sort(unique(clustered_genes[, 2]))
#'   
#'   all_db_results <- vector("list", length(all_term2genes))
#'   
#'   for (i in seq_along(all_db_results)) {
#'     all_db_results[[i]] <- data.frame(BioProcess=character())
#'   }
#'   
#'   raw_results <- list()
#'   
#'   for (cluster in unique_clusters) {
#'     enrichment_results <- list()
#'     
#'     for (i in seq_along(all_db_results)) {
#'       column_name <- paste("Cluster", cluster, sep = "_") 
#'       count_column_name <- paste0(column_name, "_odds_ratios")
#'       
#'       # Initialize the column as an empty vector for the df if it has no rows
#'       if (nrow(all_db_results[[i]]) == 0) {
#'         all_db_results[[i]][[column_name]] <- vector("logical", length = 0) 
#'         all_db_results[[i]][[count_column_name]] <- 
#'           vector("logical", length = 0)
#'       } else {
#'         all_db_results[[i]][[column_name]] <- NA
#'         all_db_results[[i]][[count_column_name]] <- NA
#'       }
#'     }
#'     
#'     # Select rows where 'Cluster' equals the current cluster number
#'     selected_rows <- clustered_genes[clustered_genes$Cluster == cluster, ]
#'     
#'     cluster_genes <- selected_rows$Gene
#'     gene_list <- as.list(cluster_genes)
#'     gene_list <- unlist(gene_list)
#'     gene_list <- toupper(gene_list)
#'     
#'     at_least_one_result = FALSE
#'     
#'     ## Run clusterProfiler and process output ----------------------------------
#'     for (database_name in names(all_term2genes)) {
#'       term2gene <- all_term2genes[[database_name]]
#'       
#'       enrichment <- enricher(gene          = gene_list,
#'                              pvalueCutoff  = params$adj_p_value,
#'                              pAdjustMethod = params$pAdjustMethod,
#'                              universe      = background,
#'                              minGSSize     = params$minGSSize,
#'                              maxGSSize     = params$maxGSSize,
#'                              qvalueCutoff  = params$qvalueCutoff,
#'                              gson          = NULL,
#'                              TERM2GENE     = term2gene, 
#'                              TERM2NAME     = NA)
#' 
#'       enrichment <- as.data.frame(enrichment)  
#'       
#'       if (is.null(enrichment) || (nrow(enrichment) == 0 & 
#'                                   ncol(enrichment) == 0)) {
#'         next
#'       }
#'       
#'       at_least_one_result = TRUE
#'       
#'       enrichment_results[[length(enrichment_results) + 1]] <- enrichment
#'       
#'       # Store all for returning in the end.
#'       name <- sprintf("cluster: %s, database: %s", cluster, database_name)
#'       raw_results[[name]] <- enrichment
#'     }
#'     
#'     if (at_least_one_result) {
#'       all_db_results <- process_enrichment_results(all_db_results,
#'                                                    enrichment_results,
#'                                                    params$adj_p_value,
#'                                                    column_name,
#'                                                    count_column_name,
#'                                                    background = use_background
#'       )
#'     } else {
#'       print(paste0("No enrichment results for cluster", cluster))
#'     }
#'   }
#'   
#'   # Make dotplot ---------------------------------------------------------------
#'   any_result <- sapply(all_db_results, function(df) nrow(df) > 0)
#'   has_true <- any(any_result)
#'   
#'   if (has_true) {
#'     result <- make_enrich_dotplot(all_db_results, names(all_term2genes), 
#'                                   plot_title)
#'   } else {
#'     stop("No database led to an enrichment result!")
#'   }
#' 
#'   return(list(dotplot = result[[1]], dotplot_nrows = result[[2]], 
#'               full_enrich_results = result[[3]], raw_results = raw_results))
#' 
#' }