#' Download Gene Set Annotations from Bioconductor Organism Databases
#' 
#' @noRd
#' 
#' @description
#' This function extracts gene-to-ontology mappings from a specified 
#' Bioconductor organism annotation package (e.g., `org.Hs.eg.db`, 
#' `org.Mm.eg.db`) and saves the gene sets to a `.tsv` file in a 
#' standardized format. The output includes mappings for Gene Ontology 
#' (GO) Biological Process (BP), Molecular Function (MF), Cellular 
#' Component (CC), and KEGG pathways. The resulting file can be used 
#' directly with enrichment functions such as `clusterProfiler::enricher()` 
#' with `TERM2GENE`.
#'
#' @param organism_db A string specifying the Bioconductor organism 
#' annotation database to use (e.g., `"org.Hs.eg.db"` for human or 
#' `"org.Mm.eg.db"` for mouse).
#' @param output_dir A string specifying the output directory where the 
#' `.tsv` file will be saved. Defaults to the current project directory 
#' as defined by `here::here()`.
#' @param filename An optional string specifying the filename for the 
#' output file. If `NULL` (default), a filename is generated automatically 
#' with a timestamp.
#'
#' @return Invisibly returns the full file path to the saved `.tsv` file.
#'
#' @details The output file contains three columns: `DB` (ontology source, 
#' e.g., "GO_BP", "KEGG"), `Geneset` (ontology term ID or pathway ID), and 
#' `Gene` (gene symbol). This format mirrors outputs from Enrichr and is 
#' suitable for use in custom enrichment workflows.
#'
#' @examples
#' \dontrun{
#' # Download mouse gene sets and save as .tsv
#' download_bioconductor_database(
#'     organism_db = "org.Mm.eg.db",
#'     output_dir = "results/annotation"
#' )
#' }
#'
download_bioconductor_database <- function(
    organism_db = "org.Hs.eg.db",
    output_dir = here::here(),
    filename = NULL
) {
  
  # Load required packages
  if (!requireNamespace(organism_db, quietly = TRUE)) {
    stop_call_false("Organism database ", organism_db, " is not installed.")
  }
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop_call_false("Package 'AnnotationDbi' is required but not installed.")
  }
  
  orgdb <- get(organism_db, envir = asNamespace(organism_db))
  
  # Fetch GO BP, MF, CC
  go_df <- AnnotationDbi::select(
    orgdb,
    keys = AnnotationDbi::keys(orgdb, keytype = "ENTREZID"),
    columns = c("SYMBOL", "GO", "ONTOLOGY"),
    keytype = "ENTREZID"
  )
  
  go_df <- go_df[!is.na(go_df$GO), ]
  go_df <- go_df |>
    dplyr::mutate(
      DB = dplyr::case_when(
        ONTOLOGY == "BP" ~ "GO_BP",
        ONTOLOGY == "MF" ~ "GO_MF",
        ONTOLOGY == "CC" ~ "GO_CC",
        TRUE ~ NA_character_
      ),
      Geneset = GO,
      Gene = SYMBOL
    ) |>
    dplyr::select(DB, Geneset, Gene) |>
    dplyr::filter(!is.na(DB) & !is.na(Gene))
  
  # Fetch KEGG
  kegg_df <- AnnotationDbi::select(
    orgdb,
    keys = AnnotationDbi::keys(orgdb, keytype = "ENTREZID"),
    columns = c("SYMBOL", "PATH"),
    keytype = "ENTREZID"
  )
  
  kegg_df <- kegg_df[!is.na(kegg_df$PATH), ]
  kegg_df <- kegg_df |>
    dplyr::mutate(
      DB = "KEGG",
      Geneset = PATH,
      Gene = SYMBOL
    ) |>
    dplyr::select(DB, Geneset, Gene) |>
    dplyr::filter(!is.na(Gene))
  
  # Combine all
  genesets <- dplyr::bind_rows(go_df, kegg_df)
  
  # Create output dir
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Generate filename if NULL
  if (is.null(filename)) {
    timestamp <- format(Sys.time(), "%d_%m_%Y-%H_%M_%S")
    filename <- paste0(
      "bioconductor_",
      organism_db,
      "_",
      timestamp, 
      ".tsv"
      )
  }
  
  filename_path <- here::here(output_dir, filename)
  
  # Write output
  utils::write.table(
    genesets,
    file = filename_path,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
  
  message(
    "\nDownload complete! The file has been saved as: ",
    filename_path
  )
  
  return(invisible(filename_path))
}
