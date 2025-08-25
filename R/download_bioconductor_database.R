#' Download Gene Set Annotations from Bioconductor Organism Databases
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
#' 
#' @param output_dir A string specifying the output directory where the 
#' `.tsv` file will be saved. Defaults to the current project directory 
#' as defined by `here::here()`.
#' 
#' @param filename An optional string specifying the filename for the 
#' output file. If `NULL` (default), a filename is generated automatically 
#' with a timestamp.
#'
#' @return
#' A `data.frame` of gene set annotations with three columns:
#' \describe{
#'   \item{DB}{Ontology/database source, e.g. `"GO_BP"`, `"GO_MF"`, `"GO_CC"`,
#'   or `"KEGG"` (if available).}
#'   \item{Geneset}{Ontology term ID or pathway ID (e.g. GO ID, KEGG ID).}
#'   \item{Gene}{Gene symbol (\code{SYMBOL}).}
#' }
#'
#' @details
#' The TSV has three columns:
#' \describe{
#'   \item{DB}{Ontology/database source, e.g., \code{"GO_BP"}, \code{"GO_MF"},
#'   \code{"GO_CC"}, or \code{"KEGG"} (if available).}
#'   \item{Geneset}{Ontology term ID or pathway ID (e.g., GO ID, KEGG ID).}
#'   \item{Gene}{Gene symbol (\code{SYMBOL}).}
#' }
#' Note: Some \code{org.*.eg.db} packages no longer include KEGG mappings; in
#' such cases the KEGG section will be empty.
#' 
#' In addition to returning the `data.frame`, the function also writes the same
#' table to disk as a `.tsv` file in the specified `output_dir`.
#' 
#' @importFrom here here
#' 
#' @export
#'
download_bioconductor_database <- function(
    organism_db = "org.Hs.eg.db",
    output_dir  = here::here(),
    filename    = NULL
) {
  # Dependencies (Not core of SplineOmics, must be downloaded for this function)
  if (!requireNamespace(organism_db, quietly = TRUE)) {
    stop_call_false(
      "Organism database ",
      organism_db,
      " is not installed."
      )
  }
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop_call_false("Package 'AnnotationDbi' is required but not installed.")
  }
  
  # Load the org.* namespace and get the OrgDb object
  orgdb <- get(
    organism_db,
    envir = asNamespace(organism_db)
    )
  
  entrez_keys <- AnnotationDbi::keys(
    orgdb,
    keytype = "ENTREZID"
    )
  
  # GO (BP/MF/CC)
  go_df <- AnnotationDbi::select(
    orgdb,
    keys    = entrez_keys,
    columns = c(
      "SYMBOL",
      "GO",
      "ONTOLOGY"
      ),
    keytype = "ENTREZID"
  )
  
  go_df <- go_df[!is.na(go_df$GO), , drop = FALSE]
  go_df <- go_df |>
    dplyr::mutate(
      DB = dplyr::case_when(
        ONTOLOGY == "BP" ~ "GO_BP",
        ONTOLOGY == "MF" ~ "GO_MF",
        ONTOLOGY == "CC" ~ "GO_CC",
        TRUE ~ NA_character_
      ),
      Geneset = GO,
      Gene    = SYMBOL
    ) |>
    dplyr::select(DB, Geneset, Gene) |>
    dplyr::filter(!is.na(DB), !is.na(Gene))
  
  # KEGG (guard for missing PATH)
  kegg_df <- NULL
  if ("PATH" %in% AnnotationDbi::columns(orgdb)) {
    kegg_raw <- AnnotationDbi::select(
      orgdb,
      keys    = entrez_keys,
      columns = c("SYMBOL", "PATH"),
      keytype = "ENTREZID"
    )
    if (!is.null(kegg_raw) && nrow(kegg_raw)) {
      kegg_raw <- kegg_raw[!is.na(kegg_raw$PATH), , drop = FALSE]
      if (nrow(kegg_raw)) {
        kegg_df <- kegg_raw |>
          dplyr::mutate(DB = "KEGG", Geneset = PATH, Gene = SYMBOL) |>
          dplyr::select(DB, Geneset, Gene) |>
          dplyr::filter(!is.na(Gene))
      }
    }
  }
  
  # Combine & de-duplicate
  genesets <- if (!is.null(kegg_df)) {
    dplyr::bind_rows(go_df, kegg_df)
  } else {
    go_df
  } |>
    dplyr::distinct(
      DB,
      Geneset,
      Gene
      )
  
  # Output path
  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
    )
  
  if (is.null(filename)) {
    timestamp <- format(
      Sys.time(),
      "%Y_%m_%d-%H_%M_%S"
      )
    safe_org  <- gsub(
      "[^A-Za-z0-9]+",
      "_",
      organism_db
      )
    filename  <- paste0(
      "bioconductor_",
      safe_org,
      "_",
      timestamp,
      ".tsv"
      )
  }
  
  filename_path <- file.path(
    output_dir,
    filename
    )  
  
  # Write TSV
  utils::write.table(
    genesets,
    file      = filename_path,
    sep       = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote     = FALSE
  )
  
  message(
    "\nDownload complete! The file has been saved as: ",
    filename_path
    )
  genesets
}