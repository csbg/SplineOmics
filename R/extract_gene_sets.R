#' Extract gene set annotations from Bioconductor organism databases
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
#' @param organism_db `character(1)`: A string specifying the Bioconductor 
#' organism annotation database to use (e.g., `"org.Hs.eg.db"` for human or 
#' `"org.Mm.eg.db"` for mouse).
#'
#' @param output_dir `character(1)`: A string specifying the output directory 
#' where the `.tsv` file will be saved.
#'
#' @param filename `character(1)` | `NULL`: An optional string specifying the 
#' filename for the output file. If `NULL` (default), a filename is generated 
#' automatically with a timestamp.
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
#' such cases the KEGG mappings may be absent.
#'
#' In addition to returning the `data.frame`, the function also writes the same
#' table to disk as a `.tsv` file in the specified `output_dir`.
#'
#' @importFrom rlang .data abort
#' @importFrom dplyr mutate case_when select filter bind_rows distinct
#'
#' @examples
#' # Minimal real example (runs only if org package is installed)
#' tmp <- tempdir()
#' if (requireNamespace("org.Mm.eg.db", quietly = TRUE) &&
#'     requireNamespace("AnnotationDbi", quietly = TRUE)) {
#'     gs <- extract_gene_sets(
#'         organism_db = "org.Mm.eg.db",
#'         output_dir  = tmp,
#'         filename    = "mm_genesets.tsv"
#'     )
#'     head(gs)
#'     # The file path:
#'     file.path(tmp, "mm_genesets.tsv")
#' }
#'
#' # If the organism package is not installed, you can still see the TSV format:
#' tiny <- data.frame(
#'     DB = c("GO_BP", "GO_MF"),
#'     Geneset = c("GO:0008150", "GO:0003674"),
#'     Gene = c("Trp53", "Egfr"),
#'     stringsAsFactors = FALSE
#' )
#' utils::write.table(
#'     tiny,
#'     file = file.path(tmp, "example_genesets.tsv"),
#'     sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
#' )
#'
#' @export
#'
extract_gene_sets <- function(
    organism_db = "org.Hs.eg.db",
    output_dir = tempdir(),
    filename = NULL
    ) {
    # Dependencies (Not core of SplineOmics, must be downloaded for this
    # function)
    if (!requireNamespace(organism_db, quietly = TRUE)) {
        rlang::abort(
            paste0(
                "Organism database '",
                organism_db,
                "' is not installed."
            )
        )
    }
    
    if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
        rlang::abort(
            "Package 'AnnotationDbi' is required but not installed."
        )
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
    if (!length(entrez_keys)) {
        rlang::abort(
            "No ENTREZID keys were retrieved from the organism database."
            )
    }

    # GO (BP/MF/CC)
    go_df <- AnnotationDbi::select(
        orgdb,
        keys = entrez_keys,
        columns = c(
            "SYMBOL",
            "GO",
            "ONTOLOGY"
        ),
        keytype = "ENTREZID",
        check = FALSE
    )
    if (!nrow(go_df)) {
        rlang::abort(
            "No GO mappings were retrieved from the organism database."
            )
    }

    go_df <- go_df[!is.na(go_df$GO), , drop = FALSE]
    go_df <- dplyr::mutate(
        go_df,
        DB = dplyr::case_when(
            .data$ONTOLOGY == "BP" ~ "GO_BP",
            .data$ONTOLOGY == "MF" ~ "GO_MF",
            .data$ONTOLOGY == "CC" ~ "GO_CC",
            TRUE ~ NA_character_
        ),
        Geneset = .data$GO,
        Gene = .data$SYMBOL
    )
    keep <- intersect(c(
        "DB",
        "Geneset",
        "Gene"
        ), names(go_df))
    go_df <- dplyr::select(
        go_df,
        dplyr::all_of(keep)
    )
    go_df <- dplyr::filter(
        go_df,
        !is.na(.data$DB), 
        !is.na(.data$Gene)
        )

    kegg_df <- NULL
    if ("PATH" %in% AnnotationDbi::columns(orgdb)) {
        kegg_raw <- AnnotationDbi::select(
            orgdb,
            keys    = entrez_keys,
            columns = c("SYMBOL", "PATH"),
            keytype = "ENTREZID",
            check = FALSE
        )
        if (!is.null(kegg_raw) && nrow(kegg_raw)) {
            kegg_raw <- kegg_raw[!is.na(kegg_raw$PATH), , drop = FALSE]
            if (nrow(kegg_raw)) {
                kegg_df <- kegg_raw |>
                    dplyr::mutate(
                        DB      = "KEGG",
                        Geneset = .data$PATH,
                        Gene    = .data$SYMBOL
                    ) |>
                    (\(d) {
                        keep <- intersect(c(
                            "DB",
                            "Geneset",
                            "Gene"
                            ), names(d))
                        dplyr::select(d, dplyr::all_of(keep))
                    })() |>
                    dplyr::filter(!is.na(.data$Gene))
            }
        }
    }

    # Combine
    genesets <- if (!is.null(kegg_df)) {
        dplyr::bind_rows(go_df, kegg_df)
    } else {
        go_df
    }
    
    # De-duplicate (always)
    genesets <- dplyr::distinct(
        genesets,
        .data$DB,
        .data$Geneset,
        .data$Gene
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
        safe_org <- gsub(
            "[^A-Za-z0-9]+",
            "_",
            organism_db
        )
        filename <- paste0(
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
        "\nGene set extraction complete complete! The file has been saved as: ",
        filename_path
    )
    invisible(genesets)
}
