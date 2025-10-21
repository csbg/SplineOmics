#' Download the Enrichr databases
#'
#' @description
#' This function downloads gene sets from specified Enrichr databases and saves
#' them to a specified output directory as a .tsv file per default. The file is
#' named with a timestamp per default to ensure uniqueness (all databases are
#' stored in a single file). This file has 3 columns: DB containing the database
#' name, Geneset, with the genesets, and Gene, with the gene names.
#'
#' @param gene_set_lib A character vector of database names to download from
#'                     Enrichr, for example: c("WikiPathways_2019_Human",
#'                     "NCI-Nature_2016",)
#' @param output_dir A string specifying the output directory where the
#' `.tsv` file will be saved. Defaults to the current project directory
#' as defined by `here::here()`.
#'
#' @param filename Name of the output file (with file extension. Due to commas
#'                 present in some terms, .tsv is recommended). When left out,
#'                 the file is named all_databases_timestamp.tsv.
#'
#' @return
#' A `data.frame` of gene set annotations with three columns:
#' \describe{
#'   \item{DB}{Database name (e.g. `"WikiPathways_2019_Human"`,
#'   `"NCI-Nature_2016"`).}
#'   \item{Geneset}{The gene set or pathway term from that database.}
#'   \item{Gene}{A gene contained in the gene set.}
#' }
#'
#' In addition to returning the `data.frame`, the function also writes the same
#' table to disk as a `.tsv` file in the specified `output_dir`.
#'
#' @importFrom here here
#' @importFrom rlang .data
#' @importFrom tibble tibble
#' @importFrom dplyr mutate
#' @importFrom utils write.table
#'
#' @export
#'
download_enrichr_databases <- function(
    gene_set_lib,
    output_dir = here::here(),
    filename = NULL) {
    # Control the user inputs
    if (!is.character(gene_set_lib) || length(gene_set_lib) == 0) {
        stop_call_false(
            "gene_set_lib must be a character vector with length > 0"
            )
    }

    # Control the filename input (must be NULL or a valid string)
    if (!is.null(filename) &&
        (!is.character(filename) ||
            length(filename) != 1)) {
        stop_call_false("filename must be a single string or NULL.")
    }

    # Control the inputs
    args <- lapply(as.list(match.call()[-1]), eval, parent.frame())
    check_null_elements(args)
    input_control <- InputControl$new(args)
    input_control$auto_validate()

    genesets <- enrichr_get_genesets(databases = gene_set_lib)

    genesets <- do.call(rbind, lapply(names(genesets), function(db.nam) {
        do.call(rbind, lapply(names(genesets[[db.nam]]), function(set.nam) {
            tibble::tibble(
                DB = db.nam,
                Geneset = set.nam,
                Gene = genesets[[db.nam]][[set.nam]]
            )
        }))
    }))

    # When the user is offline, genesets won't be available, raising a weird
    # error
    if (!exists("genesets") || is.null(genesets)) {
        stop_call_false(
            "Object `genesets` is missing or NULL - are you online?"
            )
    }
    genesets <- genesets |>
        dplyr::mutate(Gene = gsub(",.+$", "", .data$Gene))

    dir.create(
        output_dir,
        recursive = TRUE,
        showWarnings = FALSE
    )

    # Create filename if not specified
    if (is.null(filename)) {
        timestamp <- format(
            Sys.time(),
            "%d_%m_%Y-%H_%M_%S"
        )
        filename <- paste0(
            "enrichr_databases_",
            timestamp, ".tsv"
        )
    }

    filename_path <- file.path(
        output_dir,
        filename
    )

    utils::write.table(
        x = genesets,
        file = filename_path,
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE # Do not quote strings
    )

    message(
        "\nDownload complete! The file has been saved as: ",
        filename_path
    )

    genesets
}


# Level 1 internal functions ---------------------------------------------------


#' Get Enrichr Gene Sets
#'
#' @noRd
#'
#' @description
#' This function downloads gene sets from specified Enrichr databases.
#' It returns a list where each element is a list corresponding to a database,
#' with each element containing a vector of human gene symbols for a gene set.
#'
#' @param databases A character vector of database names to download from
#'                  Enrichr.
#'
#' @return A named list of gene sets from the specified Enrichr databases. Each
#'         database is represented as a list, with gene set names as list
#'         names and vectors of human gene symbols as list elements.
#'
enrichr_get_genesets <- function(databases) {
    pb <- create_progress_bar(databases, message = "Downloading")

    failed <- character(0)
    success <- character(0) # <-- track successes
    results <- setNames(vector("list", length(databases)), databases)

    for (dbx in databases) {
        pb$tick()

        url <- paste0(
            "https://maayanlab.cloud/Enrichr/geneSetLibrary?",
            "mode=text&libraryName=", dbx
        )

        dblines <- tryCatch(
            readLines(url, warn = FALSE),
            error = function(e) NULL
        )

        if (!is.null(dblines) && length(dblines) > 0) {
            res <- strsplit(dblines, "\t")
            names(res) <- vapply(res, `[`, 1, FUN.VALUE = character(1))
            res <- lapply(res, function(x) x[-c(1, 2)]) # keep only genes
            results[[dbx]] <- res
            success <- c(success, dbx) # <-- record success
        } else {
            failed <- c(failed, dbx)
        }
    }

    # Keep only successfully downloaded libraries
    results <- results[success]

    # User-friendly summary
    if (length(success) || length(failed)) {
        msg <- c("\nEnrichr download summary:")
        if (length(success)) {
            msg <- c(
                msg,
                "  [OK] downloaded : ",
                paste(success, collapse = ", ")
                )
        }
        if (length(failed)) {
            msg <- c(
                msg, "  [!!] not found  : ", paste(failed, collapse = ", "),
                "\n      (Check spelling or run enrichR::listEnrichrDbs() ",
                "for valid names)"
            )
        }
        message(paste(msg, collapse = "\n"))
    }

    invisible(results)
}
