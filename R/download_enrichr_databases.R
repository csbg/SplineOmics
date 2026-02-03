#' Download Enrichr databases
#'
#' @description
#' Download gene sets from one or more Enrichr libraries and return them as a
#' table with one row per gene-set membership.
#'
#' If `output_dir` is provided, the same table is also written to disk as a
#' tab-separated file (`.tsv`). If `filename` is `NULL`, a timestamped name is
#' generated to ensure uniqueness.
#'
#' @param gene_set_lib
#' `character()`. A character vector of Enrichr library names to download,
#' e.g. `c("WikiPathways_2019_Human", "NCI-Nature_2016")`.
#'
#' @param output_dir
#' `character(1)` or `NULL`. Output directory for writing a `.tsv` file. If
#' `NULL`, no file is written.
#'
#' @param filename
#' `character(1)` or `NULL`. Output file name (not a path). If `NULL` and
#' `output_dir` is not `NULL`, a default name of the form
#' `enrichr_databases_YYYYmmdd-HHMMSS.tsv` is used. Due to commas in some
#' terms, `.tsv` is recommended.
#'
#' @return
#' A `data.frame` with three columns:
#' \describe{
#'   \item{DB}{Enrichr library name.}
#'   \item{Geneset}{Gene set or pathway term within that library.}
#'   \item{Gene}{Gene symbol contained in the gene set.}
#' }
#'
#' If a requested library cannot be downloaded, it may be omitted from the
#' result. The function errors if no gene sets can be retrieved.
#'
#' @examples
#' if (interactive()) {
#'   libs <- c("WikiPathways_2019_Human")
#'   out <- download_enrichr_databases(
#'     gene_set_lib = libs,
#'     output_dir = tempdir(),
#'     filename = "enrichr_demo.tsv"
#'   )
#'   head(out)
#' }
#'
#' @importFrom rlang .data abort
#' @importFrom tibble tibble
#' @importFrom dplyr mutate bind_rows
#' @importFrom utils write.table
#'
#' @export
#'
download_enrichr_databases <- function(
        gene_set_lib,
        output_dir = NULL,
        filename = NULL
) {
    check_dwld_enrichdb_args(
        gene_set_lib = gene_set_lib,
        output_dir = output_dir,
        filename = filename
    )
    genesets_list <- tryCatch(
        enrichr_get_genesets(databases = gene_set_lib),
        error = function(e) {
            rlang::abort(
                c(
                    "Failed to download Enrichr gene set libraries.",
                    "i" = paste(
                        "This can happen if you are offline or Enrichr is",
                        "unreachable."
                    ),
                    "x" = conditionMessage(e)
                )
            )
        }
    )
    
    if (is.null(genesets_list) || length(genesets_list) == 0L) {
        rlang::abort(
            c(
                "No Enrichr gene sets were retrieved.",
                "i" = paste(
                    "You may be offline, the endpoint may be down, or the",
                    "library names may be invalid."
                )
            )
        )
    }
    
    genesets <- do.call(
        dplyr::bind_rows,
        lapply(names(genesets_list), function(db.nam) {
            db <- genesets_list[[db.nam]]
            if (is.null(db) || length(db) == 0L) {
                return(NULL)
            }
            
            dplyr::bind_rows(lapply(names(db), function(set.nam) {
                genes <- db[[set.nam]]
                if (is.null(genes) || length(genes) == 0L) {
                    return(NULL)
                }
                tibble::tibble(
                    DB = db.nam,
                    Geneset = set.nam,
                    Gene = genes
                )
            }))
        })
    )
    
    if (is.null(genesets) || nrow(genesets) == 0L) {
        rlang::abort(
            c(
                "Downloaded libraries contained no gene sets after parsing.",
                "i" = paste(
                    "This can happen if the libraries are empty or the",
                    "download response format changed."
                )
            )
        )
    }
    
    genesets <- genesets |>
        dplyr::mutate(
            Gene = gsub(",.+$", "", .data$Gene)
        )
    
    if (!is.null(output_dir)) {
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        
        if (is.null(filename)) {
            timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
            filename <- paste0("enrichr_databases_", timestamp, ".tsv")
        }
        
        filename_path <- file.path(output_dir, filename)
        
        utils::write.table(
            x = genesets,
            file = filename_path,
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE
        )
        
        message(
            "Download complete! The file has been saved as: ",
            filename_path
        )
    }
    
    genesets
}


# Level 1 internal functions ---------------------------------------------------


#' Validate inputs for Enrichr database download
#'
#' @noRd
#'
#' @description
#' Validate and sanitize inputs used by the Enrichr database download
#' functions.
#'
#' This function performs defensive checks on the provided arguments to
#' ensure they are well-formed and safe to use for downloading data and
#' optionally writing results to disk. Invalid inputs result in an error.
#'
#' @param gene_set_lib
#' A character vector of Enrichr database names.
#'
#' @param output_dir
#' A length-one character string specifying an output directory, or
#' `NULL` to disable writing results to disk.
#'
#' @param filename
#' A length-one character string specifying an output file name, or
#' `NULL` to use a default name.
#'
#' @return
#' Invisibly returns a list containing the validated input arguments.
#' 
check_dwld_enrichdb_args <- function(
        gene_set_lib,
        output_dir = NULL,
        filename = NULL
) {
    if (!is.character(gene_set_lib) || length(gene_set_lib) < 1L) {
        rlang::abort("`gene_set_lib` must be a non-empty character vector.")
    }
    if (anyNA(gene_set_lib) || any(!nzchar(gene_set_lib))) {
        rlang::abort("`gene_set_lib` must not contain NA or empty strings.")
    }
    if (any(grepl("[\r\n\t]", gene_set_lib))) {
        rlang::abort(
            "`gene_set_lib` must not contain tabs or newlines."
        )
    }
    
    if (!is.null(output_dir)) {
        if (!is.character(output_dir) || length(output_dir) != 1L) {
            rlang::abort(
                "`output_dir` must be NULL or a length-1 character string."
            )
        }
        if (is.na(output_dir) || !nzchar(output_dir)) {
            rlang::abort(
                "`output_dir` must not be NA or an empty string."
            )
        }
        if (grepl("[\r\n\t]", output_dir)) {
            rlang::abort(
                "`output_dir` must not contain tabs or newlines."
            )
        }
    }
    
    if (!is.null(filename)) {
        if (!is.character(filename) || length(filename) != 1L) {
            rlang::abort(
                "`filename` must be NULL or a length-1 character string."
            )
        }
        if (is.na(filename) || !nzchar(filename)) {
            rlang::abort("`filename` must not be NA or an empty string.")
        }
        if (grepl("[\r\n\t]", filename)) {
            rlang::abort("`filename` must not contain tabs or newlines.")
        }
        
        if (grepl("[/\\\\]", filename)) {
            rlang::abort(
                c(
                    "`filename` must be a file name, not a path.",
                    "i" = "Do not include '/' or '\\\\' in `filename`."
                )
            )
        }
        if (grepl("^\\.+$", filename)) {
            rlang::abort(
                "`filename` must not be '.' or '..'."
            )
        }
    }
    
    invisible(
        list(
            gene_set_lib = gene_set_lib,
            output_dir = output_dir,
            filename = filename
        )
    )
}


#' Get Enrichr Gene Sets
#'
#' @noRd
#'
#' @description
#' Download gene sets from specified Enrichr databases.
#'
#' The function returns a named list where each element corresponds to an
#' Enrichr database. Each database element is itself a named list of gene
#' sets, with gene set names as list names and character vectors of gene
#' symbols as values.
#'
#' If a database cannot be downloaded (e.g., due to connectivity issues or an
#' invalid database name), it is omitted from the result.
#'
#' @param databases
#' A character vector of Enrichr database names.
#'
#' @return
#' A named list of gene sets from the successfully downloaded Enrichr
#' databases. Each database is represented as a named list of gene sets, with
#' gene symbols stored as character vectors.
#'
enrichr_get_genesets <- function(databases) {
    if (!requireNamespace("curl", quietly = TRUE)) {
        rlang::abort("Package `curl` is required for Enrichr downloads.")
    }
    
    if (!curl::has_internet()) {
        rlang::abort(
            c(
                "No internet connection detected.",
                "i" = "Enrichr gene set libraries cannot be downloaded."
            )
        )
    }
    
    show_pb <- interactive()
    pb <- NULL
    if (show_pb) {
        pb <- create_progress_bar(databases, message = "Downloading")
    }
    
    failed <- character(0)
    success <- character(0)
    results <- setNames(vector("list", length(databases)), databases)
    
    for (dbx in databases) {
        if (show_pb) {
            pb$tick()
        }
        
        lib <- utils::URLencode(dbx, reserved = TRUE)
        url <- paste0(
            "https://maayanlab.cloud/Enrichr/geneSetLibrary?",
            "mode=text&libraryName=",
            lib
        )
        
        txt <- tryCatch(
            {
                h <- curl::new_handle()
                curl::handle_setopt(
                    h,
                    timeout = 30,
                    useragent = "Bioconductor package"
                )
                res <- curl::curl_fetch_memory(url, handle = h)
                rawToChar(res$content)
            },
            error = function(e) NULL
        )
        
        if (!is.null(txt) && nzchar(txt)) {
            dblines <- strsplit(txt, "\n", fixed = TRUE)[[1]]
            dblines <- dblines[nzchar(dblines)]
        } else {
            dblines <- character(0)
        }
        
        if (length(dblines) > 0L) {
            res <- strsplit(dblines, "\t", fixed = TRUE)
            
            ok <- lengths(res) >= 3L
            res <- res[ok]
            
            if (length(res) > 0L) {
                nm <- vapply(res, `[`, 1L, FUN.VALUE = character(1))
                res <- lapply(res, function(x) x[-c(1L, 2L)])
                names(res) <- nm
                results[[dbx]] <- res
                success <- c(success, dbx)
            } else {
                failed <- c(failed, dbx)
            }
        } else {
            failed <- c(failed, dbx)
        }
    }
    
    results <- results[success]
    
    if ((length(success) || length(failed)) && interactive()) {
        msg <- c("Enrichr download summary:")
        if (length(success)) {
            msg <- c(
                msg,
                paste0("  [OK] downloaded: ", paste(success, collapse = ", "))
            )
        }
        if (length(failed)) {
            msg <- c(
                msg,
                paste0("  [!!] not found: ", paste(failed, collapse = ", "))
            )
            msg <- c(
                msg,
                "      (Check spelling or run enrichR::listEnrichrDbs()"
            )
        }
        message(paste(msg, collapse = "\n"))
    }
    
    invisible(results)
}