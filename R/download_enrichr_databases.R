#' download_enrichr_databases()
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
#' @param output_dir A character string specifying the output directory
#'                   where the .tsv file will be saved. Defaults to the current
#'                   working directory.
#' @param filename Name of the output file (with file extension. Due to commas
#'                 present in some terms, .tsv is recommended). When left out,
#'                 the file is named all_databases_{timestamp}.tsv.
#'
#' @return This function does not return a value but saves a .tsv file in the
#'         specified directory containing the gene sets from the specified
#'         Enrichr databases.
#'
#' @importFrom here here
#' @importFrom rlang .data
#'
#' @export
#'
download_enrichr_databases <- function(
    gene_set_lib,
    output_dir = here::here(),
    filename = NULL
    ) {
  
  # Control the user inputs
  if (!is.character(gene_set_lib) || length(gene_set_lib) == 0) {
    stop_call_false("gene_set_lib must be a character vector with length > 0")
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

  genesets <- genesets |>
    mutate(Gene = gsub(",.+$", "", .data$Gene))

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
      "all_databases_",
      timestamp, ".tsv"
    )
  }

  filename_path <- here::here(
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
    "Download complete! The file has been saved as: ",
    filename_path
    )

  return(invisible(filename_path))
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
  
  pb <- create_progress_bar(
    databases,
    message = "Downloading"
  )

  setNames(lapply(databases, function(dbx) {
    # Update the progress bar
    pb$tick()

    fpath <- paste0(
      "http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?",
      "mode=text&libraryName=",
      dbx
    )

    fhandle <- file(fpath)
    dblines <- tryCatch(
      {
        readLines(con = fhandle)
      },
      error = function(e) {
        message(e, "\nFailed reading database: ", dbx)
        NULL
      }
    )
    close(fhandle)

    if (is.null(dblines)) {
      return(list())
    } else {
      res <- strsplit(dblines, "\t")
      names(res) <- vapply(res, function(x) x[1], character(1))
      res <- lapply(res, function(x) x[3:length(x)])
      return(res)
    }
  }), databases)
}
