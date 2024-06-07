#' The exported function download_enrichr_databases allows to selects & 
#' download Enrichr databases for local use with clusterProfiler.



# Exported function: download_enrichr_databases() ------------------------------


#' Download Enrichr Databases
#'
#' @description
#' This function downloads gene sets from specified Enrichr databases and saves
#'  them to a specified output directory as a .tsv file. The file is named with 
#'  a timestamp to ensure uniqueness.
#'
#' @param gene_set_lib A character vector of database names to download from 
#'                     Enrichr.
#' @param output_dir A character string specifying the output directory
#'                   where the .tsv file will be saved. Defaults to the current
#'                   working directory.
#'
#' @return This function does not return a value but saves a .tsv file in the 
#'         specified directory containing the gene sets from the specified 
#'         Enrichr databases.
#'
#' @importFrom data.table data.table
#' @importFrom fs dir_create
#' @importFrom here here
#' @importFrom readr write_tsv
#' 
#' @export
#'
download_enrichr_databases <- function(gene_set_lib,
                                       output_dir = here::here()) {
  
  # Control the user inputs
  if (!is.character(gene_set_lib) || length(gene_set_lib) == 0) {
    stop(paste("gene_set_lib must be a character vector with length > 0"))
  }
  
  # Control the inputs
  args <- lapply(as.list(match.call()[-1]), eval, parent.frame())
  check_null_elements(args)
  input_control <- InputControl$new(args)
  input_control$auto_validate()

  
  debug(enrichr_get_genesets)
  genesets <- enrichr_get_genesets(databases = gene_set_lib)
  
  genesets <- do.call(rbind, lapply(names(genesets), function(db.nam){
    do.call(rbind,lapply(names(genesets[[db.nam]]), function(set.nam){
      data.table::data.table(DB = db.nam, Geneset = set.nam, 
                             Gene = genesets[[db.nam]][[set.nam]])
    }))
  }))
  
  genesets[,Gene := gsub(",.+$", "", Gene)]
  
  fs::dir_create("databases")
  
  timestamp <- format(Sys.time(), "%d_%m_%Y-%H_%M_%S")
  filename <- paste0("all_databases_", timestamp, ".tsv")
  filename_path <- here::here(output_dir, filename)
  readr::write_tsv(x = genesets, filename_path)
  
  return(filename_path)
}



# Level 1 internal functions ---------------------------------------------------


#' Get Enrichr Gene Sets
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
  
  pb <- create_progress_bar(databases)
  
  setNames(lapply(databases, function(dbx) {
    
    # Update the progress bar
    pb$tick()
    
    fpath <- paste0("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?", 
                    "mode=text&libraryName=", dbx)
    
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
      
    } else {
      
      res <- strsplit(dblines, "\t")
      names(res) <- sapply(res, function(x) x[1])
      res <- lapply(res, function(x) x[3:length(x)])
      return(res)
    }
  }), databases)
}
