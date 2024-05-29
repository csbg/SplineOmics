# ## ---------------------------
# ##
# ## Script name: download_enrichr_databases
# ##
# ## Purpose of script: selects & downloads Enrichr databases for fGSEA analysis
# ##
# ## Author: Dr. Veronika Schäpertöns & Prof. Nikolaus Fortelny
# ##
# ## Date Created: 24.03.2023
# ##
# ## Copyright (c) Veronika Schäpertöns, 2023
# ## Email: veronika.schaepertoens@plus.ac.at
# ##
# ## ---------------------------
# ##
# ## Notes:
# ##   
# ##
# ## 
# 
# 
# # libraries ----------------------------------------------------------------
# 
# library(tidyverse)
# library(data.table)
# library(fs)
# 
# 
# # specify which databases -------------------------------------------------
# 
# gene_set_lib <- c("WikiPathways_2019_Human",
#                   "NCI-Nature_2016",
#                   "TRRUST_Transcription_Factors_2019",
#                   "MSigDB_Hallmark_2020",
#                   "GO_Cellular_Component_2018",
#                   "CORUM",
#                   "KEGG_2019_Human",
#                   "TRANSFAC_and_JASPAR_PWMs",
#                   "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
#                   "GO_Biological_Process_2018",
#                   "GO_Molecular_Function_2018",
#                   "Human_Gene_Atlas"
# )
# 
# # get enrichr gene sets ---------------------------------------------------
# # This function downloads databases from Enrichr. Returns:
# # - A list, where each element is
# # -- another list (names correspond to the individual databases), where each element is
# # ---- a vector (names correspond to gene sets) where each element is a 
# # ------ human gene symbol (character)
# enrichrGetGenesets <- function(databases) {
#   
#   setNames(lapply(databases, function(dbx){
#     fpath <- paste0("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=",dbx)
#     fhandle <- file(fpath)
#     dblines <- tryCatch({
#       readLines(con = fhandle)
#     }, error = function(e){
#       message(e, "\nFailed reading database: ", dbx)
#       NULL
#     })
#     close(fhandle)
#     if (is.null(dblines)) {
#       return(list())
#     }else {
#       res <- strsplit(dblines, "\t")
#       names(res) <- sapply(res, function(x) x[1])
#       res <- lapply(res, function(x) x[3:length(x)])
#       return(res)
#     }
#   }), databases)
# }
# 
# # download genesets -------------------------------------------------------
# 
# genesets <- enrichrGetGenesets(databases = gene_set_lib)
# 
# genesets <- do.call(rbind, lapply(names(genesets), function(db.nam){
#   do.call(rbind,lapply(names(genesets[[db.nam]]), function(set.nam){
#     data.table(DB = db.nam, Geneset = set.nam, Gene = genesets[[db.nam]][[set.nam]])
#   }))
# }))
# 
# genesets[,Gene := gsub(",.+$", "", Gene)]
# 
# fs::dir_create("databases")
# 
# timestamp <- format(Sys.time(), "%d_%m_%Y-%H_%M_%S")
# filename <- paste0("all_databases_", timestamp, ".tsv")
# write_tsv(x = genesets, paste0("databases/", filename))