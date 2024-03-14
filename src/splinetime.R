# This R script contains the functions for the package splinetime


#----- load libraries -------------------------------------


# ==============================================================================
# Section 1: Internal functions
# ==============================================================================


# ==============================================================================
# Section 2: Exported package functions
# ==============================================================================


# There needs to be an .RData file containing the dfs data, meta, and genes.
st_load <- function () {
  
}


# This allows to "screen" which hyperparameters are best and reports it
st_screen <- function () {
  
} 


# This performs the final limma spline with the determined hyperparameters
st_fit <- function () {
  
} 


# This performs the clustering based on the coefficients and reports it
st_cluster <- function () {
  
} 


# This performs the F statistic ranking as preparation for fgsea
st_rank <- function () {
  
} 

# ------------------------------ Enrichment ------------------------------------
# These package functions perform the enrichment and generate a report


st_enrichr <- function () {
  
} 


st_clusterprofiler <- function () {
  
} 


st_fgsea <- function () {
  
} 