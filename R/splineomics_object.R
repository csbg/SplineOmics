#' Create and update the SplineOmics object
#' =====
#'
#' Description
#' -----------
#' Contains the functions to create and update a SplineOmics object. This object
#' is used to collect function arguments, that are equivalent for more than one 
#' exported function of the SplineOmics package. Additionally
#'
#' Functions
#' ---------
#' - create_splineomics: Create a SplineOmics object
#' - update_splineomics: Add additional arguments to the SplineOmics 
#' object or overwrite existing arguments.
#'
#' Classes
#' -------
#' None
#'
#' Notes
#' -----
#' None


# Exported functions -----------------------------------------------------------


#' Create a SplineOmics object
#'
#' @description
#' Creates a SplineOmics object containing variables that are commonly used 
#' across multiple functions in the package.
#' 
#' @param data The actual omics data.
#' @param meta Metadata associated with the omics data.
#' @param condition A condition variable.
#' @param annotation A dataframe with the feature descriptions of data 
#' (optional).
#' @param report_info A list containing report information such as omics data 
#' type, data description, data collection date, analyst name, contact info, 
#' and project name (optional). 
#' @param meta_batch_column Column for meta batch information (optional).
#' @param meta_batch2_column Column for secondary meta batch information 
#' (optional).
#' @param design A design matrix or similar object (optional).
#' @param spline_params Parameters for spline functions (optional).
#' @param preprocess_rna_seq Boolean specifying whether to preprocess RNA seq
#' @param normalization_fun Function used for normalizing RNA-seq. Must take as
#' input the y of: y <- edgeR::DGEList(counts = raw_counts) and output the y 
#' with the normalized counts.
#' @param padjust_method Method for p-value adjustment, one of "none", "BH", 
#' "BY", "holm", "bonferroni", "hochberg", or "hommel". 
#' Defaults to "BH" (Benjamini-Hochberg).
#'
#' @return A SplineOmics object.
#'
#' @export
#' 
create_splineomics <- function(
    data,
    meta,
    condition,
    annotation = NULL,
    report_info = NULL,
    meta_batch_column = NULL,
    meta_batch2_column = NULL,
    feature_name_columns = NULL,
    design = NULL,
    spline_params = NULL,
    preprocess_rna_seq = FALSE,
    normalization_fun = NULL,
    padjust_method = "BH"
    ) {
  
  splineomics <- list(
    data = data,
    preprocess_rna_seq = preprocess_rna_seq,
    meta = meta,
    condition = condition,
    annotation = annotation,
    report_info = report_info,
    meta_batch_column = meta_batch_column,
    meta_batch2_column = meta_batch2_column,
    feature_name_columns = feature_name_columns,
    design = design,
    spline_params = spline_params,
    padjust_method = padjust_method
  )
  
  class(splineomics) <- "SplineOmics"
  return(splineomics)
}


#' Update a SplineOmics object
#'
#' @description
#' Updates a SplineOmics object by modifying existing fields or adding new ones.
#'
#' @param splineomics A SplineOmics object to be updated.
#' @param ... Named arguments with new values for fields to be updated or added.
#'
#' @return The updated SplineOmics object.
#'
#' @export
#' 
update_splineomics <- function(
    splineomics,
    ...
    ) {
  
  if (!inherits(splineomics, "SplineOmics")) {
    stop("The object must be of class 'SplineOmics'")
  }
  
  allowed_fields <- c(
    "data",
    "meta",
    "condition",
    "annotation",
    "report_info",
    "meta_batch_column",
    "meta_batch2_column",
    "feature_name_columns",
    "design",
    "spline_params",
    "limma_splines_result"
    )
  
  args <- list(...)
  for (name in names(args)) {
    if (!(name %in% allowed_fields)) {
      stop(paste("Field", name, "is not allowed."))
    }
    splineomics[[name]] <- args[[name]]
  }
  
  return(splineomics)
}