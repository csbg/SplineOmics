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
#' @param report_info A list containing report information such as omics data 
#' type, data description, data collection date, analyst name, contact info, 
#' and project name.
#' @param condition A condition variable.
#' @param meta_batch_column Column for meta batch information (optional).
#' @param meta_batch2_column Column for secondary meta batch information 
#' (optional).
#' @param design A design matrix or similar object (optional).
#' @param spline_params Parameters for spline functions (optional).
#' @param data The actual omics data (optional).
#' @param meta Metadata associated with the omics data (optional).
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
    design = NULL,
    spline_params = NULL
    ) {
  
  splineomics <- list(
    data = data,
    meta = meta,
    condition = condition,
    annotation = annotation,
    report_info = report_info,
    meta_batch_column = meta_batch_column,
    meta_batch2_column = meta_batch2_column,
    design = design,
    spline_params = spline_params
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