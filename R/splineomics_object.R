# Exported functions -----------------------------------------------------------


#' create_splineomics()
#'
#' @description
#' Creates a SplineOmics object containing variables that are commonly used
#' across multiple functions in the package. This object is then passed as an
#' argument to the other functions of this package. 
#'
#' @param data The actual omics data. In the case the rna_seq_data argument is
#' used, still provide this argument. In that case, input the data matrix in
#' here (for example the $E part of the voom object). Assign your feature names
#' as row headers (otherwise, just numbers will be your feature names).
#' @param meta Metadata associated with the omics data.
#' @param condition A condition variable.
#' @param rna_seq_data An object containing the preprocessed RNA-seq data,
#' such as the output from `limma::voom` or a similar preprocessing pipeline.
#' This argument is not subjected to input control.
#' Rather, in that regard it relies on the input control from the `limma::lmfit`
#' function.
#' @param annotation A dataframe with the feature descriptions of data
#' (optional).
#' @param report_info A list containing report information such as omics data
#' type, data description, data collection date, analyst name, contact info,
#' and project name (optional).
#' @param meta_batch_column Column for meta batch information (optional).
#' @param meta_batch2_column Column for secondary meta batch information
#' (optional).
#' @param feature_name_columns Character vector containing the column names of
#'                             the annotation info that describe the features.
#'                             This argument is used to specify in the HTML
#'                             report how exactly the feature names displayed
#'                             above each individual spline plot have been
#'                             created. Use the same vector that was used to
#'                             create the row headers for the data matrix!
#' @param design A design matrix or similar object (optional).
#' @param robust_fit Boolean flag indicating if the robust fit strategy to deal
#' with heteroscedasticity should be used or not. If ommited (parameter value
#' is NULL) then this is handeled implicitly based on the result of the Wilcoxon
#' signed exact test (tests if there is a significant differene in the variance
#' across samples between at least two levels of the experiment). If this test
#' is significant for at least one pairwise comparison, then the robust strategy
#' is applied. The robust strategy uses the function voomWithQualityWeights for
#' RNA-seq data instead of the normal voom function. For other, non-count-based
#' data, the function limma::arrayWeights is used instead, combined with setting
#' the robust argument to TRUE for the limma::eBayes function. In summary, the
#' strategy employed by those functions is to downweights samples with higher
#' variance. This can be neccessary, because linear models have the assumption
#' of homoscedasticity, which means that the variance is (approx.) the same 
#' across all datapoints where the linear model is fitted. If this is violated,
#' then the resulting p-values cannot be trusted (common statistical wisdom).
#' @param dream_params #' A named list or NULL. When not NULL, it can contain
#'  the following named elements:
#' - `dof`: An integer greater than 1, specifying the degrees of freedom for 
#'   the dream topTable.
#' - `KenwardRoger`: A boolean indicating whether to use the Kenward-Roger 
#'   approximation for mixed models.
#' Note that random effects are now directly specified in the design formula 
#' and not in `dream_params`.
#' @param mode For the design formula, you must specify either 'isolated' or
#' 'integrated'. Isolated means limma determines the results for each level
#' using only the data from that level. Integrated means limma determines the
#'  results for all levels using the full dataset (from all levels).
#' @param spline_params Parameters for spline functions (optional). Must contain
#' the named elements spline_type, which must contain either the string "n" for
#' natural cubic splines, or "b", for B-splines, the named element degree in the
#' case of B-splines, that must contain only an integer, and the named element
#' dof, specifying the degree of freedom, containing an integer and required
#' both for natural and B-splines.
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
    rna_seq_data = NULL,
    annotation = NULL,
    report_info = NULL,
    meta_batch_column = NULL,
    meta_batch2_column = NULL,
    feature_name_columns = NULL,
    design = NULL,
    robust_fit = NULL,
    dream_params = NULL,
    mode = NULL,
    spline_params = NULL,
    padjust_method = "BH"
    ) {
  
  splineomics <- list(
    data = data,
    rna_seq_data = rna_seq_data,
    meta = meta,
    condition = condition,
    annotation = annotation,
    report_info = report_info,
    meta_batch_column = meta_batch_column,
    meta_batch2_column = meta_batch2_column,
    feature_name_columns = feature_name_columns,
    design = design,
    robust_fit = robust_fit,
    dream_params = dream_params,
    mode = mode,
    spline_params = spline_params,
    padjust_method = padjust_method
  )

  class(splineomics) <- "SplineOmics"
  return(splineomics)
}


#' update_splineomics()
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
    stop("The passed object must be of class 'SplineOmics'")
  }

  allowed_fields <- c(
    "data",
    "rna_seq_data",
    "meta",
    "condition",
    "annotation",
    "report_info",
    "meta_batch_column",
    "meta_batch2_column",
    "feature_name_columns",
    "design",
    "robust_fit",
    "dream_params",
    "mode",
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


#' Print function for SplineOmics objects
#'
#' @description
#' This function provides a summary print of the SplineOmics object, showing
#' relevant information such as the number of features, samples, metadata,
#' RNA-seq data, annotation, and spline parameters.
#'
#' @param x A SplineOmics object created by the `create_splineomics` function.
#' @param ... Additional arguments passed to or from other methods.
#'
#' @details
#' This function is automatically called when a SplineOmics object is printed.
#' It provides a concise overview of the object's contents and attributes,
#' including the dimensions of the data, available metadata, and other relevant
#' information such as annotations and spline parameters.
#'
#' @return
#' The function does not return a value. It prints a summary of
#' the SplineOmics object.
#'
#' @importFrom utils head
#'
#' @export
#'
print.SplineOmics <- function(x, ...) {
  cat("data:")
  cat("SplineOmics Object\n")
  cat("-------------------\n")

  # Print summary information
  cat("Number of features (rows):", nrow(x$data), "\n")
  cat("Number of samples (columns):", ncol(x$data), "\n")

  cat("Meta data columns:", ncol(x$meta), "\n")
  cat("First few meta columns:\n")
  print(utils::head(x$meta, 3))

  cat("Condition:", x$condition, "\n")

  if (!is.null(x$rna_seq_data)) {
    cat("RNA-seq data is provided.\n")
  } else {
    cat("No RNA-seq data provided.\n")
  }

  if (!is.null(x$annotation)) {
    cat("Annotation provided with", nrow(x$annotation), "entries.\n")
  } else {
    cat("No annotation provided.\n")
  }

  if (!is.null(x$spline_params)) {
    cat("Spline parameters are set:\n")
    print(x$spline_params)
  } else {
    cat("No spline parameters set.\n")
  }

  cat("P-value adjustment method:", x$padjust_method, "\n")
}
