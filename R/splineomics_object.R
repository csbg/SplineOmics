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
#' @param report_info A named list describing the experiment.
#'   Must include the following fields:
#'     - \code{"omics_data_type"}
#'     - \code{"data_description"}
#'     - \code{"data_collection_date"}
#'     - \code{"analyst_name"}
#'     - \code{"contact_info"}
#'     - \code{"project_name"}
#'
#'   May also include the following optional fields:
#'     - \code{"method_description"}
#'     - \code{"results_summary"}
#'     - \code{"conclusions"}
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
#' @param use_array_weights Boolean flag indicating if the robust fit strategy
#' to deal
#' with heteroscedasticity should be used or not. If set to NULL, then this is
#' handeled implicitly based on the result of the Levene test. If this test is
#' significant for at least 10% of the features, then the robust strategy is
#' used. The robust strategy uses the function voomWithQualityWeights for
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
#' 'integrated' for the mode. Isolated means limma determines the results for
#' each level
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
#' @param bp_cfg A named numeric vector specifying the parallelization
#'   configuration, with expected names `"n_cores"` and `"blas_threads"`.
#'
#'   This controls how many **R worker processes** (`n_cores`) and how many
#'   **BLAS/OpenBLAS threads per process** (`blas_threads`) should be used
#'   during parallel computation.
#'
#'   If `bp_cfg` is `NULL`, missing, or any of its required fields is
#'   `NA`, both `n_cores` and `blas_threads` default to `1`. This effectively
#'   disables parallelization and avoids oversubscription of CPU threads.
#'
#' @return A SplineOmics object.
#'
#' @examples
#' set.seed(1)
#'
#' # 6 samples, 4 features
#' toy_data <- matrix(
#'     rnorm(4 * 6, mean = 0, sd = 1),
#'     nrow = 4, ncol = 6,
#'     dimnames = list(
#'         paste0("gene", 1:4),
#'         paste0("S", 1:6)
#'     )
#' )
#'
#' # Sample metadata
#' toy_meta <- data.frame(
#'     SampleID = colnames(toy_data),
#'     Time = c(0, 0, 1, 1, 2, 2),
#'     Condition = factor(c("Ctrl", "Ctrl", "Ctrl", "Trt", "Trt", "Trt"),
#'         levels = c("Ctrl", "Trt")
#'     ),
#'     Batch = factor(c("B1", "B1", "B1", "B2", "B2", "B2")),
#'     stringsAsFactors = FALSE,
#'     row.names = colnames(toy_data)
#' )
#'
#' # Condition vector (must align with samples)
#' cond <- toy_meta$Condition
#'
#' # Minimal annotation (feature-level info)
#' toy_anno <- data.frame(
#'     feature_id = rownames(toy_data),
#'     symbol = c("G1", "G2", "G3", "G4"),
#'     stringsAsFactors = FALSE,
#'     row.names = rownames(toy_data)
#' )
#'
#' # Spline parameters (natural splines with df = 3)
#' toy_spline <- list(spline_type = "n", dof = 3)
#'
#' # Parallel config (single-threaded for examples)
#' toy_bp <- c(n_cores = 1, blas_threads = 1)
#'
#' # Dream params example (optional)
#' toy_dream <- list(dof = 3L, KenwardRoger = FALSE)
#'
#' # Simple design matrix (intercept + condition + time)
#' toy_design <- stats::model.matrix(~ Condition + Time, data = toy_meta)
#'
#' # Required report fields
#' toy_report <- list(
#'     omics_data_type = "RNA-seq (toy)",
#'     data_description = "Simulated expression matrix (4x6)",
#'     data_collection_date = "2025-10-07",
#'     analyst_name = "Analyst A",
#'     contact_info = "analyst@example.org",
#'     project_name = "SplineOmics Demo",
#'     method_description = "Toy example to construct a SplineOmics object"
#' )
#'
#' so <- create_splineomics(
#'     data                 = toy_data,
#'     meta                 = toy_meta,
#'     condition            = cond,
#'     rna_seq_data         = NULL, # not used in this toy
#'     annotation           = toy_anno,
#'     report_info          = toy_report,
#'     meta_batch_column    = "Batch",
#'     meta_batch2_column   = NULL,
#'     feature_name_columns = c("feature_id", "symbol"),
#'     design               = toy_design,
#'     use_array_weights    = FALSE,
#'     dream_params         = toy_dream,
#'     mode                 = "isolated",
#'     spline_params        = toy_spline,
#'     padjust_method       = "BH",
#'     bp_cfg               = toy_bp
#' )
#'
#' class(so)
#' str(so, max.level = 1)
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
    use_array_weights = FALSE,
    dream_params = NULL,
    mode = "isolated",
    spline_params = NULL,
    padjust_method = "BH",
    bp_cfg = NULL) {
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
        use_array_weights = use_array_weights,
        dream_params = dream_params,
        mode = mode,
        spline_params = spline_params,
        padjust_method = padjust_method,
        bp_cfg = bp_cfg
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
#' @examples
#' set.seed(1)
#' toy_data <- matrix(rnorm(12),
#'     nrow = 3,
#'     dimnames = list(paste0("gene", 1:3), paste0("S", 1:4))
#' )
#' toy_meta <- data.frame(
#'     SampleID = colnames(toy_data),
#'     Condition = c("Ctrl", "Ctrl", "Trt", "Trt"),
#'     stringsAsFactors = FALSE,
#'     row.names = colnames(toy_data)
#' )
#'
#' so <- create_splineomics(
#'     data = toy_data,
#'     meta = toy_meta,
#'     condition = toy_meta$Condition
#' )
#'
#' # Update the mode and add a new design matrix
#' new_design <- model.matrix(~Condition, data = toy_meta)
#' so_updated <- update_splineomics(so,
#'     mode = "integrated",
#'     design = new_design
#' )
#'
#' str(so_updated, max.level = 1)
#'
#' @export
#'
update_splineomics <- function(
    splineomics,
    ...) {
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
        "fit",
        "homosc_violation_result",
        "use_array_weights",
        "dream_params",
        "mode",
        "spline_params",
        "limma_splines_result",
        "bp_cfg"
    )

    args <- list(...)

    for (name in names(args)) {
        if (!(name %in% allowed_fields)) {
            stop(paste0("Field '", name, "' is not allowed."))
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
#' @examples
#' # Example: create and print a SplineOmics object
#' set.seed(1)
#' toy_data <- matrix(rnorm(12),
#'     nrow = 3,
#'     dimnames = list(paste0("gene", 1:3), paste0("S", 1:4))
#' )
#' toy_meta <- data.frame(
#'     SampleID = colnames(toy_data),
#'     Condition = c("Ctrl", "Ctrl", "Trt", "Trt"),
#'     stringsAsFactors = FALSE,
#'     row.names = colnames(toy_data)
#' )
#'
#' so <- create_splineomics(
#'     data = toy_data,
#'     meta = toy_meta,
#'     condition = toy_meta$Condition,
#'     spline_params = list(spline_type = "n", dof = 3),
#'     padjust_method = "BH"
#' )
#'
#' # The print method is automatically called:
#' so
#'
#' # Or explicitly:
#' print(so)
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
