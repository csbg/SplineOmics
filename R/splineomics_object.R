#' Create an object containing variables often used by SplineOmics functions
#'
#' @description
#' Creates a SplineOmics object containing variables that are commonly used
#' across multiple functions in the package. This object is then passed as an
#' argument to the other functions of this package.
#'
#' @param data `matrix`: The actual omics data. If the `rna_seq_data` argument
#'  is
#' used, still provide this argument. In that case, input the data matrix here
#' (for example, the `$E` part of a `voom` object). Assign your feature names as
#' row headers; otherwise, numeric indices will be used.
#'
#' @param meta `data.frame`: Metadata associated with the omics data.
#'
#' @param condition `character(1)`: Condition variable describing the
#' experimental groups.
#'
#' @param rna_seq_data `list` | `NULL`: An object containing the preprocessed
#' RNA-seq data, such as the output from `limma::voom` or a similar pipeline.
#' This argument is not validated directly; input checks rely on
#' `limma::lmFit()`.
#'
#' @param annotation `data.frame` | `NULL`: Feature annotations (optional)
#' providing descriptive information about each feature in `data`.
#'
#' @param report_info `list`: Named list describing the experiment. Must include
#' the following fields (all `character(1)`):
#'   - `"omics_data_type"`
#'   - `"data_description"`
#'   - `"data_collection_date"`
#'   - `"analyst_name"`
#'   - `"contact_info"`
#'   - `"project_name"`
#'
#'   Optional fields (all `character(1)`):
#'   - `"method_description"`
#'   - `"results_summary"`
#'   - `"conclusions"`
#'
#' @param meta_batch_column `character(1)` | `NULL`: Column name in `meta`
#' specifying batch information (optional).
#'
#' @param meta_batch2_column `character(1)` | `NULL`: Column name in `meta`
#' specifying secondary batch information (optional).
#'
#' @param feature_name_columns `character()`: Vector of column names from
#' `annotation` that describe the features. Used in the HTML report to define
#' how feature names displayed above each spline plot were created. Use the same
#' vector that was used to create the row headers for the data matrix.
#'
#' @param design `matrix` | `NULL`: Design matrix or similar object (optional).
#'
#' @param use_array_weights `logical(1)`: Boolean flag indicating whether to use
#' the robust fitting strategy to handle heteroskedasticity. If `NULL`, this is
#' determined automatically via the Levene test: if at least 10% of features are
#' significant, the robust strategy is enabled. For RNA-seq data, this uses
#' `limma::voomWithQualityWeights()`, otherwise `limma::arrayWeights()` with
#' `robust = TRUE` in `limma::eBayes()`. These approaches down-weight samples
#' with higher variance, improving validity of statistical inference.
#'
#' @param dream_params `list` | `NULL`: Optional named list controlling
#' mixed-model fitting. When not `NULL`, may include:
#'   - `dof` `integer(1)` Degrees of freedom for the DREAM `topTable`.
#'   - `KenwardRoger` `logical(1)` Whether to use the Kenward-Roger correction.
#'
#'   Random effects are specified directly in the design formula, not here.
#'
#' @param mode `character(1)`: Either `"isolated"` or `"integrated"`. Determines
#' whether conditions are analysed independently (`"isolated"`) or jointly
#' (`"integrated"`). The integrated mode fits a single model across all levels.
#'
#' @param spline_params `list` | `NULL`: Parameters for spline functions.
#' Must contain:
#'   - `spline_type`: `character(1)` `"n"` for natural cubic or `"b"` for
#'     B-splines.
#'   - `dof`: `integer(1)` Degrees of freedom. If set to `0`,
#'     `SplineOmics` automatically determines the best value using
#'     leave-one-out cross-validation.
#'   - `degree`: `integer(1)` Degree of the spline (B-splines only).
#'
#' @param padjust_method `character(1)`: Method for p-value adjustment. One of
#' `"none"`, `"BH"`, `"BY"`, `"holm"`, `"bonferroni"`, `"hochberg"`, or
#' `"hommel"`. Defaults to `"BH"` (Benjamini–Hochberg).
#'
#' @param bp_cfg `numeric()` | `NULL`: Named numeric vector specifying
#' parallelization settings, with expected names `"n_cores"` and
#' `"blas_threads"`. Controls the number of R worker processes (`n_cores`) and
#' BLAS/OpenBLAS threads per process (`blas_threads`). If `bp_cfg` is `NULL` or
#' missing, both default to `1`, disabling parallelization and avoiding thread
#' oversubscription.
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


#' Update the variables in a SplineOmics object
#'
#' @description
#' Updates a SplineOmics object by modifying existing fields or adding new ones.
#'
#' @param splineomics `SplineOmics`: A SplineOmics object to be updated.
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
#' @param x `SplineOmics`: A SplineOmics object created by the
#'  `create_splineomics` function.
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
        cat(
            "Annotation provided with",
            nrow(x$annotation),
            "entries.\n"
            )
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
