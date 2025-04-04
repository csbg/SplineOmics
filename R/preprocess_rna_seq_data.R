#' preprocess_rna_seq_data()
#'
#' @description
#' The `preprocess_rna_seq_data()` function performs essential preprocessing
#' steps for raw RNA-seq counts. This includes creating a `DGEList` object,
#' normalizing the counts using the default TMM (Trimmed Mean of M-values)
#' normalization via the `edgeR::calcNormFactors` function, and applying the
#' `voom` transformation from the `limma` package to obtain log-transformed
#' counts per million (logCPM) with associated precision weights. If you
#' require a different normalization method, you can supply your own
#' custom normalization function.
#'
#' @param splineomics An S3 object of class `SplineOmics` that must contain the 
#' following elements:
#' \itemize{
#'   \item \code{data}: A matrix of the omics dataset, with feature names 
#'   optionally as row headers (genes as rows, samples as columns).
#'   \item \code{meta}: A dataframe containing metadata corresponding to the 
#'   \code{data}. The dataframe must include a 'Time' column and a column 
#'   specified by the \code{condition}.
#'   \item \code{design}: A character string representing the design formula 
#'   for the limma analysis (e.g., \code{'~ 1 + Phase*X + Reactor'}).
#'   \item \code{spline_params}: A list of spline parameters used in the 
#'   analysis. This can include:
#'     \itemize{
#'       \item \code{spline_type}: A character string specifying the type of 
#'       spline. Must be either \code{'n'} for natural cubic splines or 
#'       \code{'b'} for B-splines.
#'       \item \code{dof}: An integer specifying the degrees of freedom. 
#'       Required for both natural cubic splines and B-splines.
#'       \item \code{degree}: An integer specifying the degree of the spline 
#'       (for B-splines only).
#'       \item \code{knots}: Positions of the internal knots (for B-splines).
#'       \item \code{bknots}: Boundary knots (for B-splines).
#'     }
#'   \item \code{dream_params}: A named list or \code{NULL}. When not 
#'   \code{NULL}, it can contain:
#'     \itemize{
#'       \item \code{dof}: An integer greater than 1, specifying the degrees 
#'       of freedom for the dream topTable.
#'       \item \code{KenwardRoger}: A boolean indicating whether to use the 
#'       Kenward-Roger method.
#'     }
#' }
#' @param normalize_func An optional normalization function. If provided, this
#' function will be used to normalize the `DGEList` object. If not provided,
#' TMM normalization (via `edgeR::calcNormFactors`) will be used by default.
#' Must take as
#' input the y of: y <- edgeR::DGEList(counts = raw_counts) and output the y
#' with the normalized counts.
#' @return A `voom` object, which includes the log2-counts per million (logCPM)
#'  matrix and observation-specific weights.
#'
#' @seealso \code{\link[edgeR]{DGEList}}, \code{\link[edgeR]{calcNormFactors}}
#' @importFrom limma voom
#' @importFrom variancePartition voomWithDreamWeights
#'
#' @export
#'
preprocess_rna_seq_data <- function(
    splineomics,
    normalize_func = NULL
    ) {
  
  check_splineomics_elements(
    splineomics = splineomics,
    func_type = "preprocess_rna_seq_data"
  )
  
  args <- lapply(
    as.list(match.call()[-1]),
    eval,
    parent.frame()
  )
  
  check_null_elements(args)
  input_control <- InputControl$new(args)
  input_control$auto_validate()
  
  raw_counts <- splineomics[["data"]]
  meta <- splineomics[["meta"]]
  spline_params <- splineomics[["spline_params"]]
  design <- splineomics[["design"]]
  robust_fit <- splineomics[["robust_fit"]]

  # Because at first I enforced that X in the design formula stands for the time
  # and I heavily oriented my code towards that. But then I realised that it is
  # nonsense to encode the time as X, and now it is explicitly "Time" (because
  # meta must contain the exact name "Time" for this respective column).
  design <- gsub("Time", "X", design)  
  
  effects <- extract_effects(design)

  message("Preprocessing RNA-seq data (normalization + voom)...")

  # Check if edgeR is installed; if not, inform the user
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop_call_false(
      "The 'edgeR' package is not installed.\n",
      "Please install it manually (to your custom_lib_path) using the command", 
      " below and re-run the function:\n\n",
      "  BiocManager::install('edgeR')\n\n",
      "This is an optional dependency of the SplineOmics package, ",
      "only needed when working with RNA-seq data."
    )
  }
  
  # voomWithDreamWeights wants it like this
  colnames(raw_counts) <- rownames(meta)  
  
  # Step 1: Create DGEList object from raw counts
  y <- edgeR::DGEList(counts = raw_counts)

  # Step 2: Apply the normalization function (either user-provided or default)
  if (!is.null(normalize_func) && is.function(normalize_func)) {
    y <- normalize_func(y) # user provided normalisation function
  } else {
    # Default: Normalize the counts using TMM normalization
    y <- edgeR::calcNormFactors(y)
  }

  # Step 3: Create design matrix
  result <- design2design_matrix(
    meta = meta,
    spline_params = spline_params,
    level_index = 1,
    design = effects[["fixed_effects"]]
  )

  # Step 4: Apply voom transformation to get logCPM values and weights
  if (effects[["random_effects"]] != "") {
    if (!is.null(robust_fit) && robust_fit == TRUE) {
      message(
        "⚠️  robust_fit = TRUE is ignored for mixed model RNA-seq.\n",
        "voomWithDreamWeights already handles heteroscedasticity internally."
        )
    }
    
    voom_obj <- variancePartition::voomWithDreamWeights(
      counts = y,
      formula = stats::as.formula(design),
      data = result[["meta"]]
    )
  }
  else {
    design_matrix <- result[["design_matrix"]]
    
    if (is.null(robust_fit)) {     # means fallback to implicit handling
      # Run voom normally first
      voom_obj <- limma::voom(
        counts = y,
        design = design_matrix
      )
      
      violation <- check_homoscedasticity_violation(
        data = voom_obj$E,
        meta = meta,
        data_type = "rna-seq"
      )
      
      # Step 4: If any pair was violated, rerun with robust weights
      if (violation) {
        message(
        "Using voomWithQualityWeights() due to detected violation of the 
        assumption of homoscedasticity."
        )
        voom_obj <- limma::voomWithQualityWeights(
          counts = y,
          design = design_matrix
        )
      }
    } else if (robust_fit == TRUE) {
      voom_obj <- limma::voomWithQualityWeights(
        counts = y,
        design = design_matrix
      )
    } else {   # robust_fit == FALSE
      voom_obj <- limma::voom(
        counts = y,
        design = design_matrix
      )
    }
  }

  return(voom_obj)
}
