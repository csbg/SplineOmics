# Exported function: preprocess_rna_seq_data() ---------------------------------


#' Perform default preprocessing of raw RNA-seq counts
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
#' @param raw_counts A matrix of raw RNA-seq counts (genes as rows, samples as
#'  columns).
#' @param meta A dataframe containing the metadata for data.
#' @param spline_params Parameters for spline functions (optional). Must contain
#' the named elements spline_type, which must contain either the string "n" for
#' natural cubic splines, or "b", for B-splines, the named element degree in the
#' case of B-splines, that must contain only an integer, and the named element 
#' dof, specifying the degree of freedom, containing an integer and required
#' both for natural and B-splines.
#' @param design A design formula for the limma analysis, such as 
#' '~ 1 + Phase*X + Reactor'.
#' @param normalize_func An optional normalization function. If provided, this 
#' function will be used to normalize the `DGEList` object. If not provided,
#' TMM normalization (via `edgeR::calcNormFactors`) will be used by default.
#' Must take as
#' input the y of: y <- edgeR::DGEList(counts = raw_counts) and output the y 
#' with the normalized counts.
#' @return A `voom` object, which includes the log2-counts per million (logCPM)
#'  matrix and observation-specific weights.
#' 
#' @importFrom limma voom
#' 
#' @export
#'   
preprocess_rna_seq_data <- function(
    raw_counts,
    meta,
    spline_params,
    design,
    normalize_func = NULL
) {
  
  message("Preprocessing RNA-seq data (normalization + voom)...")
  
  # Check if edgeR is installed; if not, prompt the user
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    message("The 'edgeR' package is not installed.")
    
    # Prompt user for action
    repeat {
      user_input <- readline(
        prompt = 
          "What would you like to do?\n
        1: Automatically install edgeR\n
        2: Manually install edgeR\n
        3: Cancel\n
        Please enter 1, 2, or 3: "
      )
      
      if (user_input == "1") {
        # Try to install edgeR automatically from Bioconductor
        message("Attempting to install 'edgeR' automatically 
                from Bioconductor...")
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          utils::install.packages("BiocManager")
        }
        tryCatch(
          {
            BiocManager::install("edgeR", update = FALSE)
          },
          error = function(e) {
            stop(
              "Automatic installation of 'edgeR' failed.
            Please install it manually and try again.",
              call. = FALSE
            )
          }
        )
        break  # Exit the loop if installation is successful
      } else if (user_input == "2") {
        stop(
          "Please install 'edgeR' manually using 
         BiocManager::install('edgeR') and then re-run the function.",
          call. = FALSE
        )
      } else if (user_input == "3") {
        stop("Operation canceled by the user.", call. = FALSE)
      } else {
        message("Invalid input. Please enter 1, 2, or 3.")
      }
    }
  }
  
  design_matrix <- design2design_matrix(
    meta = meta,
    spline_params = spline_params,
    level_index = 1,
    design = design
  )
  
  # Step 1: Create DGEList object from raw counts
  y <- edgeR::DGEList(counts = raw_counts)
  
  # Step 2: Apply the normalization function (either user-provided or default)
  if (!is.null(normalize_func) && is.function(normalize_func)) {
    y <- normalize_func(y)   # user provided normalisation function
  } else {
    # Default: Normalize the counts using TMM normalization
    y <- edgeR::calcNormFactors(y)
  }
  
  # Step 3: Apply voom transformation to get logCPM values and weights
  voom_obj <- limma::voom(
    y,
    design_matrix
  )
  
  return(voom_obj)
}