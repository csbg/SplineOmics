#' Find Peaks and Valleys in Time-Series Omics Data
#'
#' @description
#' Identifies significant local peaks or valleys (excursions) in time-series 
#' omics data using a Union-Intersection Test (UIT)-based approach. This 
#' function wraps the detection and plotting steps, returning visualizations of 
#' all features with at least one excursion.
#'
#' @param splineomics A list containing the preprocessed time-series input data. 
#' Must include named elements `"data"` (a numeric matrix), `"meta"` (a metadata 
#' data frame with a `"Time"` column), `"meta_batch_column"` (name of the column 
#' in `meta` identifying replicates or batches), and `"padjust_method"` (a 
#' string specifying the method for p-value adjustment).
#' @param alphas A single numeric value or a named list of numeric thresholds 
#' used to identify significant excursion points. If a single value is provided
#' (either as a numeric scalar or a list of length 1), the same threshold is 
#' applied to all condition levels. If a named list is provided, it must contain
#' one numeric 
#' value per condition level, with names matching the condition levels exactly. 
#' This input is normalized internally to ensure consistent per-level access.
#' @param padjust_method A character string specifying the method for multiple 
#' testing correction. Defaults to `"BH"` (Benjamini-Hochberg).
#' @param support Minimum amount of non-NA values in each timepoint that 
#' influenced a PVC-test result. For example, when the timepoints are 10, 15,
#' 20, etc. and support = 1, then for timepoint 15 for a given feature, the 
#' timepoints 10, 15, and 20 each must have at least 1 non-NA value. If one or
#' more of those timepoints for that feature don't meet this criterium, then the
#' p-value for that timepoint 15 and feature is set to NA.
#' @param plot_info List containing the elements y_axis_label (string),
#'                  time_unit (string), treatment_labels (character vector),
#'                  treatment_timepoints (integer vector). All can also be NA.
#'                  This list is used to add this info to the spline plots.
#'                  time_unit is used to label the x-axis, and treatment_labels
#'                  and -timepoints are used to create vertical dashed lines,
#'                  indicating the positions of the treatments (such as
#'                  feeding, temperature shift, etc.).
#' @param report_dir Character string specifying the directory path where the
#' HTML report and any other output files should be saved.
#'
#' @return A named list of ggplot objects, where each element corresponds to a 
#' feature with at least one detected peak or valley. Each plot shows expression 
#' profiles across timepoints, highlights excursion points in red, and annotates 
#' statistically significant excursions with significance stars.
#'
#' @details
#' A peak or valley is defined as a timepoint whose expression value is 
#' significantly different from both its immediate neighbors and deviates in 
#' the same direction - i.e., it is either significantly higher than both 
#' (a peak) or significantly lower than both (a valley).  
#' 
#' Statistically, this is tested using a compound contrast in limma:  
#' (T - T_prev) + (T - T_next) = 2T - T_prev - T_next
#' This compound contrast has power only when the timepoint `T` is an outlier 
#' compared to both neighbors in the same direction. The resulting p-value is 
#' FDR-adjusted and compared to the `alpha` threshold.
#' 
#' - Performs internal input validation via `check_splineomics_elements()` and 
#'   `InputControl`.
#' - Detects local excursions using the `pvc_test()` function.
#' - Displays the number of total excursion hits found.
#' - Generates plots using `plot_pvc()`, with each excursion point 
#'   evaluated for significance based on the provided `alpha`.
#'
#' @export
#' 
find_pvc <- function(
    splineomics,
    alphas = 0.05,
    padjust_method = "BH",
    support = 1,
    plot_info = list(
      y_axis_label = "Value",
      time_unit = "min",
      treatment_labels = NA,
      treatment_timepoints = NA
    ),
    report_dir = here::here()
) {
  
  check_splineomics_elements(
    splineomics = splineomics,
    func_type = "find_peaks_valleys"
  )
  
  args <- lapply(
    as.list(match.call()[-1]),
    eval,
    parent.frame()
  )
  
  check_null_elements(args)
  input_control <- InputControl$new(args)
  input_control$auto_validate()
  
  data <- splineomics[["data"]]
  annotation <- splineomics[["annotation"]]
  meta <- splineomics[["meta"]]
  condition <- splineomics[["condition"]]
  meta_batch_column <- splineomics[["meta_batch_column"]]
  report_info <- splineomics[["report_info"]]
  feature_name_columns <- splineomics[["feature_name_columns"]]
  
  # Get all unique condition levels
  condition_levels <- unique(meta[[condition]])
  alphas <- normalize_alphas(
    alphas = alphas,
    condition_levels = condition_levels
  )
  
  # Initialize output list
  results <- list()
  
  for (level_index in seq_along(condition_levels)) {
    level <- condition_levels[[level_index]]
    
    # Select samples for this condition level
    selected_samples <- meta[[condition]] == level
    sub_data <- data[, selected_samples, drop = FALSE]
    sub_meta <- meta[selected_samples, , drop = FALSE]
    
    alpha <- alphas[[level]]
    
    # Run peak/valley detection on the subset
    pvc_pvals <- pvc_test(
      data = sub_data,
      meta = sub_meta,
      padjust_method = padjust_method
    )
    
    pvc_pvals <- filter_pvals_by_support(
      pvc_pvals,
      sub_data,
      sub_meta,
      support
    )
    
    labels <- classify_excursions(
      data = sub_data,
      meta = sub_meta,
      pvc_pvals = pvc_pvals,
      alpha = alpha
    )
    
    # Count each label type per column (i.e., per timepoint)
    pattern_counts <- apply(
      labels,
      2,
      function(col) table(factor(col, levels = c("p", "v", "b", "t")))
    )
    
    # Convert to data frame for easier aggregation
    pattern_df <- as.data.frame(pattern_counts)
    pattern_summary <- rowSums(pattern_df)
    total_hits <- sum(pattern_summary)
    
    # Compose message
    message(
      "\nDetected ",
      total_hits,
      " total pattern hits for condition level: ",
      level,
      "\n\n",
      "Summary by pattern type:\n",
      paste(
        names(pattern_summary),
        pattern_summary,
        sep = ": ", 
        collapse = ", "
      ), 
      "\n\n",
      "Breakdown by timepoint:\n",
      paste(
        colnames(pattern_df),
        apply(
          pattern_df,
          2,
          function(x) paste(
            names(x),
            x,
            sep = "=",
            collapse = "; "
          )),
        sep = ": ",
        collapse = "\n"
      ),
      "\n"
    )
    
    # Plot results
    plots <- plot_pvc(
      pvc_pvals = pvc_pvals,
      data = sub_data,
      meta = sub_meta,
      meta_batch_column = meta_batch_column,
      alpha = alpha,
      plot_info = plot_info,
      level = level
    )
    
    results[[as.character(level)]][["plots"]] <- plots
    results[[as.character(level)]][["pvc_adj_pvals"]] <- pvc_pvals
    results[[as.character(level)]][["pvc_pattern_summary"]] <- pattern_df
  }
  
  # This info is passed like this so that it can be written in the HTML report
  level_headers_info <- list(
    pvc_settings = list(
      adj_value_thresh = alphas,
      padjust_method = padjust_method,
      support = support
    )
  )
  
  generate_report_html(
    plots = results,
    # required, but statically handled downstream for the pvc HTML report.
    plots_sizes = NULL,   
    report_info = report_info,
    data = bind_data_with_annotation(data, annotation),
    meta = meta,
    level_headers_info = level_headers_info,
    report_type = "find_pvc",
    feature_name_columns = feature_name_columns,
    filename = "pvc_report",
    report_dir = report_dir
  )
  
  print_info_message(
    message_prefix = "PVC report generation",
    report_dir = report_dir
  )
  
  return(results)
}


# Level 1 function definitions -------------------------------------------------


#' Normalize alpha thresholds for multiple condition levels
#'
#' @noRd
#' 
#' @description
#' This helper function standardizes the `alphas` input to ensure it can be 
#' consistently used across multiple condition levels. The input can either be 
#' a single numeric value, which will be replicated and named for each condition 
#' level, or a named list containing a numeric threshold for each condition 
#' level.
#'
#' @param alphas Either a single numeric value (e.g., `0.05`) or a named list of 
#' numeric values. If a single value is provided, it is applied to all 
#' `condition_levels`. If a list is provided, it must be named, with each name 
#' corresponding to a value in `condition_levels`.
#' 
#' @param condition_levels A character vector of condition levels for which 
#' alpha thresholds are required. These must match the names in the `alphas` 
#' list if a list is provided.
#'
#' @return A named list of alpha thresholds, one per condition level.
#'
normalize_alphas <- function(
    alphas,
    condition_levels
) {
  
  if (is.numeric(alphas) && length(alphas) == 1) {
    alphas <- setNames(as.list(rep(
      alphas,
      length(condition_levels))),
      condition_levels
    )
  } else if (is.list(alphas)) {
    if (!all(condition_levels %in% names(alphas))) {
      stop_call_false("Alpha list must be named for all condition levels.")
    }
  } else {
    stop_call_false("Invalid alphas: must be a single numeric or a named list.")
  }
  return(alphas)
}


#' Detect peaks/valleys in time-series omics using compound contrasts in limma
#'
#' @noRd
#'
#' @description
#' This function identifies excursions in time-series omics data using a 
#' rolling window Union-Intersection Test (UIT). A timepoint is flagged as 
#' an excursion if it is significantly different from both its adjacent 
#' neighbors in the same direction (higher or lower). The significance is 
#' determined using limma's moderated t-test with FDR correction.
#'
#' @param data A numeric matrix, where rows correspond to features (e.g., genes, 
#' proteins, metabolites) and columns correspond to samples.
#' @param meta A data frame containing metadata for the samples. Must include 
#' a column named `"Time"` that specifies the timepoint for each sample.
#' @param padjust_method Method to correct the p-values for multiple hypothesis
#'                       testing.
#'
#' @return A list containing two elements:
#' \describe{
#'   \item{`results_df`}{A data frame with binary indicators for excursions 
#'   at each timepoint for each feature. Rows correspond to features, columns 
#'   correspond to timepoints, and values are `1` (excursion) or `0` 
#'   (no excursion).}
#'   \item{`pairwise_pvals`}{A matrix of one-tailed adjusted p-values from 
#'   limma's moderated t-test for pairwise comparisons between consecutive 
#'   timepoints.}
#' }
#'
#' @importFrom limma lmFit makeContrasts contrasts.fit eBayes topTable
#' @importFrom stats model.matrix
#'
pvc_test <- function(
    data,
    meta,
    padjust_method = "BH"
) {
  
  # Extract unique timepoints in sorted order
  unique_timepoints <- sort(unique(meta$Time))
  num_timepoints <- length(unique_timepoints)
  
  # Ensure there are at least 3 timepoints for meaningful comparisons
  if (num_timepoints < 3) {
    stop("Not enough timepoints for compound-contrast-based comparisons.")
  }
  
  batch_effects <- c("Reactor")
  valid_timepoints <- make.names(as.character(unique_timepoints))
  # batch_effects <- NULL
  
  if (is.null(batch_effects)) {
    time <- factor(meta$Time, levels = unique_timepoints)
    design <- stats::model.matrix(~ 0 + time)
    colnames(design) <- valid_timepoints
  } else{
    meta$Time <- factor(meta$Time, levels = unique_timepoints)
    # Dynamically convert all batch effect columns to factors
    for (batch_var in batch_effects) {
      if (!batch_var %in% colnames(meta)) {
        stop(paste("Batch variable", batch_var, "not found in meta"))
      }
      meta[[batch_var]] <- factor(meta[[batch_var]])
    }

    # Build formula string: e.g. "~ 0 + Time + Reactor + HPLC_column"
    batch_formula_part <- paste(batch_effects, collapse = " + ")
    design_formula_str <- paste("~ 0 + Time +", batch_formula_part)
    design_formula <- as.formula(design_formula_str)

    # Create design matrix
    design <- model.matrix(design_formula, data = meta)
    colnames(design)[1:num_timepoints] <- valid_timepoints
  }
  
  # Fit base limma model
  fit <- limma::lmFit(data, design)
  
  # Initialize matrix to store compound contrast p-values
  peak_valley_pvals <- matrix(
    NA,
    nrow = nrow(data),
    ncol = num_timepoints - 2
  )
  rownames(peak_valley_pvals) <- rownames(data)
  colnames(peak_valley_pvals) <- paste0(
    unique_timepoints[-c(1, num_timepoints)],
    "_vs_neighbors"
  )
  
  # Loop over internal timepoints to compute compound contrast p-values
  for (t in 2:(num_timepoints - 1)) {
    compound_contrast <- paste0(
      "(", valid_timepoints[t], "-", valid_timepoints[t - 1], ") + ",
      "(", valid_timepoints[t], "-", valid_timepoints[t + 1], ")"
    )
    
    contrast_matrix <- limma::makeContrasts(
      contrasts = compound_contrast,
      levels = design
    )
    
    fit2 <- limma::contrasts.fit(fit, contrast_matrix)
    fit2 <- limma::eBayes(fit2)
    
    peak_valley_pvals[, t - 1] <- limma::topTable(
      fit2,
      coef = 1,
      number = Inf,
      adjust.method = padjust_method,
      sort.by = "none"
    )[,"adj.P.Val"]
  }
  
  return(peak_valley_pvals)
}


#' Filter p-values Based on Support Across Neighboring Timepoints
#'
#' @noRd
#'
#' @description
#' This function filters a p-value matrix by evaluating the data support
#' for each feature (row) and timepoint (column). A p-value is considered
#' valid only if the corresponding feature has at least a specified number
#' of non-NA values (`support`) across the current, previous, and next
#' timepoints.
#'
#' @details
#' Timepoint information is taken from the `Time` column in `sub_meta`,
#' which must be numeric. For each column in `pvc_pvals`, the function:
#'
#' * Extracts the timepoint from the column name (before the first underscore).
#' * Identifies the previous and next unique timepoints from `sub_meta$Time`.
#' * Finds all columns in `sub_data` associated with each of these timepoints.
#' * Checks if the feature has at least `support` non-NA values in each group.
#'
#' If any of the three groups (previous, current, next) fails the support
#' threshold, the corresponding p-value is replaced with `1`.
#'
#' Timepoints at the beginning or end of the time series (i.e., with no
#' previous or next timepoint) are treated as invalid and all their p-values
#' are set to `1`.
#'
#' @param pvc_pvals A numeric matrix of p-values. Each column corresponds to
#' a specific timepoint (with replicates), and each row to a feature.
#' @param sub_data A numeric matrix of feature data. Rows are features and
#' columns are samples. Column order must align with rows in `sub_meta`.
#' @param sub_meta A `data.frame` containing metadata for `sub_data` columns.
#' Must include a numeric `Time` column indicating sample timepoints.
#' @param support A non-negative integer specifying the minimum number of
#' non-NA values required for each of the three timepoint groups.
#'
#' @return A modified version of `pvc_pvals` where p-values are replaced with
#' `1` if support criteria are not met.
#' 
filter_pvals_by_support <- function(
    pvc_pvals,
    sub_data,
    sub_meta,
    support
) {
  
  # Ensure Time is numeric and sorted
  time_values <- sub_meta$Time
  if (!is.numeric(time_values)) {
    stop("sub_meta$Time must be numeric for sorting and comparison.")
  }
  
  # Get unique sorted timepoints
  unique_times <- sort(unique(time_values))
  
  # Create a lookup for each timepoint -> column indices in sub_data
  time_to_indices <- lapply(unique_times, function(t) which(time_values == t))
  names(time_to_indices) <- as.character(unique_times)
  
  # Initialize result
  filtered_pvals <- pvc_pvals
  
  # Loop over each column in pvc_pvals (each timepoint)
  for (col_idx in seq_len(ncol(pvc_pvals))) {
    # Extract timepoint t from column name (before first underscore)
    col_name <- colnames(pvc_pvals)[col_idx]
    t <- as.numeric(sub("_.*", "", col_name))
    if (is.na(t)) next  # skip if invalid timepoint
    
    # Identify previous and next timepoints
    t_pos <- which(unique_times == t)
    if (length(t_pos) == 0) next  # timepoint not found
    
    # Skip if t is the first or last (no prev/next)
    if (t_pos == 1 || t_pos == length(unique_times)) {
      filtered_pvals[, col_idx] <- 1
      next
    }
    
    t_prev <- unique_times[t_pos - 1]
    t_next <- unique_times[t_pos + 1]
    
    # Get column indices in sub_data for each group
    idx_prev <- time_to_indices[[as.character(t_prev)]]
    idx_curr <- time_to_indices[[as.character(t)]]
    idx_next <- time_to_indices[[as.character(t_next)]]
    
    # Combine data per feature (i.e., per row)
    for (row_idx in seq_len(nrow(pvc_pvals))) {
      vals_prev <- sub_data[row_idx, idx_prev]
      vals_curr <- sub_data[row_idx, idx_curr]
      vals_next <- sub_data[row_idx, idx_next]
      
      non_na_prev <- sum(!is.na(vals_prev))
      non_na_curr <- sum(!is.na(vals_curr))
      non_na_next <- sum(!is.na(vals_next))
      
      if (non_na_prev < support 
          || non_na_curr < support 
          || non_na_next < support) {
        filtered_pvals[row_idx, col_idx] <- 1
      }
    }
  }
  
  return(filtered_pvals)
}


#' Classify Peaks, Valleys, and Cliffs from Compound Contrast P-values
#'
#' @noRd
#'
#' @description
#' Assigns a label to each internal timepoint (T2 to Tn-1) indicating
#' whether it is a peak (\code{"p"}), valley (\code{"v"}), top of a 
#' cliff (\code{"t"}), or bottom of a cliff (\code{"b"}). Timepoints 
#' that are not statistically significant or cannot be classified are 
#' labeled as \code{"."}.
#'
#' The classification is based on compound contrast p-values from the 
#' \code{peak_valley_test()} function, as well as the relative symmetry 
#' of changes between adjacent timepoints.
#'
#' @param data A numeric matrix or data frame of expression values, with 
#' rows as features and columns matching \code{meta$Time}.
#' @param meta A data frame containing metadata with at least a column 
#' \code{Time} indicating timepoint assignment.
#' @param peak_valley_pvals A matrix of adjusted p-values from a compound 
#' contrast test, with rows matching \code{data} and columns corresponding 
#' to internal timepoints (i.e., excluding first and last timepoints).
#' @param alpha Significance threshold for p-value filtering. Only 
#' timepoints with \code{p < alpha} are considered for classification.
#' @param symmetry_ratio Minimum ratio of the smaller to larger slope 
#' required to consider a pattern symmetric enough to be a peak or valley 
#' (default: 0.3). Values below this threshold are classified as cliffs.
#'
#' @return A data frame with the same number of rows as \code{data} and 
#' one column per timepoint. Each entry is a character label: \code{"p"}, 
#' \code{"v"}, \code{"t"}, \code{"b"}, or \code{"."} for no signal.
#' 
classify_excursions <- function(
    data,
    meta,
    pvc_pvals,
    alpha = 0.05,
    symmetry_ratio = 0.3
) {
  
  # Extract sorted timepoints
  unique_timepoints <- sort(unique(meta$Time))
  num_timepoints <- length(unique_timepoints)
  
  # Initialize label matrix with "."
  peak_valley_labels <- matrix(
    ".", 
    nrow = nrow(data),
    ncol = num_timepoints
  )
  rownames(peak_valley_labels) <- rownames(data)
  colnames(peak_valley_labels) <- unique_timepoints
  
  # Classification logic for each triplet
  for (i in seq_len(nrow(data))) {
    for (t in 2:(num_timepoints - 1)) {
      p_val <- pvc_pvals[i, t - 1]
      
      # Skip if NA or not significant
      if (is.na(p_val) || p_val >= alpha) next
      
      # Extract timepoint means
      prev_mean <- mean(
        data[i, which(meta$Time == unique_timepoints[t - 1])],
        na.rm = TRUE
      )
      curr_mean <- mean(
        data[i, which(meta$Time == unique_timepoints[t])],
        na.rm = TRUE
      )
      next_mean <- mean(
        data[i, which(meta$Time == unique_timepoints[t + 1])],
        na.rm = TRUE
      )
      
      prev_change <- curr_mean - prev_mean
      next_change <- curr_mean - next_mean
      
      if (!is.na(prev_change) && !is.na(next_change)) {
        same_direction <- sign(prev_change) == sign(next_change)
        if (same_direction) {
          abs_prev <- abs(prev_change)
          abs_next <- abs(next_change)
          smaller <- min(abs_prev, abs_next)
          larger  <- max(abs_prev, abs_next)
          
          if (smaller / larger >= symmetry_ratio) {
            if (prev_change > 0) {
              peak_valley_labels[i, t] <- "p"
            } else {
              peak_valley_labels[i, t] <- "v"
            }
          } else {
            peak_valley_labels[i, t] <- classify_cliff(
              prev_change,
              next_change
            )
          }
        } else {
          peak_valley_labels[i, t] <- classify_cliff(
            prev_change,
            next_change
          )
        }
      }
    }
  }
  
  # Return as data.frame for compatibility
  results_df <- data.frame(
    peak_valley_labels,
    check.names = FALSE
  )
  
  return(results_df)
}


#' Plot Peaks and Valleys in Time-Series Omics Data
#'
#' @noRd
#'
#' @description
#' This function generates scatter plots for features that exhibit significant 
#' local peaks or valleys (excursions) in time-series omics data. Excursion 
#' points are highlighted in red, while normal points remain grey. If the 
#' excursion is statistically significant (based on the compound contrast test), 
#' a significance star is shown directly above the excursion point.
#'
#' @param results A list returned from `detect_excursions()`, containing 
#' `results_df` (excursion matrix) and `pairwise_pvals` (one-tailed p-values for 
#' the excursion contrast).
#' @param data A numeric matrix, where rows correspond to features (e.g., genes, 
#' proteins, metabolites) and columns correspond to samples.
#' @param meta A data frame containing metadata for the samples. Must include
#' a column named `"Time"` that specifies the timepoint for each sample.
#' @param meta_replicates_column A character string specifying the column name 
#' in `meta` that indicates biological replicates.
#' @param alpha A numeric value specifying the significance threshold 
#' for displaying stars above excursion points. Defaults to `0.05`.
#' @param plot_info List containing the elements y_axis_label (string),
#'                  time_unit (string), treatment_labels (character vector),
#'                  treatment_timepoints (integer vector). All can also be NA.
#'                  This list is used to add this info to the spline plots.
#'                  time_unit is used to label the x-axis, and treatment_labels
#'                  and -timepoints are used to create vertical dashed lines,
#'                  indicating the positions of the treatments (such as
#'                  feeding, temperature shift, etc.).
#'
#' @return A named list of ggplot objects, where each element corresponds to a 
#' feature with at least one detected excursion. Each plot displays the 
#' expression levels across timepoints, with replicates distinguished by shape.
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_shape_manual 
#'             scale_color_manual geom_text labs theme_minimal
#'
plot_pvc <- function(
    pvc_pvals,
    data,
    meta,
    meta_batch_column,
    alpha = 0.05,
    plot_info,
    level
) {
  
  peak_valley_flags <- ifelse(
    is.na(pvc_pvals), 
    0, 
    ifelse(pvc_pvals < alpha, 1, 0)
  )
  
  unique_timepoints <- sort(unique(meta$Time))
  num_timepoints <- length(unique_timepoints)
  
  unique_replicates <- unique(meta[[meta_batch_column]])
  num_replicates <- length(unique_replicates)
  
  symbols_available <- c(21, 22, 23, 24, 25, 7, 8, 10, 12)  
  symbols <- symbols_available[1:num_replicates]  
  
  plots <- list()
  
  for (protein_index in which(rowSums(peak_valley_flags) > 0)) {
    feature_name <- rownames(peak_valley_flags)[protein_index]
    protein_data <- data[feature_name, ]
    
    plot_data <- data.frame(
      Time = meta$Time, 
      Feature_value = as.numeric(protein_data), 
      Replicate = as.factor(meta[[meta_batch_column]])
    )  
    
    excursion_flags <- as.numeric(peak_valley_flags[protein_index, ])
    
    # Step 1: define internal timepoints
    internal_timepoints <- unique_timepoints[2:(num_timepoints - 1)]
    
    # Step 2: map excursion flags to only internal timepoints
    named_flags <- excursion_flags
    names(named_flags) <- as.character(internal_timepoints)
    
    # Step 3: assign only to matching timepoints
    plot_data$Excursion <- named_flags[as.character(plot_data$Time)]
    plot_data$Excursion[is.na(plot_data$Excursion)] <- 0
    
    
    plot_data$Point_Type <- factor(
      ifelse(
        plot_data$Excursion == 1,
        "Excursion",
        "Normal"
      ),
      levels = c(
        "Normal",
        "Excursion"
      )
    )
    
    # For plotting significance stars directly above excursion points
    sig_df <- data.frame(
      Time = numeric(0),
      Label = character(0),
      y_pos = numeric(0)
    )
    
    get_stars <- function(p, alpha) {
      if (p < alpha / 500) return("****")
      else if (p < alpha / 50) return("***")
      else if (p < alpha / 5) return("**")
      else if (p < alpha) return("*")
      else return("")
    }
    
    # Loop over internal timepoints
    for (t in 2:(num_timepoints - 1)) {
      timepoint <- unique_timepoints[t]
      
      # excursion_flags is aligned to internal timepoints:
      # T2 -> [1], T3 -> [2], ...
      flag_index <- t - 1
      
      if (excursion_flags[flag_index] == 1) {
        p_val <- pvc_pvals[protein_index, flag_index]
        
        if (p_val < alpha) {
          stars <- get_stars(p_val, alpha)
          
          max_value <- max(
            plot_data$Feature_value[plot_data$Time == timepoint],
            na.rm = TRUE
          )
          
          sig_df <- rbind(
            sig_df,
            data.frame(
              Time = timepoint,
              Label = stars,
              max_value = max_value,
              PValue = p_val
            )
          )
        }
      }
    }
    
    # Build the p-value string if any are significant
    if (nrow(sig_df) > 0) {
      pval_str <- paste0(
        "adj.p-val: ",
        paste0(
          "T=", sig_df$Time, " -> ", 
          formatC(sig_df$PValue, format = "fg", digits = 4),
          collapse = "; "
        )
      )
      plot_title <- paste(feature_name, "\n", pval_str)
    } else {
      plot_title <- feature_name
    }
    
    p <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        x = Time,
        y = Feature_value
      )
    ) +
      ggplot2::geom_point(
        ggplot2::aes(
          shape = Replicate,
          color = Point_Type
        ),
        size = 3,
        stroke = 1.2,
        fill = "white",
        na.rm = TRUE
      ) +
      ggplot2::scale_shape_manual(values = symbols) +
      ggplot2::geom_text(
        data = sig_df,
        ggplot2::aes(
          x = Time, 
          y = max_value,
          label = Label
        ),
        size = 5,
        hjust = 0.5,
        vjust = -0.5   # Positions the significance stars a bit above datapoint
      ) +
      ggplot2::labs(
        title = plot_title,
        x = paste0("Time [", plot_info[["time_unit"]], "]"),
        y = plot_info[["y_axis_label"]],
        color = "Timepoints",
        shape = "Replicates"
      ) +
      ggplot2::theme_minimal()
    
    y_min <- min(plot_data$Feature_value, na.rm = TRUE)
    y_max <- max(plot_data$Feature_value, na.rm = TRUE)
    y_extension <- (y_max - y_min) * 0.1
    y_pos <- y_max + y_extension
    
    result <- maybe_add_dashed_lines(   
      p = p,
      plot_info = plot_info,
      level = level,
      y_pos = y_pos
    )
    
    p <- result$p
    treatment_colors <- result$treatment_colors
    
    color_values <- c(
      "Normal" = "grey40",
      "Excursion" = "red",
      treatment_colors
    )
    
    p <- p + ggplot2::scale_color_manual(
      values = color_values,
      guide = ggplot2::guide_legend(
        override.aes = list(
          size = c(
            rep(1.5, 2),  # Point sizes for Normal and Excursion
            rep(0.5, length(treatment_colors))  # Line size for treatments
          )
        )
      )
    )
    
    plots[[feature_name]] <- p
  }
  
  return(plots)
}


#' Build full PVC HTML report
#' 
#' @noRd
#'
#' @description
#' This function assembles the complete HTML report for PVC analysis, including 
#' section headers, significance metadata, subplot titles, plots, and a table of 
#' contents. It processes a nested list of plots and writes the final report to 
#' the specified output file path.
#'
#' @param header_section A character string of HTML content that represents the 
#'   top-level header section of the report (e.g., title, description, metadata)
#' @param plots A named list of plot objects, each containing a `plots` field 
#'   (itself a named list of ggplot2 objects). The names define report sections.
#' @param level_headers_info A list of metadata associated with each report 
#'   section, including p-value thresholds and padjust methods, structured under 
#'   the `pvc_settings` key.
#' @param report_info A list containing metadata about the entire report (e.g., 
#'   parameters used, timestamps, version info) used in the report footer.
#' @param output_file_path File path (character) where the final HTML report 
#'   should be written. Defaults to `here::here()` if not provided.
#'
#' @return Used for side effects only. The function writes an HTML file to disk.
#'
build_pvc_report <- function(
    header_section,
    plots,
    level_headers_info,
    report_info,
    output_file_path = here::here()
) {
  
  html_content <- paste(
    header_section,
    "<!--TOC-->",
    sep = "\n"
  )
  
  toc <- create_toc()
  
  styles <- define_html_styles()
  section_header_style <- styles$section_header_style
  toc_style <- styles$toc_style
  
  current_header_index <- 1
  
  # Flatten the nested structure into a list of (header_name, plot) pairs
  flattened_plots <- list()
  for (header_name in names(plots)) {
    subplots <- plots[[header_name]]$plots
    subplot_names <- names(subplots)
    
    for (i in seq_along(subplots)) {
      flattened_plots[[length(flattened_plots) + 1]] <- list(
        header_name = header_name,
        subplot_name = subplot_names[i],
        plot = subplots[[i]]
      )
    }
  }
  
  # Clean up header info
  levels <- names(plots)
  n_plots_per_header <- vapply(plots, function(x) length(x$plots), integer(1))
  
  # Track current section and plot offset
  current_header_index <- 1
  plots_processed_in_section <- 0
  
  # Create progress bar
  pb <- create_progress_bar(flattened_plots)
  nr_of_hits <- length(flattened_plots)
  
  for (index in seq_along(flattened_plots)) {
    if (plots_processed_in_section == 0) {
      # Insert section header
      level <- levels[[current_header_index]]
      
      html_content <- paste(
        html_content,
        generate_section_header_block(
          level = levels[[current_header_index]],
          index = index,
          section_header_style = section_header_style,
          level_headers_info = level_headers_info,
          nr_of_hits = nr_of_hits
        ),
        sep = "\n"
      )
      
      hits_info <- "<p style='text-align: center; font-size: 30px;'>"
      html_content <- paste(
        html_content,
        hits_info,
        sep = "\n"
      )
      
      toc_entry <- sprintf(
        "<li style='%s'><a href='#section%d'>%s</a></li>",
        toc_style,
        index,
        level
      )
      toc <- paste(
        toc,
        toc_entry,
        sep = "\n"
      )
    }
    
    # Add subplot name as an HTML sub-header
    subplot_name <- flattened_plots[[index]]$subplot_name
    subplot_header <- sprintf(
      paste0(
        '<h4 style="text-align: center; font-size: 20px; ',
        'margin-top: 30px;">\n%s\n</h4>'
      ),
      subplot_name
    )
    
    # Inject it into the HTML content before the plot
    html_content <- paste(
      html_content,
      subplot_header,
      sep = "\n"
    )
    
    result <- process_plots(
      plots_element = flattened_plots[[index]]$plot,
      plots_size = 2,
      html_content = html_content,
      toc = toc,
      header_index = current_header_index
    )
    html_content <- result$html_content
    toc <- result$toc
    
    # Update section counters
    plots_processed_in_section <- plots_processed_in_section + 1
    if (plots_processed_in_section >= 
        n_plots_per_header[[current_header_index]]) {
      current_header_index <- current_header_index + 1
      plots_processed_in_section <- 0
    }
    
    pb$tick()
  }
  
  generate_and_write_html(
    toc = toc,
    html_content = html_content,
    report_info = report_info,
    output_file_path = output_file_path
  )
}


# Level 2 function definitions -------------------------------------------------


#' Classify Cliff Pattern as Top or Bottom
#'
#' @noRd
#'
#' @description
#' Determines whether a significant excursion that is not a peak or valley 
#' represents the top or bottom of a cliff, based on the direction and 
#' magnitude of changes before and after a timepoint.
#'
#' @param prev_change Numeric value representing the change from the previous 
#'                    timepoint to the current timepoint.
#' @param next_change Numeric value representing the change from the current 
#'                    timepoint to the next timepoint.
#'
#' @return A single character: \code{"t"} for top of a cliff (upward into or 
#' downward out of the timepoint),
#' or \code{"b"} for bottom of a cliff (downward into or upward out of the 
#' timepoint).
#'
#' @examples
#' classify_cliff(-2, 0.5)  # returns "b" (bottom of a cliff)
#' classify_cliff(1.5, -3)  # returns "t" (top of a cliff)
#'
classify_cliff <- function(
    prev_change,
    next_change
) {
  
  abs_prev <- abs(prev_change)
  abs_next <- abs(next_change)
  
  if (abs_prev > abs_next) {
    return(ifelse(prev_change < 0, "b", "t"))  # drop into point -> bottom
  } else {
    return(ifelse(next_change < 0, "t", "b"))  # drop after point -> top
  }
}


#' Generate HTML header block for a condition level section
#' 
#' @noRd
#'
#' @description
#' Constructs an HTML block that includes the section header and related 
#' metadata 
#' for a given condition level. This includes the adjusted p-value threshold, 
#' the p-adjustment method, and a visual legend for significance stars based on 
#' the threshold. The output is centered and styled with inline HTML for use in 
#' reports or R Markdown rendering.
#'
#' @param level A character string representing the current condition level name
#' @param index An integer index used to create a unique HTML section ID.
#' @param section_header_style A character string of inline CSS styles to apply 
#'   to the section header (`<h2>` element).
#' @param level_headers_info A list containing `pvc_settings`, with named 
#'   entries for `adj_value_thresh` (a named list of alpha values by level) and 
#'   `padjust_method` (a scalar string).
#'
#' @return A character string of HTML content to be included in the full HTML 
#'         output.
#'
generate_section_header_block <- function(
    level,
    index,
    section_header_style,
    level_headers_info,
    nr_of_hits
) {
  
  # Access per-level alpha threshold
  alpha_thresh <- 
    level_headers_info[["pvc_settings"]][["adj_value_thresh"]][[level]]
  
  # Access global padjust_method
  padjust_method <- level_headers_info[["pvc_settings"]][["padjust_method"]]
  support <- level_headers_info[["pvc_settings"]][["support"]]
  
  # Format alpha threshold smartly
  format_alpha <- function(x) {
    if (x >= 0.0001) {
      # remove trailing zeros
      sub("\\.?0+$", "", format(round(x, 4), nsmall = 0))  
    } else {
      format(
        x,
        scientific = TRUE,
        digits = 3
      )
    }
  }
  
  # Compute significance thresholds
  thresholds <- c(
    "*"    = alpha_thresh,
    "**"   = alpha_thresh / 5,
    "***"  = alpha_thresh / 50,
    "****" = alpha_thresh / 500
  )
  
  thresholds_formatted <- vapply(
    thresholds,
    format_alpha,
    character(1)
  )
  alpha_thresh_formatted <- format_alpha(alpha_thresh)
  
  # Generate section header
  section_header <- sprintf(
    "<h2 style='%s' id='section%d'>%s</h2>",
    section_header_style,
    index,
    level
  )
  
  # Generate info block
  info_block <- sprintf(
    "<div style='text-align: center; font-size: 24px; margin-bottom: 20px;'>
       <p><strong>Adjusted p-value threshold:</strong> %s<br>
       <strong>p-adjustment method:</strong> %s<br>
       <strong>Support (min nr non-NA in all timepoints of hit):</strong> %s</p>
       <strong>Number of hits:</strong> %s</p>
       <p style='margin-top: 15px; font-size: 22px;'>
         <strong>Significance stars:</strong><br>
         <span>*</span> : p &lt; %s<br>
         <span>**</span> : p &lt; %s<br>
         <span>***</span> : p &lt; %s<br>
         <span>****</span> : p &lt; %s
       </p>
       <hr style='border-top: 1px dashed #999; width: 100%%; margin-top: 15px;'>
     </div>",
    alpha_thresh_formatted,
    padjust_method,
    support,
    nr_of_hits,
    thresholds_formatted["*"],
    thresholds_formatted["**"],
    thresholds_formatted["***"],
    thresholds_formatted["****"]
  )
  
  # Combine and return the full HTML block
  return(paste(
    section_header,
    info_block,
    sep = "\n"
  ))
}