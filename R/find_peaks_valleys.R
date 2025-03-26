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
#' @param alpha A numeric significance threshold used to identify excursion 
#' points. Defaults to `0.05`.
#' @param padjust_method A character string specifying the method for multiple 
#' testing correction. Defaults to `"BH"` (Benjamini-Hochberg).
#'
#' @return A named list of ggplot objects, where each element corresponds to a 
#' feature with at least one detected peak or valley. Each plot shows expression 
#' profiles across timepoints, highlights excursion points in red, and annotates 
#' statistically significant excursions with significance stars.
#'
#' @details
#' A **peak** or **valley** is defined as a timepoint whose expression value is 
#' significantly different from both its immediate neighbors and deviates in 
#' the same direction — i.e., it is either significantly higher than both 
#' (a peak) or significantly lower than both (a valley).  
#' 
#' Statistically, this is tested using a compound contrast in limma:  
#' \deqn{(T - T_{prev}) + (T - T_{next}) = 2T - T_{prev} - T_{next}}  
#' This compound contrast has power only when the timepoint `T` is an outlier 
#' compared to both neighbors in the same direction. The resulting p-value is 
#' FDR-adjusted and compared to the `alpha` threshold.
#' 
#' - Performs internal input validation via `check_splineomics_elements()` and 
#'   `InputControl`.
#' - Detects local excursions using the `peaks_valleys_uit()` function.
#' - Displays the number of total excursion hits found.
#' - Generates plots using `plot_peaks_valleys()`, with each excursion point 
#'   evaluated for significance based on the provided `alpha`.
#'
#' @export
#' 
find_peaks_valleys <- function(
    splineomics,
    alpha = 0.05,
    padjust_method = "BH"
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
  meta <- splineomics[["meta"]]
  condition <- splineomics[["condition"]]
  meta_batch_column <- splineomics[["meta_batch_column"]]
  
  # Get all unique condition levels
  condition_levels <- unique(meta[[condition]])
  
  # Initialize output list
  results <- list()
  
  for (level in condition_levels) {
    # Select samples for this condition level
    selected_samples <- meta[[condition]] == level
    sub_data <- data[, selected_samples, drop = FALSE]
    sub_meta <- meta[selected_samples, , drop = FALSE]
    
    # Run peak/valley detection on the subset
    peak_valley_pvals <- peak_valley_test(
      data = sub_data,
      meta = sub_meta,
      alpha = alpha,
      padjust_method = padjust_method
    )

    labels <- classify_excursions(
      data = sub_data,
      meta = sub_meta,
      peak_valley_pvals = peak_valley_pvals,
      alpha = alpha
    )
    
    results[[as.character(level)]][["peak_valley_pvals"]] <- peak_valley_pvals
    results[[as.character(level)]][["labels"]] <- labels

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
      "\n\nDetected ",
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
      )
    )

    # Plot results
    plots <- plot_peaks_valleys(
      peak_valley_pvals = peak_valley_pvals,
      data = sub_data,
      meta = sub_meta,
      meta_batch_column = meta_batch_column,
      alpha = alpha
    )

    results[[as.character(level)]][["plots"]] <- plots
  }
  
  return(results)
}


# Level 1 function definitions -------------------------------------------------


#' Detect peaks/valleys in time-series omics using compound contrasts in limma
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
#' @param alpha A numeric value specifying the significance threshold for 
#' excursion detection. Defaults to `0.05`.
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
peak_valley_test <- function(
    data,
    meta,
    alpha = 0.05,
    padjust_method = "BH"
) {
  
  # Extract unique timepoints in sorted order
  unique_timepoints <- sort(unique(meta$Time))
  num_timepoints <- length(unique_timepoints)
  
  # Ensure there are at least 3 timepoints for meaningful comparisons
  if (num_timepoints < 3) {
    stop("Not enough timepoints for compound-contrast-based comparisons.")
  }
  
  # Create valid timepoint names for model matrix
  valid_timepoints <- make.names(as.character(unique_timepoints))
  
  # Design matrix for linear modeling
  time_factor <- factor(meta$Time, levels = unique_timepoints)
  design <- stats::model.matrix(~ 0 + time_factor)
  colnames(design) <- valid_timepoints
  
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


#' Classify Peaks, Valleys, and Cliffs from Compound Contrast P-values
#'
#' Assigns a label to each internal timepoint (T₂ to Tₙ₋₁) indicating 
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
    peak_valley_pvals,
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
  for (i in 1:nrow(data)) {
    for (t in 2:(num_timepoints - 1)) {
      p_val <- peak_valley_pvals[i, t - 1]
      
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
#'
#' @return A named list of ggplot objects, where each element corresponds to a 
#' feature with at least one detected excursion. Each plot displays the 
#' expression levels across timepoints, with replicates distinguished by shape.
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_shape_manual 
#'             scale_color_manual geom_text labs theme_minimal
#'
plot_peaks_valleys <- function(
    peak_valley_pvals,
    data,
    meta,
    meta_batch_column,
    alpha = 0.05  
) {
  
  peak_valley_flags <- ifelse(
    is.na(peak_valley_pvals), 
    0, 
    ifelse(peak_valley_pvals < alpha, 1, 0)
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
      Expression = as.numeric(protein_data), 
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
      # T₂ → [1], T₃ → [2], ...
      flag_index <- t - 1
      
      if (excursion_flags[flag_index] == 1) {
        p_val <- peak_valley_pvals[protein_index, flag_index]
        
        if (p_val < alpha) {
          stars <- get_stars(p_val, alpha)
          
          max_expr <- max(
            plot_data$Expression[plot_data$Time == timepoint],
            na.rm = TRUE
          )
          
          sig_df <- rbind(
            sig_df,
            data.frame(
              Time = timepoint,
              Label = stars,
              y_pos = max_expr + 0.01 * max_expr
            )
          )
        }
      }
    }

    p <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        x = Time,
        y = Expression
        )
    ) +
      ggplot2::geom_point(
        ggplot2::aes(
          shape = Replicate,
          color = Point_Type
          ),
        size = 3,
        stroke = 1.2,
        fill = "white"
      ) +
      ggplot2::scale_shape_manual(values = symbols) +
      ggplot2::scale_color_manual(values = c(
        "Normal" = "grey40",
        "Excursion" = "red"
        )) +
      ggplot2::geom_text(
        data = sig_df,
        ggplot2::aes(
          x = Time, 
          y = y_pos,
          label = Label
          ),
        size = 5,
        hjust = 0.5
      ) +
      ggplot2::labs(
        title = paste(
          "Feature:",
          feature_name
          ),
        x = "Time",
        y = "Feature Value",
        color = "Pattern",
        shape = "Replicate"
      ) +
      ggplot2::theme_minimal()
    
    plots[[feature_name]] <- p
  }
  
  return(plots)
}


# Level 2 function definitions -------------------------------------------------


#' Classify Cliff Pattern as Top or Bottom
#'
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
    return(ifelse(prev_change < 0, "b", "t"))  # drop into point → bottom
  } else {
    return(ifelse(next_change < 0, "t", "b"))  # drop after point → top
  }
}
