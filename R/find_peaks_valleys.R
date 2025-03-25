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
  all_plots <- list()
  
  for (cond_level in condition_levels) {
    # Select samples for this condition level
    selected_samples <- meta[[condition]] == cond_level
    sub_data <- data[, selected_samples, drop = FALSE]
    sub_meta <- meta[selected_samples, , drop = FALSE]
    
    # Run peak/valley detection on the subset
    uit_output <- peaks_valleys_uit(
      data = sub_data,
      meta = sub_meta,
      alpha = alpha,
      padjust_method = padjust_method
    )
    
    message(
      "Found ",
      sum(uit_output$results_df[, -1]),
      " peak/valley hits for condition level: ",
      cond_level
    )
    
    # Plot results
    plots <- plot_peaks_valleys(
      uit_output = uit_output,
      data = sub_data,
      meta = sub_meta,
      meta_batch_column = meta_batch_column,
      alpha = alpha
    )
    
    all_plots[[as.character(cond_level)]] <- plots
  }
  
  return(all_plots)
}


# Level 1 function definitions -------------------------------------------------


#' Detect Excursions in Time-Series Omics Data Using Union-Intersection Testing
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
peaks_valleys_uit <- function(
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
    stop("Not enough timepoints for UIT-based comparisons.")
  }
  
  # Create valid timepoint names for model matrix
  valid_timepoints <- make.names(as.character(unique_timepoints))
  
  # Design matrix for linear modeling
  time_factor <- factor(
    meta$Time,
    levels = unique_timepoints
  )
  design <- stats::model.matrix(~ 0 + time_factor)
  colnames(design) <- valid_timepoints

  # Fit limma model
  fit <- limma::lmFit(
    data,
    design
    )
  
  # Initialize matrix to store pairwise p-values
  pairwise_pvals <- matrix(
    NA,
    nrow = nrow(data),
    ncol = num_timepoints - 2
  )
  rownames(pairwise_pvals) <- rownames(data)
  colnames(pairwise_pvals) <- paste0(
    unique_timepoints[-c(1, num_timepoints)],
    "_vs_neighbors"
  )
  
  # Compute rolling-window UIT-based comparisons
  for (t in 2:(num_timepoints - 1)) {  # Rolling window, skipping first & last
    
    contrast_name <- paste0(
      "(", valid_timepoints[t], "-", valid_timepoints[t - 1], ") + ",
      "(", valid_timepoints[t], "-", valid_timepoints[t + 1], ")"
    )
    
    contrast_matrix <- limma::makeContrasts(
      contrasts = contrast_name,
      levels = design
    )
    
    fit2 <- limma::contrasts.fit(fit, contrast_matrix)
    fit2 <- limma::eBayes(fit2)
    
    pairwise_pvals[, t - 1] <- limma::topTable(
      fit2,
      coef = 1,
      number = Inf,
      adjust.method = padjust_method,
      sort.by = "none"
    )[, "adj.P.Val"]
  }
  
  # Initialize excursion matrix
  excursion_matrix <- matrix(
    0,
    nrow = nrow(data),
    ncol = num_timepoints
  )
  rownames(excursion_matrix) <- rownames(data)
  colnames(excursion_matrix) <- unique_timepoints
  
  # Detect excursions based on UIT
  for (i in 1:nrow(data)) {
    for (t in 2:(num_timepoints - 1)) {
      
      p_uit <- pairwise_pvals[i, t - 1]  # UIT p-value for T2 vs. neighbors
      
      prev_mean <- mean(data[i, which(meta$Time == unique_timepoints[t - 1])])
      curr_mean <- mean(data[i, which(meta$Time == unique_timepoints[t])])
      next_mean <- mean(data[i, which(meta$Time == unique_timepoints[t + 1])])
      
      prev_change <- curr_mean - prev_mean
      next_change <- next_mean - curr_mean
      
      # Excursion condition: UIT significant & T2 is strictly higher or lower
      if (!is.na(p_uit) && p_uit < alpha) {
        if (!is.na(prev_change) && !is.na(next_change)) {
          if ((prev_change > 0 && next_change < 0) ||
              (prev_change < 0 && next_change > 0)) {
            excursion_matrix[i, t] <- 1
          }
        }
      }
    }
  }
  
  results_df <- data.frame(
    feature_nr = rownames(excursion_matrix),
    excursion_matrix
  )
  
  return(list(
    results_df = results_df,
    pairwise_pvals = pairwise_pvals
  ))
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

plot_peaks_valleys <- function(
    uit_output,
    data,
    meta,
    meta_batch_column,
    alpha = 0.05  
) {
  
  results_df <- uit_output$results_df
  pairwise_pvals <- uit_output$pairwise_pvals
  
  unique_timepoints <- sort(unique(meta$Time))
  num_timepoints <- length(unique_timepoints)
  
  unique_replicates <- unique(meta[[meta_batch_column]])
  num_replicates <- length(unique_replicates)
  
  symbols_available <- c(21, 22, 23, 24, 25, 7, 8, 10, 12)  
  symbols <- symbols_available[1:num_replicates]  
  
  plots <- list()
  
  for (protein_index in which(rowSums(results_df[, -1]) > 0)) {
    
    feature_name <- results_df$feature_nr[protein_index]
    protein_data <- data[feature_name, ]
    
    plot_data <- data.frame(
      Time = meta$Time, 
      Expression = as.numeric(protein_data), 
      Replicate = as.factor(meta[[meta_batch_column]])
    )  
    
    excursion_flags <- as.numeric(results_df[protein_index, -1])
    timepoints_numeric <- sort(unique(meta$Time))
    plot_data$Excursion <- excursion_flags[match(
      plot_data$Time,
      timepoints_numeric
    )]
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
    
    for (t in 2:(num_timepoints - 1)) {
      if (excursion_flags[t] == 1) {
        p_val <- pairwise_pvals[protein_index, t - 1]
        if (p_val < alpha) {
          stars <- get_stars(
            p_val,
            alpha
            )
          
          # Compute the max expression at this excursion timepoint
          max_expr <- max(
            plot_data$Expression[plot_data$Time == unique_timepoints[t]],
            na.rm = TRUE
            )
          
          sig_df <- rbind(
            sig_df,
            data.frame(
              Time = unique_timepoints[t],
              Label = stars,
              y_pos = max_expr + 0.3
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