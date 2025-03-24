# Level 1 function definitions -------------------------------------------------


#' Detect Excursions in Time-Series Omics Data
#'
#' @description
#' This function identifies excursions in time-series omics data by performing
#' adjacency-based pairwise comparisons between consecutive timepoints using
#' limma's moderated t-test. Excursions are defined as timepoints that are
#' significantly different from both their adjacent neighbors, with their mean
#' being either higher or lower than both neighbors. All p-values are converted
#' to one-tailed values to enhance sensitivity.
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
#' @details
#' - The function first fits a linear model to the data using limma.
#' - It then computes pairwise comparisons between adjacent timepoints.
#' - All p-values are adjusted using the Benjamini-Hochberg method for false
#' discovery rate (FDR) correction and subsequently divided by `2` to
#' reflect a one-tailed test.
#' - Excursions are identified if a timepoint exhibits significant differences
#' from both its neighbors and is either higher or lower than both.
#'
#' @importFrom limma lmFit makeContrasts contrasts.fit eBayes topTable
#' @importFrom stats model.matrix
#'
detect_excursions <- function(
    data,
    meta,
    alpha = 0.05
) {

  # Extract unique timepoints in sorted order
  unique_timepoints <- sort(unique(meta$Time))
  num_timepoints <- length(unique_timepoints)

  # Ensure there are at least 3 timepoints for meaningful comparisons
  if (num_timepoints < 3) {
    stop("Not enough timepoints for pairwise comparisons.")
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

  # Compute pairwise comparisons
  pairwise_pvals <- matrix(
    NA,
    nrow = nrow(data),
    ncol = num_timepoints - 1
  )
  rownames(pairwise_pvals) <- rownames(data)
  colnames(pairwise_pvals) <- paste0(
    unique_timepoints[-1],
    "_vs_",
    unique_timepoints[-num_timepoints]
  )

  for (t in 1:(num_timepoints - 1)) {
    contrast_name <- paste(
      valid_timepoints[t + 1],
      "-", valid_timepoints[t]
    )
    contrast_matrix <- limma::makeContrasts(
      contrasts = contrast_name,
      levels = design
    )

    fit2 <- limma::contrasts.fit(
      fit,
      contrast_matrix
    )
    fit2 <- limma::eBayes(fit2)

    # Store adjusted p-values (two-tailed, but we'll convert them)
    pairwise_pvals[, t] <- limma::topTable(
      fit2,
      coef = 1,
      number = Inf,
      adjust.method = "fdr",
      sort.by = "none"
    )[, "adj.P.Val"] / 2  # Convert to one-tailed
  }

  # Initialize excursion matrix
  excursion_matrix <- matrix(
    0,
    nrow = nrow(data),
    ncol = num_timepoints
  )
  rownames(excursion_matrix) <- rownames(data)
  colnames(excursion_matrix) <- unique_timepoints

  # Detect excursions
  for (i in 1:nrow(data)) {
    for (t in 2:(num_timepoints - 1)) {

      p_prev <- pairwise_pvals[i, t - 1]  # T1 vs. T2 (now one-tailed)
      p_next <- pairwise_pvals[i, t]  # T2 vs. T3 (now one-tailed)

      prev_mean <- mean(data[i, which(meta$Time == unique_timepoints[t - 1])])
      curr_mean <- mean(data[i, which(meta$Time == unique_timepoints[t])])
      next_mean <- mean(data[i, which(meta$Time == unique_timepoints[t + 1])])

      prev_change <- curr_mean - prev_mean
      next_change <- next_mean - curr_mean

      # Excursion condition: Significant changes in both directions &
      # both higher or lower
      if (p_prev < alpha & p_next < alpha) {
        if ((prev_change > 0 & next_change < 0) |
            (prev_change < 0 & next_change > 0)) {
          excursion_matrix[i, t] <- 1
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







detect_patterns <- function(
    data,
    meta,
    patterns,
    alpha = 0.05,
    delta = 0.2 # Threshold for practical equivalence when pattern contains 0
) {
  
  # Extract unique timepoints
  unique_timepoints <- sort(unique(meta$Time))
  num_timepoints <- length(unique_timepoints)
  
  # Ensure at least 3 timepoints
  if (num_timepoints < 3) {
    stop("Not enough timepoints for meaningful pattern detection.")
  }
  
  # Create valid timepoint names
  valid_timepoints <- make.names(as.character(unique_timepoints))
  
  # Design matrix for limma
  time_factor <- factor(
    meta$Time,
    levels = unique_timepoints
    )
  design <- stats::model.matrix(~ 0 + time_factor)
  colnames(design) <- valid_timepoints  
  
  # Fit limma model
  fit <- limma::lmFit(data, design)
  
  # Store pairwise p-values & logFCs
  pairwise_pvals <- matrix(
    NA,
    nrow = nrow(data),
    ncol = num_timepoints - 1
    )
  pairwise_logFC <- matrix(
    NA,
    nrow = nrow(data),
    ncol = num_timepoints - 1
    )
  rownames(pairwise_pvals) <- rownames(data)
  colnames(pairwise_pvals) <- paste0(
    unique_timepoints[-1],
    "_vs_",
    unique_timepoints[-num_timepoints]
  )
  
  # Compute pairwise comparisons
  for (t in 1:(num_timepoints - 1)) {
    contrast_name <- paste(
      valid_timepoints[t + 1],
      "-", valid_timepoints[t]
      )
    contrast_matrix <- limma::makeContrasts(
      contrasts = contrast_name,
      levels = design
      )
    
    fit2 <- limma::contrasts.fit(
      fit,
      contrast_matrix
      )
    fit2 <- limma::eBayes(fit2)
    
    # Compute raw (two-tailed) p-values
    raw_pvals <- limma::topTable(
      fit2,
      coef = 1,
      number = Inf,
      adjust.method = "none",
      sort.by = "none"
    )
    
    # Convert to one-tailed p-values
    one_tailed_pvals <- raw_pvals[, "P.Value"] / 2  # Divide by 2 first
    
    # Apply multiple testing correction AFTER one-tailed conversion
    adj_pvals <- p.adjust(one_tailed_pvals, method = "fdr")
    
    # Store adjusted p-values and log fold change
    pairwise_pvals[, t] <- adj_pvals
    pairwise_logFC[, t] <- raw_pvals[, "logFC"]
    
  }
  
  # Initialize excursion matrix
  excursion_matrix <- matrix(0, nrow = nrow(data), ncol = num_timepoints)
  rownames(excursion_matrix) <- rownames(data)
  colnames(excursion_matrix) <- unique_timepoints
  
  # **Pattern Search**
  for (pattern_idx in seq_along(patterns)) {
    pattern <- patterns[[pattern_idx]]
    pattern_length <- length(pattern)
    
    for (start_idx in 1:(num_timepoints - pattern_length)) {
      
      window_pvals <- pairwise_pvals[, start_idx:(start_idx + pattern_length - 1), drop = FALSE]
      window_logFCs <- pairwise_logFC[, start_idx:(start_idx + pattern_length - 1), drop = FALSE]
      
      is_significant <- rep(TRUE, nrow(data))
      
      for (step in seq_along(pattern)) {
        if (pattern[step] == 0) {
          # If pattern specifies "0", check if logFC is within delta range
          is_significant <- is_significant & 
            (window_pvals[, step] > alpha) & 
            (window_logFCs[, step] >= -delta & window_logFCs[, step] <= delta)
        } else {
          # One-tailed test based on pattern direction
          if (pattern[step] < 0) {
            is_significant <- is_significant & 
              (window_pvals[, step] < alpha) & (window_logFCs[, step] < 0)
          } else {
            is_significant <- is_significant & 
              (window_pvals[, step] < alpha) & (window_logFCs[, step] > 0)
          }
        }
      }
      
      # If the pattern is detected, mark the middle timepoint
      excursion_matrix[, start_idx + 1] <- ifelse(is_significant, 1, excursion_matrix[, start_idx + 1])
    }
  }
  
  results_df <- data.frame(feature_nr = rownames(excursion_matrix), excursion_matrix)
  
  return(list(
    results_df = results_df,
    pairwise_pvals = pairwise_pvals,
    pairwise_logFC = pairwise_logFC
  ))
}



# plot_excursions2 <- function(
#     results,
#     data,
#     meta,
#     meta_replicates_column
# ) {
#   results_df <- results$results_df
#   pairwise_pvals <- results$pairwise_pvals
#   
#   unique_timepoints <- sort(unique(meta$Time))
#   num_timepoints <- length(unique_timepoints)
#   
#   unique_replicates <- unique(meta[[meta_replicates_column]])
#   num_replicates <- length(unique_replicates)
#   
#   symbols_available <- c(21, 22, 23, 24, 25, 7, 8, 10, 12)  
#   symbols <- symbols_available[1:num_replicates]  
#   
#   plots <- list()
#   
#   for (protein_index in which(rowSums(results_df[, -1]) > 0)) {
#     
#     feature_name <- results_df$feature_nr[protein_index]
#     protein_data <- data[feature_name, ]
#     
#     plot_data <- data.frame(
#       Time = meta$Time, 
#       Expression = as.numeric(protein_data), 
#       Replicate = as.factor(meta[[meta_replicates_column]])
#     )  
#     
#     excursion_flags <- as.numeric(results_df[protein_index, -1])
#     timepoints_numeric <- sort(unique(meta$Time))
#     
#     # Mark all timepoints participating in a detected pattern
#     plot_data$Excursion <- ifelse(
#       plot_data$Time %in% timepoints_numeric[which(excursion_flags == 1)],
#       1, 0
#     )
#     
#     plot_data$Excursion[is.na(plot_data$Excursion)] <- 0  
#     
#     plot_data$Point_Type <- factor(
#       ifelse(plot_data$Excursion == 1, "Excursion", "Normal"),
#       levels = c("Normal", "Excursion")
#     )
#     
#     sig_df <- data.frame(Time = numeric(0), Label = character(0), y_pos = numeric(0))
#     sig_lines <- data.frame(x_start = numeric(0), x_end = numeric(0), y = numeric(0))
#     
#     # Identify all time windows where an excursion was found
#     excursion_indices <- which(excursion_flags == 1)
#     
#     for (t in excursion_indices) {
#       if (t < num_timepoints) {
#         p_prev <- pairwise_pvals[protein_index, t - 1]
#         p_next <- pairwise_pvals[protein_index, t]
#         
#         get_stars <- function(p) {
#           if (p < 0.0001) return("****")
#           else if (p < 0.001) return("***")
#           else if (p < 0.01) return("**")
#           else if (p < 0.05) return("*")
#           else return("")
#         }
#         
#         sig_prev <- get_stars(p_prev)
#         sig_next <- get_stars(p_next)
#         
#         max_value <- max(plot_data$Expression[
#           plot_data$Time %in% c(unique_timepoints[t - 1], unique_timepoints[t], unique_timepoints[t + 1])
#         ], na.rm = TRUE)
#         
#         base_y <- max_value + 0.5
#         left_y <- base_y  
#         right_y <- base_y + 0.3  
#         star_offset <- 0.2  
#         
#         if (sig_prev != "") {
#           sig_df <- rbind(sig_df, data.frame(
#             Time = mean(c(unique_timepoints[t - 1], unique_timepoints[t])),
#             Label = sig_prev, 
#             y_pos = left_y + star_offset
#           ))
#           
#           sig_lines <- rbind(sig_lines, data.frame(
#             x_start = unique_timepoints[t - 1], 
#             x_end = unique_timepoints[t] - 0.05,  
#             y = left_y)
#           )
#         }
#         
#         if (sig_next != "") {
#           sig_df <- rbind(sig_df, data.frame(
#             Time = mean(c(unique_timepoints[t], unique_timepoints[t + 1])),
#             Label = sig_next, 
#             y_pos = right_y + star_offset
#           ))
#           
#           sig_lines <- rbind(sig_lines, data.frame(
#             x_start = unique_timepoints[t] + 0.05,  
#             x_end = unique_timepoints[t + 1], 
#             y = right_y
#           ))
#         }
#       }
#     }
#     
#     p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Time, y = Expression)) +
#       ggplot2::geom_point(ggplot2::aes(shape = Replicate, color = Point_Type), size = 3, stroke = 1.2, fill = "white") +  
#       ggplot2::scale_shape_manual(values = symbols) +
#       ggplot2::scale_color_manual(values = c("Normal" = "grey40", "Excursion" = "red")) +  
#       ggplot2::geom_segment(data = sig_lines, ggplot2::aes(x = x_start, xend = x_end, y = y, yend = y), size = 0.7) +
#       ggplot2::geom_text(data = sig_df, aes(x = Time, y = y_pos, label = Label), size = 5, hjust = 0.5) +
#       ggplot2::labs(
#         title = paste("Expression of", feature_name),
#         x = "Time",
#         y = "Expression Level",
#         color = "Pattern",
#         shape = "Replicate"
#       ) +
#       ggplot2::theme_minimal()
#     
#     plots[[feature_name]] <- p
#   }
#   
#   return(plots)
# }





#' Plot Excursions in Time-Series Omics Data
#'
#' @description
#' This function generates scatter plots for features that exhibit excursions
#' in time-series omics data. Excursion points are highlighted in red, while 
#' normal points remain grey. Significance stars are added to indicate the 
#' statistical significance of adjacent pairwise comparisons.
#'
#' @param results A list returned from `detect_excursions()`, containing 
#' `results_df` (excursion matrix) and `pairwise_pvals` (one-tailed p-values).
#' @param data A numeric matrix, where rows correspond to features (e.g., genes, 
#' proteins, metabolites) and columns correspond to samples.
#' @param meta A data frame containing metadata for the samples. Must include
#' a column named `"Time"` that specifies the timepoint for each sample.
#' @param meta_replicates_column A character string specifying the column name 
#' in `meta` that indicates biological replicates.
#'
#' @return A named list of ggplot objects, where each element corresponds to a 
#' feature with detected excursions. Each plot displays the expression levels 
#' across timepoints, with replicates distinguished by different shapes.
#'
#' @details
#' - The function first extracts features with at least one excursion.
#' - Each feature's expression is plotted across time using `ggplot2`.
#' - Replicates are displayed with distinct shapes, while excursion points 
#'   are highlighted in red.
#' - Significance stars (`*`, `**`, `***`, `****`) are placed between 
#'   adjacent timepoints if their pairwise comparison is significant.
#' - The significance stars are positioned above horizontal lines, with 
#'   right-side comparisons always placed slightly higher for clarity.
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_shape_manual 
#' @importFrom ggplot2 scale_color_manual geom_segment geom_text labs 
#'                     theme_minimal
#' 
plot_peaks_valleys <- function(
    results,
    data,
    meta,
    meta_replicates_column
) {
  
  results_df <- results$results_df
  pairwise_pvals <- results$pairwise_pvals
  
  unique_timepoints <- sort(unique(meta$Time))
  num_timepoints <- length(unique_timepoints)
  
  # Extract unique replicates from meta data
  unique_replicates <- unique(meta[[meta_replicates_column]])
  num_replicates <- length(unique_replicates)
  
  # Define symbols dynamically based on the number of replicates
  symbols_available <- c(21, 22, 23, 24, 25, 7, 8, 10, 12)
  symbols <- symbols_available[1:num_replicates]
  
  plots <- list()
  
  for (protein_index in which(rowSums(results_df[, -1]) > 0)) {
    
    feature_name <- results_df$feature_nr[protein_index]
    protein_data <- data[feature_name, ]
    
    # Prepare data for plotting
    plot_data <- data.frame(
      Time = meta$Time,
      Expression = as.numeric(protein_data),
      Replicate = as.factor(meta[[meta_replicates_column]])
    )
    
    excursion_flags <- as.numeric(results_df[protein_index, -1])
    timepoints_numeric <- sort(unique(meta$Time))
    plot_data$Excursion <- excursion_flags[match(
      plot_data$Time,
      timepoints_numeric
    )]
    
    # Convert NA excursions to 0 to avoid missing values
    plot_data$Excursion[is.na(plot_data$Excursion)] <- 0
    
    # Explicitly define the point type for coloring
    plot_data$Point_Type <- ifelse(
      plot_data$Excursion == 1,
      "Excursion",
      "Normal"
    )
    
    plot_data$Point_Type <- factor(
      plot_data$Point_Type,
      levels = c("Normal", "Excursion")
    )
    
    # Get significance stars for only the relevant comparisons
    sig_df <- data.frame(
      Time = numeric(0),
      Label = character(0),
      y_pos = numeric(0)
    )
    
    sig_lines <- data.frame(
      x_start = numeric(0),
      x_end = numeric(0),
      y = numeric(0)
    )
    
    for (t in 2:(num_timepoints - 1)) {
      if (excursion_flags[t] == 1) {
        p_prev <- pairwise_pvals[protein_index, t - 1]
        p_next <- pairwise_pvals[protein_index, t]
        
        # Function to get significance stars
        get_stars <- function(p) {
          if (p < 0.0001) return("****")
          else if (p < 0.001) return("***")
          else if (p < 0.01) return("**")
          else if (p < 0.05) return("*")
          else return("")
        }
        
        sig_prev <- get_stars(p_prev)
        sig_next <- get_stars(p_next)
        
        # Determine highest value including the excursion itself
        max_value <- max(
          plot_data$Expression[plot_data$Time %in% c(
            unique_timepoints[t - 1],
            unique_timepoints[t],
            unique_timepoints[t + 1]
          )],
          na.rm = TRUE
        )
        
        # Set base y-levels slightly above the highest value found
        base_y <- max_value + 0.5
        # Same reference for left and right, ensuring consistent placement
        left_y <- base_y
        right_y <- base_y + 0.3  # Right always slightly higher than left
        star_offset <- 0.2  # Small offset to separate stars from lines
        
        # Separate horizontal lines for left and right comparison
        if (sig_prev != "") {
          sig_df <- rbind(
            sig_df,
            data.frame(
              Time = mean(c(
                unique_timepoints[t - 1],
                unique_timepoints[t]
              )),
              Label = sig_prev,
              y_pos = left_y + star_offset
            )
          )
          
          sig_lines <- rbind(
            sig_lines,
            data.frame(
              x_start = unique_timepoints[t - 1],
              x_end = unique_timepoints[t] - 0.05,
              y = left_y)
          )  # Left comparison
        }
        
        if (sig_next != "") {
          sig_df <- rbind(
            sig_df,
            data.frame(
              Time = mean(c(
                unique_timepoints[t],
                unique_timepoints[t + 1]
              )),
              Label = sig_next,
              y_pos = right_y + star_offset
            )
          )
          
          sig_lines <- rbind(
            sig_lines,
            data.frame(
              x_start = unique_timepoints[t] + 0.05,
              x_end = unique_timepoints[t + 1],
              y = right_y
            )
          )  # Right comparison, slightly higher
        }
      }
    }
    
    # Create scatter plot
    p <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = Time, y = Expression)) +
      ggplot2::geom_point(
        ggplot2::aes(
          shape = Replicate,
          color = Point_Type
        ),
        size = 3,
        stroke = 1.2,
        fill = "white"
      ) +  # White fill ensures visibility
      ggplot2::scale_shape_manual(values = symbols) +
      ggplot2::scale_color_manual(values = c(
        "Normal" = "grey40",
        "Excursion" = "red"
      )) +  # Now explicitly maps colors
      ggplot2::geom_segment(
        data = sig_lines,
        ggplot2::aes(
          x = x_start,
          xend = x_end,
          y = y,
          yend = y
        ),
        size = 0.7
      ) +
      ggplot2::geom_text(
        data = sig_df,
        aes(
          x = Time,
          y = y_pos,
          label = Label
        ),
        size = 5,
        hjust = 0.5
      ) +
      ggplot2::labs(
        title = paste(
          "Feature: ",
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

