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
  fit <- lmFit(data, design)
  
  # Store pairwise p-values
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
    
    # Store adjusted p-values
    pairwise_pvals[, t] <- limma::topTable(
      fit2,
      coef = 1,
      number = Inf,
      adjust.method = "fdr",
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
  
  # Detect excursions
  for (i in 1:nrow(data)) {
    for (t in 2:(num_timepoints - 1)) {
      p_prev <- pairwise_pvals[i, t - 1]  # t-1 vs t
      p_next <- pairwise_pvals[i, t]  # t vs t+1
      
      prev_mean <- mean(data[i, which(meta$Time == unique_timepoints[t - 1])])
      curr_mean <- mean(data[i, which(meta$Time == unique_timepoints[t])])
      next_mean <- mean(data[i, which(meta$Time == unique_timepoints[t + 1])])
      
      prev_change <- curr_mean - prev_mean
      next_change <- next_mean - curr_mean
      
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


plot_excursions <- function(
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
  symbols_available <- c(21, 22, 23, 24, 25, 7, 8, 10, 12)  # Various empty shapes
  symbols <- symbols_available[1:num_replicates]  # Assign only as many as needed
  
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
          plot_data$Expression[plot_data$Time %in% c(unique_timepoints[t - 1], 
          unique_timepoints[t], 
          unique_timepoints[t + 1])], 
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
          "Expression of",
          feature_name
          ),
        x = "Time",
        y = "Expression Level",
        color = "Excursion Status",
        shape = "Replicate"
        ) +
      ggplot2::theme_minimal()
    
    plots[[feature_name]] <- p
  }
  
  return(plots)
}










