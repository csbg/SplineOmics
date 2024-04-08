# Import libraries ---------------------------------------------
library(limma)
library(splines)
library(ggplot2)
library(tidyr)
library(dplyr)
library(tibble)
library(dendextend)
library(RColorBrewer)
library(patchwork)



# Internal functions level 4 ----------------------------------------------


convertPlotToBase64ImgTag <- function(plot, fn5_plot_nrows, width = 7, 
                                      base_height_per_row = 2.5, units = "in", 
                                      html_img_width="100%") {
  
  additional_height_per_row <- 2.1
  height <- base_height_per_row + (fn5_plot_nrows - 1) * 
    additional_height_per_row
  
  # Temporarily save the plot as an image
  img_file <- tempfile(fileext = ".png")
  ggsave(img_file, plot=plot, width = width, height = height, units = units, 
         limitsize = FALSE)
  
  # Convert the image to a Base64 string
  img_base64 <- base64enc::dataURI(file = img_file, mime = "image/png")
  
  # Delete the temporary image file
  unlink(img_file)
  
  # Return the HTML img tag with the Base64 string and a fixed width
  return(sprintf('<img src="%s" alt="Plot" style="width:%s;">', 
                 img_base64, html_img_width))
}

# Internal functions level 3 ---------------------------


buildPlotReportHtml <- function(header_section, plots, rowcounts, path) {
  

  
  index <- 1
  for (plot in plots) {
    plot_nrows <- rowcounts[[index]]
    index <- index + 1
    img_tag <- convertPlotToBase64ImgTag(plot, plot_nrows)
    header_section <- paste(header_section, img_tag, sep="\n")
  }
  
  
  # Close the HTML document
  html_content <- paste(header_section, "</body></html>", sep="\n")
  
  
  dir_path <- dirname(path)
  
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  writeLines(html_content, path)
}



# Internal functions level 2 ---------------------------


plot_dendrogram <- function(hc, k, title) {
  dend <- as.dendrogram(hc)
  
  clusters <- cutree(hc, k)
  
  palette_name <- "Set3" # This can be changed to another palette if desired
  max_colors_in_palette <- brewer.pal.info[palette_name, "maxcolors"]
  colors <- brewer.pal(min(max_colors_in_palette, k), palette_name)
  if (k > max_colors_in_palette) {
    colors <- rep(colors, length.out = k)
  }
  
  dend_colored <- color_branches(dend, k = k, labels_colors = colors)
  dend_colored <- set(dend_colored, "labels", value = NULL)
  
  ggdend <- as.ggdend(dend_colored)
  p_dend <- ggplot(ggdend) + 
    labs(title = title, 
         x = "", y = "") +
    theme_minimal() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank())
}


plot_all_shapes <- function(curve_values, title) {
  time <- as.numeric(colnames(curve_values)[-length(colnames(curve_values))])
  
  clusters <- unique(curve_values$cluster)
  average_curves <- data.frame()
  
  # Loop through each unique cluster value to calculate the average curve
  for (current_cluster in clusters) {
    # Filter rows for the current cluster
    subset_hits <- curve_values[curve_values$cluster == current_cluster, ]
    
    
    last_timepoint <- (which(names(curve_values) == "cluster")) - 1
    
    average_curve <- colMeans(subset_hits[,1:last_timepoint])
    
    # Create a data frame for the average curve with an additional 'Cluster' column
    curve_df <- data.frame(Time = time, Value = average_curve, 
                           cluster = as.factor(current_cluster))
    
    # Bind the curve data frame to the cumulative data frame
    average_curves <- rbind(average_curves, curve_df)
  }
  
  average_curves$cluster <- 
    factor(average_curves$cluster, 
           levels = sort(unique(as.numeric(average_curves$cluster))))
  
  p_curves <- ggplot(average_curves, aes(x = Time, y = Value, color = 
                                           factor(cluster))) +
    geom_line() + 
    ggtitle(title) +
    xlab("Timepoints") + ylab("Values") +
    scale_color_brewer(palette = "Dark2", name = "Cluster") + 
    theme_minimal()
}


plot_single_and_consensus_splines <- function(time_series_data, title) {
  # Transform the dataframe to a long format for ggplot2
  df_long <- as.data.frame(t(time_series_data)) %>%
    rownames_to_column(var = "time") %>%
    pivot_longer(cols = -time, names_to = "feature", 
                 values_to = "intensity") %>%
    arrange(feature) %>%
    mutate(time = as.numeric(time))
  
  # Compute consensus (mean of each column)
  consensus <- colMeans(time_series_data, na.rm = TRUE)
  consensus_df <- data.frame(time = as.numeric(colnames(time_series_data)), 
                             consensus = consensus)
  
  p <- ggplot() +
    geom_line(data = df_long, aes(x = time, y = intensity, group = feature,
                                  colour = "Single Shapes"),
              alpha = 0.3, linewidth = 0.5) +
    geom_line(data = consensus_df, aes(x = time, y = consensus,
                                       colour = "Consensus Shape"),
              linewidth = 1.5) +
    scale_colour_manual("", values = c("Consensus Shape" = "darkblue",
                                       "Single Shapes" = "#6495ED")) +
    theme_minimal() +
    ggtitle(title) +
    xlab("time to feeding [min]") +
    ylab("Y")

  return(p)
}


plot_consensus_shapes <- function(curve_values, title) {
  clusters <- sort(unique(curve_values$cluster))
  time <- as.numeric(colnames(curve_values)[-length(colnames(curve_values))])
  
  plots <- list()
  for (current_cluster in clusters) {
    current_title <- paste(title, current_cluster, sep = "_")
    subset_df <- subset(curve_values, cluster == current_cluster)
    subset_df$cluster <- NULL 
    
    plots[[length(plots) + 1]] <- plot_single_and_consensus_splines(
      subset_df, current_title)
  }
  return(plots)
}


plot_splines <- function(top_table, data, meta, main_title) {
  
  DoF <- which(names(top_table) == "AveExpr") - 1
  time_points <- meta$Time
  
  titles <- data.frame(
    FeatureID = top_table$feature_index,
    feature_names = top_table$feature_names
  )
  
  plot_list <- list()
  
  ## Generate individual plots ----
  for (hit in 1:nrow(top_table)) {
    hit_index <- as.numeric(top_table$feature_index[hit])
    y_values <- data[hit_index, ]
    
    intercept <- top_table$intercept[hit]
    
    spline_coeffs <- as.numeric(top_table[hit, 1:DoF])
    
    Time <- seq(meta$Time[1], meta$Time[length(meta$Time)], length.out = 100)
    X <- ns(Time, df = DoF, intercept = FALSE)
    
    fitted_values <- X %*% spline_coeffs + intercept
    
    plot_data <- data.frame(Time = time_points, 
                            Y = y_values)
    
    plot_spline <- data.frame(Time = Time,
                              Fitted = fitted_values)
    
    x_max <- as.numeric(max(time_points))
    x_extension <- x_max * 0.05 
    
    p <- ggplot() +
      geom_point(data = plot_data, aes(x = Time, y = Y), color = 'blue') +
      geom_line(data = plot_spline, aes(x = Time, y = Fitted), 
                color = 'red') +
      theme_minimal() +
      scale_x_continuous(limits = c(min(time_points), x_max + x_extension)) +
      labs(x = "Time [min]", y = "Intensity")
    
    matched_row <- subset(titles, FeatureID == hit_index)
    title <- as.character(matched_row$feature_name)
    if (is.na(title)) {
      title <- paste("feature:", hit_index)
    }
    
    p <- p + labs(title = title, 
                  x = "Time [min]", y = "Intensity") +
      theme(plot.title = element_text(size = 4),
            axis.title.x = element_text(size = 8), 
            axis.title.y = element_text(size = 8)) +
      annotate("text", x = x_max + (x_extension / 2), y = 
                 max(fitted_values, na.rm = TRUE),
               label = "",
               hjust = 0.5, vjust = 1, size = 3.5, angle = 0, color = "black")
    
    plot_list[[length(plot_list) + 1]] <- p
  }
  
  ## Generate the combined plot ----
  if(length(plot_list) > 0) {           
    num_plots <- length(plot_list)
    ncol <- 3
    nrows <- ceiling(num_plots / ncol)
    
    composite_plot <- patchwork::wrap_plots(plot_list, ncol = 3) + 
      plot_annotation(title = paste(main_title, "| DoF:", DoF),
                      theme = theme(plot.title = element_text(hjust = 0.5, 
                                                              size = 14)))
    return(list(composite_plot = composite_plot, nrows = nrows))
  } else {
    stop("plot_list in function plot_splines splinetime package has length 0!")
  }
}


generate_report_html <- function(plot_list, plot_list_nrows, report_dir) {
  timestamp <- format(Sys.time(), "%d_%m_%Y-%H_%M_%S")
  
  omics_data_type <- "PTX"
  header_text <- paste("limma clustered hits", 
                       omics_data_type, timestamp, sep=" | ")
  
  if (omics_data_type == "PTX") {
    design_text <- "splines | X,data = exp AND stat | limma design: 1 + 
      Phase*X + Reactor | timepoints: E12_TP05_Exponential & 
      E10_TP10_Stationary removed<br>(Note: batch-corrected data used for 
      plotting
      individual blue datapoints (plots with the red spline, below))"
  } else if (omics_data_type == "PPTX") {
    design_text <- "splines | X,data = exp AND stat | limma design: 1 + 
      Phase*X + Reactor | timepoints: all<br>(Note: batch-corrected data used 
      for plotting
      individual blue datapoints (plots with the red spline, below))"
  }
  
  html_content <- paste(
    "<html><head><title>My Plots</title></head><body>",
    "<h1 style='color:red;'>", header_text, "</h1>",
    "<h2>Design</h2>",
    "<p>", design_text, "</p>",
    "</body></html>",
    sep=""
  )
  
  file_name <- sprintf("report_clustered_splines_%s_%s.html",
                       omics_data_type, timestamp)
  
  output_file_path <- here::here(report_dir, file_name)
  
  buildPlotReportHtml(html_content, plot_list, plot_list_nrows, 
                      output_file_path)
}



# Internal functions level 1 ---------------------------


make_clustering_report <- function(all_groups_clustering, group_factors, data, 
                                   meta, p_values, report_dir) {
  
  plot_list <- list()
  plot_list_nrows <- list()
  
  for (group_factor in group_factors) {
    i <- 0
    for (group_clustering in all_groups_clustering) {
      i <- i + 1
      
      curve_values <- group_clustering$curve_values
      
      title <- 
        "Exponential phase\n\nHierarchical Clustering Dendrogram 
         (colors = clusters)"
      dendrogram <- plot_dendrogram(group_clustering$hc, clusters[i], title)
      
      title <- 
        "Average Curves by Cluster (colors not matching with plot above!)"
      p_curves <- plot_all_shapes(curve_values, title)
      
      title <- paste("min-max normalized shapes | ", "cluster", sep = " ")
      consensus_shape_plots <- plot_consensus_shapes(curve_values, title)
      
      main_title <- paste(",cluster:", 1, sep = " ")
      titles <- annotation$First.Protein.Description
      
      top_table <- group_clustering$top_table
      p_value <- p_values[i]
      levels <- as.character(unique(meta[[group_factor]]))
      meta_level <- meta %>% filter(.data[[group_factor]] == levels[i])
      sample_names <- as.character(meta_level$Sample)
      data_level <- data[, colnames(data) %in% sample_names]
      
      composite_plots <- list()
      nrows <- list()
      
      for (nr_cluster in unique(na.omit(top_table$cluster))) {
        main_title <- paste(",cluster:", nr_cluster, sep = " ")
        
        top_table_filt <- top_table %>%
          dplyr::filter(adj.P.Val < p_value, .data$cluster == nr_cluster)
        
        plot_splines_result <- plot_splines(top_table_filt, data_level, 
                                            meta_level, main_title)
        
        composite_plots[[length(composite_plots) + 1]] <- 
          plot_splines_result$composite_plot 
        
        nrows[[length(nrows) + 1]] <- plot_splines_result$nrows 
      }
      
      plot_list <- c(plot_list, list(dendrogram, p_curves), 
                     consensus_shape_plots, composite_plots)
      
      plot_list_nrows <- c(plot_list_nrows, 2, 2, 
                           rep(1, length(consensus_shape_plots)), unlist(nrows))
      
    }
  }
  
  generate_report_html(plot_list, plot_list_nrows, report_dir)
  
}



# Section 5: exported package functions -----------------------------------


cluster_hits <- function(top_tables, data, meta, group_factors, p_values, 
                         clusters, report_dir) {
  
  if (!is.list(top_tables) || !all(sapply(top_tables, is.data.frame))) {
    stop("top_tables must be a list of dataframes")
  } else if ((is.data.frame(data)) && (!is.matrix(data))) {
    data <- as.matrix(data)
  } else if (!is.data.frame(meta) || !"Time" %in% names(meta)) {
    stop("meta must be a dataframe")
  } else if (!(is.character(group_factors))) {
    stop("group_factors must be a character vector")
  } else if (!is.numeric(p_values)) {
    stop("p_values must be a float vector")
  } else if (!is.integer(clusters)) {
    stop("clusters must be a integer vector")
  } else if (!is.character(report_dir)) {
    stop("clusters must be a integer vector")
  }
  
  if (!dir.exists(report_dir)) {
    dir.create(report_dir)
  }
  
  groups <- unique(meta[group_factors])
  all_groups_clustering <- list()
  
  for (i in seq_along(top_tables)) {
    top_table <- top_tables[[i]]
    
    p_value <- p_values[i]
    if (is.na(p_value)) {
      p_value <- p_values[1]
    }
    
    k <- clusters[i]
    
    spline_results_hits <- subset(top_table, adj.P.Val < p_value)
    
    DoF <- which(names(top_table) == "AveExpr") - 1
    
    # Get smooth curves 
    selected_group <- groups[i, ]
    subset_meta <- meta[apply(meta[, group_factors], 1, function(row) {
      all(row == unlist(selected_group))
    }), ]
    
    
    smooth_timepoints <- seq(subset_meta$Time[1], 
                             subset_meta$Time[length(subset_meta$Time)], 
                             length.out = 100)
    
    X <- ns(smooth_timepoints, df = DoF, intercept = FALSE)
    
    columns_to_select <- 1:DoF
    
    splineCoeffs <- spline_results_hits %>% 
      select(all_of(1:DoF)) %>%
      as.matrix()
    
    curve_values <- matrix(nrow = nrow(splineCoeffs), 
                           ncol = length(smooth_timepoints))
    
    for(i in 1:nrow(splineCoeffs)) {
      current_coeffs <- matrix(splineCoeffs[i, ], ncol = ncol(splineCoeffs), 
                               byrow = TRUE)
      
      curve_values[i, ] <- current_coeffs %*% t(X)
    }
    
    curve_values <- as.data.frame(curve_values)
    rownames(curve_values) <- rownames(splineCoeffs)
    
    ## Normalize smooth curves between 0 and 1 --------------------------------- 
    # Apply min-max normalization to each row (curve)
    normalized_curves <- apply(curve_values, 1, function(row) {
      (row - min(row)) / (max(row) - min(row))
    })
    
    # Transpose the result to match the original dataframe structure
    normalized_curves <- t(normalized_curves)
    
    # Replace the original dataframe with the normalized values
    curve_values[,] <- normalized_curves
    
    ## Hierarchical clustering -------------------------------------------------
    distance_matrix <- dist(curve_values, method = "euclidean")
    hc <- hclust(distance_matrix, method = "complete")
    cluster_assignments <- cutree(hc, k = k)
    
    clustered_hits <- data.frame(cluster = cluster_assignments)
    clustered_hits$feature <- rownames(clustered_hits)
    clustered_hits <- clustered_hits[, c("feature", "cluster")]
    
    colnames(curve_values) <- smooth_timepoints
    curve_values$cluster <- cluster_assignments
    
    top_table$cluster <- NA
    top_table$cluster[1:nrow(clustered_hits)] <- 
      as.integer(clustered_hits$cluster)
    
    
    group_clustering <- list(clustered_hits = clustered_hits, 
                             hc = hc, 
                             curve_values = curve_values,
                             top_table = top_table)
    
    all_groups_clustering[[length(all_groups_clustering) + 1]] <- 
      group_clustering 
  }
  
  # Generate HTML report with the clustering results
  make_clustering_report(all_groups_clustering, group_factors, data, meta, 
                         p_values, report_dir)
  
  return(all_groups_clustering)
}
