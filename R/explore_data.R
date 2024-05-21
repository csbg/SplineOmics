# Exported function: explore_data () -------------------------------------------


#' Generate Exploratory Plots
#'
#' @description
#' This function takes a data matrix, checks its validity, and generates a list 
#' of exploratory plots including density plots, boxplots, PCA plots, MDS plots, 
#' variance explained plots, and violin plots.
#'
#' @param data A numeric matrix containing the data.
#'
#' @return A list of ggplot objects representing various exploratory plots.
#'
#' @examples
#' # Example usage
#' data <- matrix(rnorm(1000), ncol = 10)
#' plots <- generate_plots(data)
#' for (plot in plots) {
#'   print(plot)
#' }
#'
#' @export
#' 
explore_data <- function(data,
                         meta,
                         condition,
                         report_info,
                         meta_batch_column = NA,
                         report_dir = here::here()) {
  
  # Validate inputs
  check_data_and_meta(data,
                      meta,
                      condition,
                      meta_batch_column)
  
  check_report_info(report_info)
  check_and_create_report_dir(report_dir)
  
  
  meta[[condition]] <- as.factor(meta[[condition]])
  
  density_plots <- make_density_plots(data,
                                      meta,
                                      condition)
  
  box_plots <- make_box_plots(data,
                              meta,
                              condition)
  
  violin_plots <- make_violin_plots(data,
                                    meta,
                                    condition) 
  
  pca_plot <- make_pca_plot(data, 
                            meta,
                            condition)
  
  pca_variance_explained_plot <- make_pca_variance_explained_plot(data)
  
  # mds_plot <- make_mds_plot(data)
  
  corr_heatmaps <- make_correlation_heatmaps(data, 
                                             meta, 
                                             condition)
  
  plots <- c(density_plots, 
             box_plots, 
             violin_plots,
             list(pca_plot, pca_variance_explained_plot),
             corr_heatmaps)
  
  plots_sizes <- c(rep(1, length(density_plots)),
                   rep(1.5, length(box_plots)),
                   rep(1.5, length(violin_plots)),
                   1.5,
                   1.5,
                   rep(1.5, length(corr_heatmaps)))
  
  generate_report_html(plots = plots,
                       plots_sizes = plots_sizes,
                       report_info = report_info,
                       filename = "explore_data",
                       report_dir = report_dir)
  
  return(plots)
}


# Level 1 internal functions ---------------------------------------------------


#' Generate Density Plot
#'
#' @description
#' This function generates a density plot for a given data matrix. The density
#' plot shows the distribution of the values in the data matrix.
#'
#' @param data A numeric matrix containing the data.
#'
#' @return A ggplot object representing the density plot.
#'
#' @examples
#' # Example usage
#' data <- matrix(rnorm(1000), ncol = 10)
#' density_plot <- make_density_plot(data)
#' print(density_plot)
#'
#' @importFrom ggplot2 ggplot geom_density ggtitle aes
#' @importFrom reshape2 melt
#'
make_density_plots <- function(data, 
                               meta, 
                               condition) {
  
  density_plots <- list()
  
  # Melt the data to long format
  data_long <- reshape2::melt(as.data.frame(data), id.vars = NULL)
  
  if (length(unique(meta[[condition]])) > 1) {    # Only when > 2 levels
    # Create the overall density plot for all data
    overall_plot <- ggplot2::ggplot(data_long, 
                                    ggplot2::aes(x = value)) +
      ggplot2::geom_density(fill = "blue", alpha = 0.5) +
      ggplot2::ggtitle("All Levels")
    
    density_plots <- c(density_plots, list(overall_plot))
  }

  
  # Create density plots for each level of the condition
  levels <- unique(meta[[condition]])
  for (level in levels) {
    # Filter the data for the current level
    indices <- which(meta[[condition]] == level)
    data_level <- data[, indices, drop = FALSE]
    data_level_long <- reshape2::melt(as.data.frame(data_level), id.vars = NULL)
    
    # Create the density plot for the current level
    level_plot <- ggplot2::ggplot(data_level_long, 
                                  ggplot2::aes(x = value)) +
      ggplot2::geom_density(fill = "blue", alpha = 0.5) +
      ggplot2::ggtitle(paste("Level:", level))
    
    # Add the level plot to the list
    density_plots <- c(density_plots, list(level_plot))
  }
  
  return(density_plots)
}


#' Generate Boxplot
#'
#' @description
#' This function generates a boxplot for a given data matrix. The boxplot shows 
#' the distribution of the values in the data matrix across different variables, 
#' with log2 intensity on the y-axis.
#'
#' @param data A numeric matrix containing the data.
#'
#' @return A ggplot object representing the boxplot.
#'
#' @examples
#' # Example usage
#' data <- matrix(rnorm(1000), ncol = 10)
#' boxplot <- make_boxplot(data)
#' print(boxplot)
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot labs theme element_text
#' @importFrom reshape2 melt
#'
make_box_plots <- function(data,
                           meta,
                           condition) {

  boxplots <- list()
  
  # Create boxplots for each level of the condition
  levels <- unique(meta[[condition]])
  for (level in levels) {
    # Filter the data for the current level
    indices <- which(meta[[condition]] == level)
    data_level <- data[, indices, drop = FALSE]
    data_level_long <- reshape2::melt(as.data.frame(data_level), id.vars = NULL)
    
    # Create the boxplot for the current level
    level_plot <- ggplot2::ggplot(data_level_long, 
                                  ggplot2::aes(x = variable, y = value)) +
      ggplot2::geom_boxplot() +
      ggplot2::labs(x = "Variables", y = "log2 intensity", 
                    title = paste("Level:", level)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1,
                                                         size = 6))
    
    # Add the level plot to the list
    boxplots <- c(boxplots, list(level_plot))
  }
  
  return(boxplots)
}


#' Generate Violin Plot
#'
#' @description
#' This function generates a violin plot for a given data matrix. The violin 
#' plot shows the distribution of the values in the data matrix across different 
#' variables, with each variable's distribution displayed as a separate violin.
#'
#' @param data A numeric matrix containing the data.
#'
#' @return A ggplot object representing the violin plot.
#'
#' @examples
#' # Example usage
#' data <- matrix(rnorm(1000), ncol = 10)
#' violin_plot <- make_violin_plot(data)
#' print(violin_plot)
#' 
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes geom_violin theme labs
#' @importFrom grid unit
#' 
make_violin_plots <- function(data,
                              meta,
                              condition) {
  
  violin_plots <- list()
  
  # Create violin plots for each level of the condition
  levels <- unique(meta[[condition]])
  for (level in levels) {
    # Filter the data for the current level
    indices <- which(meta[[condition]] == level)
    data_level <- data[, indices, drop = FALSE]
    data_level_long <- reshape2::melt(as.data.frame(data_level), id.vars = NULL)
    
    # Create the violin plot for the current level
    level_plot <- ggplot2::ggplot(data_level_long, 
                                  ggplot2::aes(x = variable, y = value)) +
      ggplot2::geom_violin(trim = FALSE) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1,
                                                         size = 6),
                     plot.margin = grid::unit(c(0, 0, 0, 0), "cm")) +
      ggplot2::labs(x = "Timepoint", y = "Value", 
                    title = paste("Level:", level))
    
    # Add the level plot to the list
    violin_plots <- c(violin_plots, list(level_plot))
  }
  
  return(violin_plots)
}


#' Generate PCA Plot with Dynamic Coloring
#'
#' @description
#' This function generates a PCA plot from a data matrix, dynamically coloring
#' the points based on the levels of a specified factor in the metadata.
#'
#' @param data A numeric matrix containing the data.
#' @param meta A dataframe containing the metadata.
#' @param condition The column name in the metadata dataframe that contains the
#' factor levels for coloring the PCA plot.
#'
#' @return A ggplot object representing the PCA plot.
#'
#' @examples
#' # Example usage
#' data <- matrix(rnorm(1000), ncol = 10)
#' meta <- data.frame(Sample = colnames(data), Condition = rep(c("A", "B"), 
#' each = 5))
#' plots <- make_pca_plot(data, meta, "Condition")
#' print(plots)
#' 
#' @importFrom ggplot2 ggplot geom_point xlim xlab ylab ggtitle theme_minimal
#' theme
#' @importFrom ggrepel geom_text_repel
#' 
make_pca_plot <- function(data, 
                          meta, 
                          condition) {
  
  pc <- prcomp(t(data))
  pca_df <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2])
  
  # Get the labels and colors from the metadata
  pca_df$Labels <- colnames(data)
  # meta <- meta[match(pca_df$Labels, meta$Sample), ]
  pca_df$LabelColor <- meta[[condition]]
  
  # Calculate the variance explained
  variance_explained <- pc$sdev^2 / sum(pc$sdev^2)
  percent_variance_explained <- round(variance_explained * 100, digits = 1)
  
  # Calculate the x-axis range and extend it
  x_range <- range(pc$x[, 1])
  extended_x_max <- x_range[2] + (x_range[2] - x_range[1]) * 0.2
  
  # Create the PCA plot
  pca_plot <- ggplot2::ggplot(pca_df, aes(x = PC1, y = PC2, 
                                          color = LabelColor)) +
    ggplot2::geom_point() +
    ggrepel::geom_text_repel(aes(label = Labels), 
                    box.padding = 0.35, 
                    point.padding = 0.5, 
                    max.overlaps = Inf,
                    size = 2) +
    ggplot2::xlim(x_range[1], extended_x_max) +
    ggplot2::xlab(paste("PC1 -", percent_variance_explained[1], "% variance")) +
    ggplot2::ylab(paste("PC2 -", percent_variance_explained[2], "% variance")) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = element_text(hjust = 0.5))
}


#' Generate PCA Variance Explained Plot
#'
#' @description
#' This function generates a bar plot showing the variance explained by the 
#' first seven principal components of a given data matrix.
#'
#' @param data A numeric matrix containing the data.
#'
#' @return A ggplot object representing the variance explained by the principal 
#' components.
#'
#' @examples
#' # Example usage
#' data <- matrix(rnorm(1000), ncol = 10)
#' variance_plot <- make_pca_variance_explained_plot(data)
#' print(variance_plot)
#' 
#' @importFrom limma plotMDS
#' @importFrom ggplot2 ggplot aes geom_col geom_text xlab ylab ggtitle 
#' theme_minimal
#' 
make_pca_variance_explained_plot <- function(data) {
  # Perform PCA
  pc <- prcomp(t(data))
  
  # Calculate the variance explained
  variance_explained <- pc$sdev^2 / sum(pc$sdev^2)
  percent_variance_explained <- round(variance_explained * 100, digits = 1)
  
  # Create a dataframe for the first 7 principal components
  var_explained <- data.frame(
    PC = paste0("PC", 1:7),
    Variance = percent_variance_explained[1:7]
  )
  
  # Create the variance explained plot
  variance_plot <- ggplot2::ggplot(var_explained, 
                                   ggplot2::aes(x = PC, y = Variance)) +
    ggplot2::geom_col() +
    ggplot2::geom_text(aes(label = round(Variance, digits = 2)), 
                       vjust = -0.2) +
    ggplot2::xlab("Principal Component") +
    ggplot2::ylab("Variance Explained (%)") +
    ggplot2::theme_minimal() +
    ggplot2::ylim(0, max(var_explained$Variance) + 5)  # Space above highest bar
  
  return(variance_plot)
}


#' Generate MDS Plot
#'
#' @description
#' This function generates a multidimensional scaling (MDS) plot for a given 
#' data matrix. The MDS plot visualizes the similarities or dissimilarities 
#' between samples in the data matrix.
#'
#' @param data A numeric matrix containing the data.
#'
#' @return A ggplot object representing the MDS plot.
#'
#' @examples
#' # Example usage
#' data <- matrix(rnorm(1000), ncol = 10)
#' mds_plot <- make_mds_plot(data)
#' print(mds_plot)
#' 
#' @importFrom limma plotMDS
#' @importFrom ggplot2 ggplot geom_point ggtitle theme_minimal
#' @importFrom ggrepel geom_text_repel
#' 
make_mds_plot <- function(data) {
  
  mds <- limma::plotMDS(x = data, plot = FALSE)
  mds_df <- data.frame(Dim1 = mds$x[, 1], Dim2 = mds$x[, 2], 
                       Labels = colnames(data))
  
  mds_plot <- ggplot2::ggplot(mds_df, 
                              ggplot2::aes(x = Dim1, 
                                           y = Dim2, 
                                           label = Labels)) +
    ggplot2::geom_point() +
    ggrepel::geom_text_repel() +
    ggplot2::theme_minimal()
}


#' Generate Correlation Heatmaps
#'
#' @description
#' This function generates correlation heatmaps using Spearman correlation for 
#' a given data matrix. It creates a combined heatmap for all levels and 
#' individual heatmaps for each level specified in the condition column of the 
#' metadata.
#'
#' @param data A numeric matrix containing the data.
#' @param meta A dataframe containing the metadata.
#' @param condition The column name in the metadata dataframe that contains the 
#' factor levels for generating individual heatmaps.
#'
#' @return A list of `ComplexHeatmap` heatmap objects representing the 
#' correlation 
#' heatmaps.
#'
#' @examples
#' # Example usage
#' data <- matrix(rnorm(1000), ncol = 10)
#' colnames(data) <- paste0("Sample", 1:10)
#' meta <- data.frame(Sample = colnames(data), Condition = rep(c("A", "B"), 
#' each = 5))
#' heatmaps <- make_correlation_heatmaps(data, meta, "Condition")
#' for (heatmap in heatmaps) {
#'   draw(heatmap, heatmap_legend_side = "right")
#' }
#'
#' @importFrom circlize colorRamp2
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom RcolorBrewer
#' 
make_correlation_heatmaps <- function(data, 
                                      meta, 
                                      condition) {

  heatmaps <- list()
  
  if (length(unique(meta[[condition]])) > 1) {    # Only when > 2 levels
    # Create combined correlation heatmap
    corr_all <- cor(data, method = "spearman", use = "pairwise.complete.obs")
    # Remove perfect correlations for better visualization
    diag(corr_all) <- NA
  
    # Define the color function for the heatmap based on the range of corr_all
    breaks <- seq(min(corr_all, na.rm = TRUE), 
                  max(corr_all, na.rm = TRUE), length.out = 100)
    col_fun <- 
      circlize::colorRamp2(breaks, 
                           colorRampPalette(c("blue", "white", "red")
                                            )(length(breaks)))
    
    heatmap_all <- ComplexHeatmap::Heatmap(
      corr_all,
      col = col_fun,
      name = "Correlation",
      column_title = "All Levels",
      heatmap_legend_param = list(
        title = "Spearman Correlation",
        title_gp = gpar(fontsize = 6),  
        labels_gp = gpar(fontsize = 6),  
        title_position = "leftcenter-rot" 
      ),
      na_col = "grey",
      row_names_gp = gpar(fontsize = 6),
      column_names_gp = gpar(fontsize = 6)
    )
    heatmaps <- c(heatmaps, list(heatmap_all))
  }
  
  # Create correlation heatmaps for each level of the condition
  levels <- unique(meta[[condition]])
  for (level in levels) {
    # Filter the data for the current level
    indices <- which(meta[[condition]] == level)
    data_level <- data[, indices, drop = FALSE]
    
    # Compute Spearman correlation
    corr_level <- cor(data_level, method = "spearman", 
                      use = "pairwise.complete.obs")
    # Remove perfect correlations for better visualization
    diag(corr_level) <- NA

    # Define the color function for the heatmap based on the range of corr_level
    breaks_level <- seq(min(corr_level, na.rm = TRUE), 
                        max(corr_level, na.rm = TRUE), length.out = 100)
    
    col_fun_level <- 
      circlize::colorRamp2(breaks_level, 
                          colorRampPalette(c("blue", "white", "red")
                                           )(length(breaks_level)))
    
    # Create the correlation heatmap for the current level
    heatmap_level <- ComplexHeatmap::Heatmap(
      corr_level,
      col = col_fun_level,
      name = "Correlation",
      column_title = paste("Level:", level),
      heatmap_legend_param = list(
        title = "Spearman Correlation",
        title_gp = gpar(fontsize = 6),  
        labels_gp = gpar(fontsize = 6), 
        title_position = "leftcenter-rot"  
      ),
      na_col = "grey",
      row_names_gp = gpar(fontsize = 6),
      column_names_gp = gpar(fontsize = 6),
    )

    heatmaps <- c(heatmaps, list(heatmap_level))
  }
  
  return(heatmaps)
}




#' Build Explore Data Report
#'
#' @description
#' This function generates an HTML report containing a header section, table of 
#' contents, and a series of plots. Each plot is included in the report with 
#' specified sizes.
#'
#' @param header_section A string containing the HTML content for the header 
#' section of the report.
#' @param plots A list of ggplot objects representing the plots to be included 
#' in the report.
#' @param plots_sizes A list of sizes corresponding to each plot, defining the 
#' dimensions to be used when rendering the plots.
#' @param output_file_path A string specifying the file path where the HTML 
#' report will be saved.
#'
#' @return None. This function writes the HTML content to the specified file.
#'
#' @examples
#' # Example usage
#' header_section <- "<h1>Exploratory Data Analysis Report</h1>"
#' plots <- list(make_density_plot(data), make_boxplot(data))
#' plots_sizes <- list(c(800, 600), c(800, 600))
#' output_file_path <- "report.html"
#' build_explore_data_report(header_section, plots, plots_sizes, 
#' output_file_path)
#' 
build_explore_data_report <- function(header_section, 
                                      plots, 
                                      plots_sizes, 
                                      output_file_path) {  
  
  level_count = (length(plots)-4)/4
  
  content_with_plots <- paste(header_section, "<!--TOC-->", sep="\n")
  
  toc <- "<div id='toc' style='text-align: center; display: block; margin: auto;
          width: 80%;'> 
        <h2 style='font-size: 40px;'>Table of Contents</h2>
        <ul style='display: inline-block; text-align: left;'>"
  
  section_header_style <- "font-size: 70px; color: #001F3F; text-align: center;"
  toc_style <- "font-size: 30px;"

  pb <- create_progress_bar(plots)
  plot_names <- c("Density Plots", 
                  "Boxplots", 
                  "Violin Plots",
                  "PCA Plot", 
                  "PCA Plot Variance explained",
                  "Correlation Heatmaps")

  toc_index <- 0
  toc_index_memory <- toc_index
  
  # Generate the sections and plots
  for (index in seq_along(plots)) {
    
    # Determine when to place headers based on the provided logic
    if (length(plots) == 6) {     # Just one level exists
      
      toc_index <- toc_index + 1
      
    } else if (index == 1 || 
               index == 2 + level_count || 
               index == 2 + 2 * level_count || 
               index == 2 + 3 * level_count || 
               index == 3 + 3 * level_count || 
               index == 4 + 3 * level_count) {    # More than just one level
    
      toc_index <- toc_index + 1
    }
    
    
    if (toc_index != toc_index_memory) {
      
      section_id <- paste0("section_", toc_index)
      toc <- paste(toc, 
                   sprintf('<li style="%s"><a href="#%s">%s</a></li>', 
                           toc_style, 
                           section_id,
                           plot_names[toc_index]), 
                   sep = "\n")
      
      section_header <- sprintf('<h2 id="%s" style="%s">%s</h2>', 
                                section_id, 
                                section_header_style, 
                                plot_names[toc_index])
      
      content_with_plots <- paste(content_with_plots, section_header, 
                                  sep = "\n")
      
      toc_index_memory <- toc_index
    }
    
    
    # Process each plot
    plot <- plots[[index]]
    plot_size <- plots_sizes[[index]]
    img_tag <- plot2base64(plot, plot_size)
    
    content_with_plots <- paste(content_with_plots, img_tag, sep = "\n")
    pb$tick()
  }
  
  # Close the Table of Contents
  toc <- paste(toc, "</ul></div>", sep="\n")
  
  # Insert the Table of Contents at the placeholder
  content_with_plots <- gsub("<!--TOC-->", toc, content_with_plots)
  
  # Append the final closing tags for the HTML body and document
  content_with_plots <- paste(content_with_plots, "</body></html>", sep="\n")
  
  # Ensure the directory exists
  dir_path <- dirname(output_file_path)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  # Write the HTML content to file
  writeLines(content_with_plots, output_file_path)
}