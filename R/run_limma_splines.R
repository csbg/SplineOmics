# Import libraries --------------------------------------------------------
library(limma)
library(splines)



# Export functions -----------------------------------


run_limma_splines <- function(data, meta, design, DoFs, group_factors, 
                              feature_ids, mode) {
  # Control arguments
  if ((is.data.frame(data)) && (!is.matrix(data))) {
    data <- as.matrix(data)
  } else if (!is.data.frame(meta) || !"Time" %in% names(meta)) {
    stop("meta must be a dataframe and contain the column 'Time'.")
  } else if (!is.character(design) || length(design) != 1) {
    stop("design must be a single string.")
  } else if (!is.integer(DoFs)) {
    stop("DoFs must be a integer vector")
  } else if (!(is.character(group_factors))) {
    stop("group_factors must be a non-empty character vector")
  } else if (!(is.character(feature_ids))) {
    stop("feature_ids must be a non-empty character vector")
  } else if (!(mode == "integrated") || !(mode == "isolated")) {
    stop("mode must be either 'integrated' or 'isolated'")
  }
  
  top_tables <- list()
  
  for (factor in group_factors) {
    meta[[factor]] <- factor(meta[[factor]])
    levels <- levels(meta[[factor]])
    
    i <- 0
    for (level in levels) {
      i <- i + 1
      DoF <- DoFs[i]
      if (is.na(DoF)) {
        DoF <- DoFs[1]
      }
      
      data_copy <- data
      meta_copy <- meta
      
      if (mode == "isolated") {
        meta_copy <- subset(meta, meta[[factor]] == level)
        samples <- which(meta[[factor]] == level)
        data_copy <- data[, samples]
      }
      
      X <- ns(meta_copy$Time, df=DoF, intercept = FALSE)
      
      if (mode == "integrated") {
        # Relevel the factor with the current level as the reference
        meta_copy[[factor]] <- relevel(meta_copy[[factor]], ref = level)
      }
      
      design_matrix <- model.matrix(as.formula(design), data = meta_copy)
      
      fit <- lmFit(data_copy, design_matrix)
      fit <- eBayes(fit)
      
      coef_names <- paste0("X", seq_len(DoF))
      top_table <- topTable(fit, adjust="BH", number=Inf, coef=coef_names) %>% 
        as_tibble(rownames = "feature_index") %>% 
        select(-feature_index, everything(), feature_index) %>%
        mutate(feature_index = as.numeric(feature_index))
      
      sorted_feature_names <- feature_names[top_table$feature_index]
      top_table <- top_table %>%
        mutate(feature_names = sorted_feature_names)
      
      intercepts <- as.data.frame(coef(fit)[, "(Intercept)", drop = FALSE])
      intercepts_ordered <- intercepts[match(top_table$feature_index, 
                                             rownames(intercepts)), , 
                                       drop = FALSE]
      top_table$intercept <- intercepts_ordered[, 1]
      
      results_name <- paste(factor, level, sep = "_")
      top_tables[[results_name]] <- top_table
    }
  }
  
  return(top_tables)
}
