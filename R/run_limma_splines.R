# Import libraries --------------------------------------------------------
library(limma)
library(splines)



# Internal functions level 1 ---------------------------


run_limma_splines_integrated <- function(data, meta, design, DoFs, 
                                         group_factors, feature_ids) {
  
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
      
      X <- ns(meta$Time, df=DoF, intercept = FALSE)
      
      # Relevel the factor with the current level as the reference
      meta_releveled <- meta
      meta_releveled[[factor]] <- relevel(meta_releveled[[factor]], 
                                            ref = level)
      
      design_matrix <- model.matrix(as.formula(design), data = meta_releveled)
      
      fit <- lmFit(data, design_matrix)
      fit <- eBayes(fit)
      
      coef_names <- paste0("X", seq_len(DoF))
      top_table <- topTable(fit, adjust="BH", number=Inf, coef=coef_names)
      
      intercepts <- as.data.frame(coef(fit)[, "(Intercept)", drop = FALSE])
      intercepts_ordered <- intercepts[match(rownames(top_table), 
                                             rownames(intercepts)), , 
                                       drop = FALSE]
      top_table$intercept <- intercepts_ordered[, 1]
      
      feature_indices <- as.integer(rownames(top_table))
      sorted_feature_ids <- vector("numeric", length(feature_ids))
      for (i in 1:length(feature_indices)) {
        print(feature_indices[i])
        print(feature_ids[i])
        sorted_feature_ids[i] <- feature_ids[feature_indices[i]]
      }
      top_table$feature_id <- sorted_feature_ids
      
      results_name <- paste(factor, level, sep = "_")
      top_tables[[results_name]] <- top_table
    }
  }
  
  return(top_tables)
}


run_limma_splines_isolated <- function(data, meta, design, DoFs, 
                                       group_factors, feature_ids) {
  
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
      
      meta_level <- subset(meta, meta[[factor]] == level)
      samples <- which(meta[[factor]] == level)
      data_level <- data[, samples]
      
      X <- ns(meta_level$Time, df=DoF, intercept = FALSE)

      design_matrix <- model.matrix(as.formula(design), data = meta_level)
      
      fit <- lmFit(data_level, design_matrix)
      fit <- eBayes(fit)
      
      coef_names <- paste0("X", seq_len(DoF))
      top_table <- topTable(fit, adjust="BH", number=Inf, coef=coef_names)
      
      intercepts <- as.data.frame(coef(fit)[, "(Intercept)", drop = FALSE])
      intercepts_ordered <- intercepts[match(rownames(top_table), 
                                             rownames(intercepts)), , 
                                       drop = FALSE]
      top_table$intercept <- intercepts_ordered[, 1]
      
      feature_indices <- as.integer(rownames(top_table))
      sorted_feature_ids <- vector("numeric", length(feature_ids))
      for (i in 1:length(feature_indices)) {
        print(feature_indices[i])
        print(feature_ids[i])
        sorted_feature_ids[i] <- feature_ids[feature_indices[i]]
      }
      top_table$feature_id <- sorted_feature_ids
      
      results_name <- paste(factor, level, sep = "_")
      top_tables[[results_name]] <- top_table
    }
  }
  
  return(top_tables)
}



# Export functions -----------------------------------


run_limma_splines <- function(data, meta, design, DoFs, group_factors, feature_ids,
                              mode) {
  # Control arguments
  if ((is.data.frame(data)) && (!is.matrix(data))) {
    data <- data %>%
      lapply(function(x) as.numeric(as.character(x))) %>%
      data.frame() %>%
      as.matrix()
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
  }
  
  
  if (mode == "integrated") {
    return(run_limma_splines_integrated(data, meta, design, DoFs, group_factors, 
                                        feature_ids))
  } 
  else if (mode == "isolated") {
    return(run_limma_splines_isolated(data, meta, design, DoFs, group_factors, 
                                      feature_ids))
  } 
  else {
    stop("mode must be either 'integrated' (get hits by considering the data
         from all groups and focusing with reference level) or 'isolated' 
         (get hits only from group data)")
  }
}


