# Import libraries -------------------------------------------------------------
library(limma)
library(splines)



# Internal functions level 2 ----------------------------------------------


modify_limma_top_table <- function(top_table, feature_names) {
  top_table <- as_tibble(top_table, rownames = "feature_index") %>% 
    select(-feature_index, everything(), feature_index) %>%
    mutate(feature_index = as.numeric(feature_index))
  
  sorted_feature_names <- feature_names[top_table$feature_index]
  top_table <- top_table %>%
    mutate(feature_names = sorted_feature_names)
  
  return(top_table)
}


# Internal functions level 1 ---------------------------------------------------


between_level <- function(data, meta, design, DoF, factor, level2,
                          padjust_method, feature_names) {
  meta$X <- ns(meta$Time, df=DoF, intercept = FALSE)
  design_matrix <- model.matrix(as.formula(design), data = meta)

  fit <- lmFit(data, design_matrix)
  fit <- eBayes(fit)

  factor_only_contrast_coeff <- paste0(factor, level2)
  ttlc_factor_only <- topTable(fit, coef=factor_only_contrast_coeff,
                                   adjust.method=padjust_method, number=Inf)
  ttlc_factor_only <- modify_limma_top_table(ttlc_factor_only, feature_names)
  
  factor_time_contrast_coeffs <- paste0(factor, level2, ":X", seq_len(DoF))
  ttlc_factor_time <- topTable(fit, coef=factor_time_contrast_coeffs,
                               adjust.method=padjust_method, number=Inf)
  ttlc_factor_time <- modify_limma_top_table(ttlc_factor_time, feature_names)
  
  return(list(ttlc_factor_only = ttlc_factor_only, 
              ttlc_factor_time = ttlc_factor_time))
}


within_level <- function(data, meta, factor, level, DoF, padjust_method) {
  meta$X <- ns(meta$Time, df=DoF, intercept = FALSE)
  design_matrix <- model.matrix(as.formula(design), data = meta)
  
  fit <- lmFit(data, design_matrix)
  fit <- eBayes(fit)
  
  coef_names <- paste0("X", seq_len(DoF))

  top_table <- topTable(fit, adjust=padjust_method, number=Inf, coef=coef_names)
  
  return(list(top_table = top_table, fit = fit))
}


process_top_table <- function(top_table, top_tables, feature_names, fit, factor, 
                              level) {
  top_table <- modify_limma_top_table(top_table, feature_names)
  
  intercepts <- as.data.frame(coef(fit)[, "(Intercept)", drop = FALSE])
  intercepts_ordered <- intercepts[match(top_table$feature_index, 
                                         rownames(intercepts)), , 
                                   drop = FALSE]
  top_table$intercept <- intercepts_ordered[, 1]
  
  results_name <- paste(factor, level, sep = "_")
  top_tables[[results_name]] <- top_table
  
  return(top_tables)
}



# Export functions -------------------------------------------------------------


run_limma_splines <- function(data, meta, design, DoFs, factors, 
                              feature_names, mode, padjust_method = "BH") {
  # Control arguments
  if ((is.data.frame(data)) && (!is.matrix(data))) {
    data <- as.matrix(data)
  } else if (!is.data.frame(meta) || !"Time" %in% names(meta)) {
    stop("meta must be a dataframe and contain the column 'Time'.")
  } else if (!is.character(design) || length(design) != 1) {
    stop("design must be a single string.")
  } else if (!is.integer(DoFs)) {
    stop("DoFs must be a integer vector")
  } else if (!(is.character(factors))) {
    stop("factors must be a non-empty character vector")
  } else if (!(is.character(feature_ids))) {
    stop("feature_ids must be a non-empty character vector")
  } else if (!(mode == "integrated") || !(mode == "isolated")) {
    stop("mode must be either 'integrated' or 'isolated'")
  }
  
  top_tables <- list()
  ttslc_factor_only <- list()        # ttslc = top_tables_level_comparison
  ttslc_factor_time <- list()        # Factor AND time
  
  for (factor in factors) {             # for example fermentation phase
    meta[[factor]] <- factor(meta[[factor]])
    levels <- levels(meta[[factor]])
    
    # Within level analysis ----------------------------------------------------
    for (i in 1:(length(levels))) { 
      level <- levels[i]  # level is for example exponential or stationary phase
      
      DoF <- DoFs[i]
      if (is.na(DoF)) {
        DoF <- DoFs[1]
      }
      
      if (mode == "isolated") {
        samples <- which(meta[[factor]] == level)
        data_copy <- data[, samples]
        meta_copy <- subset(meta, meta[[factor]] == level)
      } else if (mode == "integrated") {
        data_copy <- data
        meta_copy <- meta
        meta_copy[[factor]] <- relevel(meta_copy[[factor]], ref = level)
        
        level2 <- levels[i + 1]
        if (!is.na(level2)) {
          debug(between_level)
          result <- between_level(data, meta, design, DoF, factor, level2, 
                          padjust_method, feature_names)
          
          ttslc_factor_only[[paste0(level, "_vs_", level2)]] <- 
            result$ttlc_factor_only
          ttslc_factor_time[[paste0(level, "_vs_", level2)]] <- 
            result$ttlc_factor_time
        }
      }
      
      result <- within_level(data_copy, meta_copy, factor, level, DoF, 
                             padjust_method)
      top_table <- result$top_table
      fit <- result$fit
      
      top_tables <- process_top_table(top_table, top_tables, feature_names,
                                      fit, factor, level)
    }
  }
  
  if (mode == "isolated") {
    print("mode == 'integrated' necessary for between level comparison. 
          Returning emtpy list for ttslc_factor_only and ttslc_factor_time
          (ttslc = top_tables_level_comparison).")
  }
  
  return(list(top_tables = top_tables, 
              ttslc_factor_only = ttslc_factor_only,
              ttslc_factor_time = ttslc_factor_time))
}
