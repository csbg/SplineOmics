#' Correlate and compare hits between two SplineOmics results
#'
#' @description
#' Compares adjusted p-values from two SplineOmics result objects across
#' standard categories. Calculates Spearman correlation globally and for
#' significant hits only, generates diagnostic plots, and summarizes
#' overlap statistics including Jaccard index and hit counts.
#' 
#' Note: The feature_names columns in the topTables of splineomics1 and 
#' splineomics2 must have the same names, otherwise, no comparisons are 
#' possible.
#'
#' @param splineomics1 A list containing `limma_splines_result` with topTables
#'   organized by category and subcategory.
#' @param splineomics2 Same structure as `splineomics1`; the object to compare.
#' @param splineomics1_description Character label used for x-axis annotation.
#' @param splineomics2_description Character label used for y-axis annotation.
#' @param adj_p_tresh1 Numeric p-value threshold for calling hits in analysis 1.
#' @param adj_p_tresh2 Numeric p-value threshold for calling hits in analysis 2.
#'
#' @return A list with components:
#' \describe{
#'   \item{correlation_summary}{Data frame of global and hit-level correlations}
#'   \item{plots}{Named list of ggplot objects for all and hit-only subsets.}
#'   \item{hits_summary}{Data frame with hit counts, overlap, and Jaccard index}
#' }
#' 
#' @export
#' 
compare_results <- function(
    splineomics1,
    splineomics2,
    splineomics1_description = "SplineOmics_1",
    splineomics2_description = "SplineOmics_2",
    adj_p_tresh1 = 0.05,
    adj_p_tresh2 = 0.05
) {
  
  validate_compare_inputs(
    splineomics1,
    splineomics2,
    splineomics1_description,
    splineomics2_description,
    adj_p_tresh1,
    adj_p_tresh2
  )
  
  categories <- c(    # limma result categories
    "time_effect",
    "avrg_diff_conditions",
    "interaction_condition_time"
    )
  
  # Initialize final return data structures
  correlation_summary <- data.frame(
    category = character(), subcategory = character(),
    n_common = integer(), correlation = numeric(),
    correlation_hits_only = numeric(), n_hits_common = integer(),
    stringsAsFactors = FALSE
  )
  plot_list <- list()
  hits_summary <- data.frame(
    category = character(), subcategory = character(),
    hits_1 = integer(), hits_2 = integer(),
    n_overlap = integer(), jaccard_index = numeric(),
    stringsAsFactors = FALSE
  )
  
  log10_thresh1 <- -log10(adj_p_tresh1 + 1e-16)
  log10_thresh2 <- -log10(adj_p_tresh2 + 1e-16)
  
  for (cat in categories) {    # loop over the limma result categories
    list1 <- splineomics1[["limma_splines_result"]][[cat]]
    list2 <- splineomics2[["limma_splines_result"]][[cat]]
    
    if (!is.list(list1) || !is.list(list2)) {
      message(sprintf("Category '%s' missing or not a list", cat))
      next
    }
    
    shared_names <- intersect(
      names(list1),
      names(list2)
      )
    if (length(shared_names) == 0) {
      message(sprintf("No matching subcategories in '%s'", cat))
      next
    }
    
    for (subcat in shared_names) {    # loop over the conditions/groups
      merged <- merge_limma_results(
        list1[[subcat]],
        list2[[subcat]]
        )
      if (is.null(merged)) {
        message(sprintf("Invalid or no overlap in '%s' > '%s'", cat, subcat))
        next
      }
      
      # Full correlation ----
      corr <- cor(
        merged$adj.P.Val_1,
        merged$adj.P.Val_2,
        method = "spearman",
        use = "complete.obs"
        )
      
      # Main correlation summary entry
      correlation_summary <- rbind(correlation_summary, data.frame(
        category = cat, subcategory = subcat,
        n_common = nrow(merged),
        correlation = round(corr, 3),
        correlation_hits_only = NA,
        n_hits_common = NA
      ))
      
      # Main plot
      plot_list[[paste0(cat, "_", subcat)]] <- plot_pval_correlation(
        df = merged,
        title = paste0(cat, " > ", subcat),
        subtitle = paste("Spearman rho:", round(corr, 3), "(all features)"),
        log10_thresh1 = log10_thresh1,
        log10_thresh2 = log10_thresh2,
        splineomics1_description = splineomics1_description,
        splineomics2_description = splineomics2_description,
        point_color = "black"
      )
      
      # Hits only ----
      sig_hits <- get_hits_only(
        merged,
        adj_p_tresh1,
        adj_p_tresh2
      )
      
      hit_info <- compute_hit_overlap(
        sig_hits,
        adj_p_tresh1,
        adj_p_tresh2
      )
      hits_summary <- rbind(hits_summary, data.frame(
        category = cat, subcategory = subcat,
        hits_1 = length(hit_info$sig1),
        hits_2 = length(hit_info$sig2),
        n_overlap = length(hit_info$overlap),
        jaccard_index = round(hit_info$jaccard, 3)
      ))

      if (nrow(sig_hits) >= 3) {
        corr_hits_only <- cor(
          sig_hits$adj.P.Val_1,
          sig_hits$adj.P.Val_2,
          method = "spearman",
          use = "complete.obs"
          )
        correlation_summary$correlation_hits_only[nrow(correlation_summary)] <-
          round(corr_hits_only, 3)
        correlation_summary$n_hits_common[nrow(correlation_summary)] <-
          nrow(sig_hits)
        
        plot_list[[paste0(cat, "_", subcat, "_hits_only")]] <- 
          plot_pval_correlation(
            df = sig_hits,
            title = paste0(cat, " > ", subcat, " (significant only)"),
            subtitle = paste("Spearman rho:", round(corr_hits_only, 3)),
            log10_thresh1 = log10_thresh1,
            log10_thresh2 = log10_thresh2,
            splineomics1_description = splineomics1_description,
            splineomics2_description = splineomics2_description,
            point_color = "blue"
        )
      }
    }
  }
  
  return(list(
    correlation_summary = correlation_summary,
    plots = plot_list,
    hits_summary = hits_summary
  ))
}


# Level 1 function definitions -------------------------------------------------


#' Validate input arguments for compare_results
#'
#' @noRd
#'
#' @description
#' Checks that input descriptions are single-character strings, p-value
#' thresholds are numeric values between 0 and 1 (inclusive), and both input
#' SplineOmics objects are of S3 class "SplineOmics".
#'
#' @param splineomics1 An S3 object of class "SplineOmics".
#' @param splineomics2 An S3 object of class "SplineOmics".
#' @param splineomics1_description Character. Label for first analysis 
#'                                 (length 1).
#' @param splineomics2_description Character. Label for second analysis 
#'                                 (length 1).
#' @param adj_p_tresh1 Numeric. P-value threshold between 0 and 1 (inclusive).
#' @param adj_p_tresh2 Numeric. P-value threshold between 0 and 1 (inclusive).
#'
#' @return NULL. Throws an error if any check fails.
#' 
validate_compare_inputs <- function(
    splineomics1,
    splineomics2,
    splineomics1_description,
    splineomics2_description,
    adj_p_tresh1,
    adj_p_tresh2
) {
  if (!inherits(splineomics1, "SplineOmics")) {
    stop("`splineomics1` must be an S3 object of class 'SplineOmics'")
  }
  if (!inherits(splineomics2, "SplineOmics")) {
    stop("`splineomics2` must be an S3 object of class 'SplineOmics'")
  }
  
  if (!is.character(splineomics1_description) ||
      length(splineomics1_description) != 1) {
    stop("`splineomics1_description` must be a single character string")
  }
  if (!is.character(splineomics2_description) ||
      length(splineomics2_description) != 1) {
    stop("`splineomics2_description` must be a single character string")
  }
  
  if (!is.numeric(adj_p_tresh1) || length(adj_p_tresh1) != 1 ||
      adj_p_tresh1 < 0 || adj_p_tresh1 > 1) {
    stop("`adj_p_tresh1` must be a numeric value between 0 and 1")
  }
  if (!is.numeric(adj_p_tresh2) || length(adj_p_tresh2) != 1 ||
      adj_p_tresh2 < 0 || adj_p_tresh2 > 1) {
    stop("`adj_p_tresh2` must be a numeric value between 0 and 1")
  }
  
  invisible(NULL)
}


#' Plot correlation of -log10 adjusted p-values with thresholds
#'
#' @noRd
#'
#' @description
#' Generates a ggplot showing correlation between two sets of -log10 adjusted
#' p-values with optional thresholds and axis labels based on descriptions.
#'
#' @param df A data frame with `adj.P.Val_1` and `adj.P.Val_2` columns.
#' @param title Character. Title for the plot.
#' @param subtitle Character. Subtitle for the plot.
#' @param log10_thresh1 Numeric. Threshold line for y-axis (analysis 1).
#' @param log10_thresh2 Numeric. Threshold line for x-axis (analysis 2).
#' @param splineomics1_description Character. Label for x-axis.
#' @param splineomics2_description Character. Label for y-axis.
#' @param point_color Character. Color of plotted points (default: "black").
#'
#' @return A ggplot object showing the correlation with threshold lines.
#' 
plot_pval_correlation <- function(
    df,
    title,
    subtitle,
    log10_thresh1,
    log10_thresh2,
    splineomics1_description,
    splineomics2_description,
    point_color = "black"
) {
  
  ggplot2::ggplot(df, ggplot2::aes(
    x = -log10(adj.P.Val_1 + 1e-16),
    y = -log10(adj.P.Val_2 + 1e-16)
  )) +
    ggplot2::geom_point(alpha = 0.5, color = point_color) +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "gray"
    ) +
    ggplot2::geom_vline(
      xintercept = log10_thresh1,
      linetype = "dashed",
      color = "red"
    ) +
    ggplot2::geom_hline(
      yintercept = log10_thresh2,
      linetype = "dashed",
      color = "red"
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = bquote(-log[10](adj.P.Val[1]) ~ .(splineomics1_description)),
      y = bquote(-log[10](adj.P.Val[2]) ~ .(splineomics2_description))
    ) +
    ggplot2::theme_minimal()
}


#' Merge two limma topTables by feature name
#'
#' @noRd
#'
#' @description
#' Safely merges two limma topTables on `feature_names` if both contain the
#' required columns. Returns NULL if inputs are invalid or have no overlap.
#'
#' @param tt1 A data frame from the first analysis containing `feature_names`
#'   and `adj.P.Val`.
#' @param tt2 A data frame from the second analysis with the same structure.
#'
#' @return A merged data frame with suffixed columns, or NULL if merge fails.
#' 
merge_limma_results <- function(
    tt1,
    tt2
    ) {
  
  if (!is.data.frame(tt1) || !is.data.frame(tt2)) return(NULL)
  if (!all(c("feature_names", "adj.P.Val") %in% names(tt1)) ||
      !all(c("feature_names", "adj.P.Val") %in% names(tt2))) return(NULL)
  
  merged <- merge(
    tt1[, c("feature_names", "adj.P.Val")],
    tt2[, c("feature_names", "adj.P.Val")],
    by = "feature_names",
    suffixes = c("_1", "_2")
  )
  if (nrow(merged) == 0) return(NULL)
  return(merged)
}


#' Compute hit overlap and Jaccard index between two analyses
#'
#' @noRd
#'
#' @description
#' Identifies significant features from each analysis using specified adjusted
#' p-value thresholds, then calculates the overlap and Jaccard index between
#' the two sets of hits.
#'
#' @param merged A data frame with columns `feature_names`, `adj.P.Val_1`, and
#'               `adj.P.Val_2`.
#' @param adj_p_tresh1 Numeric. Significance threshold for the first analysis.
#' @param adj_p_tresh2 Numeric. Significance threshold for the second analysis.
#'
#' @return A list with elements:
#'   \item{sig1}{Character vector of significant features from analysis 1}
#'   \item{sig2}{Character vector of significant features from analysis 2}
#'   \item{overlap}{Character vector of shared significant features}
#'   \item{jaccard}{Numeric Jaccard index of the two hit sets}
#'   
compute_hit_overlap <- function(
    merged,
    adj_p_tresh1, 
    adj_p_tresh2
    ) {
  
  sig1 <- merged$feature_names[merged$adj.P.Val_1 < adj_p_tresh1]
  sig2 <- merged$feature_names[merged$adj.P.Val_2 < adj_p_tresh2]
  overlap <- intersect(sig1, sig2)
  union_set <- union(sig1, sig2)
  jaccard <- if (length(union_set) > 0) 
    length(overlap) / length(union_set) else NA_real_
  
  list(
    sig1 = sig1,
    sig2 = sig2,
    overlap = overlap,
    jaccard = jaccard
    )
}


#' Filter significant features from merged p-value table
#'
#' @noRd
#'
#' @description
#' Returns a subset of the merged data frame containing only rows where
#' either adjusted p-value (from splineomics1 or splineomics2) is below
#' the specified threshold.
#'
#' @param merged A data frame with columns `adj.P.Val_1` and `adj.P.Val_2`
#'   representing adjusted p-values from two analyses.
#' @param adj_p_tresh1 Numeric. Significance threshold for the first analysis.
#' @param adj_p_tresh2 Numeric. Significance threshold for the second analysis.
#'
#' @return A filtered data frame with only significant features 
#'         (in either analysis).
#'
get_hits_only <- function(
    merged, 
    adj_p_tresh1,
    adj_p_tresh2
    ) {
  
  merged[
    merged$adj.P.Val_1 < adj_p_tresh1 
    | merged$adj.P.Val_2 < adj_p_tresh2, ]
}
