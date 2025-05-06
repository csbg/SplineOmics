#' Correlate adjusted p-values between two SplineOmics results
#'
#' This function compares the `adj.P.Val` values from two SplineOmics analysis
#' objects by matching features across three categories of topTables:
#' `time_effect`, `avrg_diff_conditions`, and `interaction_condition_time`.
#' For each matching subcategory, it calculates the Spearman correlation
#' of adjusted p-values (`adj.P.Val`) for shared features and creates
#' scatter plots of `-log10(adj.P.Val)` between the two analyses.
#'
#' @param splineomics1 A SplineOmics result object (typically a list) containing
#'   a `"limma_splines_result"` entry with topTables organized by category and 
#'   subcategory.
#' @param splineomics2 A second SplineOmics result object to compare against the
#'                     first.
#'
#' @return A list with two components:
#' \describe{
#'   \item{correlation_summary}{A data frame summarizing the correlation between
#'     the two objects per category and subcategory. Columns include `category`,
#'     `subcategory`, `n_common` (number of shared features), and
#'     `correlation` (Spearman rho).}
#'   \item{plots}{A named list of ggplot objects visualizing the correlation
#'     of `-log10(adj.P.Val)` values for each matched subcategory.}
#' }
#'
#' @export
#' 
correlate_results <- function(
    splineomics1,
    splineomics2
    ) {
  
  categories <- c(
    "time_effect",
    "avrg_diff_conditions",
    "interaction_condition_time"
    )
  
  plot_list <- list()
  correlation_summary <- data.frame(
    category = character(),
    subcategory = character(),
    n_common = integer(),
    correlation = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (cat in categories) {
    list1 <- splineomics1[["limma_splines_result"]][[cat]]
    list2 <- splineomics2[["limma_splines_result"]][[cat]]
    
    if (!is.list(list1) || !is.list(list2)) {
      message(sprintf(
        "Category '%s' missing or not a list in one of the objects",
        cat
        ))
      next
    }
    
    shared_names <- intersect(
      names(list1),
      names(list2)
      )
    
    for (subcat in shared_names) {
      tt1 <- list1[[subcat]]
      tt2 <- list2[[subcat]]
      
      if (!is.data.frame(tt1) || !is.data.frame(tt2)) {
        message(sprintf(
          "Subcategory '%s' under '%s' is not a data frame",
          subcat,
          cat
          ))
        next
      }
      
      if (!all(c("feature_nr", "adj.P.Val") %in% names(tt1)) ||
          !all(c("feature_nr", "adj.P.Val") %in% names(tt2))) {
        message(sprintf(
          "Missing required columns in '%s' > '%s'",
          cat,
          subcat
          ))
        next
      }
      
      merged <- merge(
        tt1[, c("feature_nr", "adj.P.Val")],
        tt2[, c("feature_nr", "adj.P.Val")],
        by = "feature_nr",
        suffixes = c("_1", "_2")
      )
      
      if (nrow(merged) == 0) {
        message(sprintf(
          "No overlapping features in '%s' > '%s'",
          cat,
          subcat
          ))
        next
      }
      
      corr <- cor(
        merged$adj.P.Val_1,
        merged$adj.P.Val_2,
        method = "spearman",
        use = "complete.obs"
        )
      
      # Add to summary
      correlation_summary <- rbind(
        correlation_summary,
        data.frame(
          category = cat,
          subcategory = subcat,
          n_common = nrow(merged),
          correlation = round(corr, 3)
        )
      )
      
      # Plot
      p <- ggplot2::ggplot(merged, ggplot2::aes(
        x = -log10(adj.P.Val_1 + 1e-16),
        y = -log10(adj.P.Val_2 + 1e-16)
      )) +
        ggplot2::geom_point(alpha = 0.5) +
        ggplot2::geom_abline(
          slope = 1,
          intercept = 0,
          linetype = "dashed",
          color = "gray"
          ) +
        ggplot2::labs(
          title = paste0(cat, " > ", subcat),
          x = expression(-log[10](adj.P.Val[1])),
          y = expression(-log[10](adj.P.Val[2])),
          subtitle = paste(
            "Spearman rho:",
            round(corr, 3)
            )
        ) +
        ggplot2::theme_minimal()
      
      plot_list[[paste0(cat, "_", subcat)]] <- p
    }
    
    # Also report if no subcategories matched
    if (length(shared_names) == 0) {
      message(sprintf(
        "No matching subcategories in '%s'",
        cat
        ))
    }
  }
  
  return(list(
    correlation_summary = correlation_summary,
    plots = plot_list
  ))
}
