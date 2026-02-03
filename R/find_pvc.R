#' Find peaks and valleys in time-series omics data
#'
#' @description
#' Identifies significant local peaks or valleys (excursions) in time-series
#' omics data using a Union-Intersection Test (UIT)-based approach. This
#' function performs statistical detection only and returns per-condition
#' excursion statistics. Plotting and report generation are handled by
#' `create_pvc_report()`.
#'
#' @param splineomics `list`: Preprocessed time-series input. Must include
#'   at least:
#' \itemize{
#'   \item \code{data}: `matrix` Feature-by-sample numeric matrix.
#'   \item \code{meta}: `data.frame` Sample metadata for \code{data} columns.
#'     Must include a \code{"Time"} column and the column specified by
#'     \code{condition}.
#'   \item \code{condition}: `character(1)` Column name in \code{meta} that
#'     defines condition levels.
#'   \item \code{meta_batch_column}: `character(1)` Column name in \code{meta}
#'     identifying replicates or batches.
#'   \item \code{meta_batch2_column}: `character(1)` Optional second batch
#'     column name in \code{meta}.
#' }
#'
#' @param alphas `numeric(1)` | `list(numeric)`: Significance threshold(s) for
#'   excursion calls. A scalar applies to all condition levels. A named list
#'   provides one threshold per condition level; names must match levels in
#'   \code{meta[[condition]]}. Internally normalized for per-level access.
#'
#' @param padjust_method `character(1)`: Multiple testing correction method
#'   passed to `pvc_test()`. Defaults to `"BH"`.
#'
#' @param support `numeric(1)`: Minimum number of non-NA values required per
#'   timepoint neighborhood to keep a PVC p-value. For a given feature at a
#'   target timepoint, the target and its neighboring timepoints must each
#'   have at least \code{support} non-NA observations; otherwise the p-value
#'   at the target timepoint is set to \code{NA}.
#'   
#' @param verbose `logical(1)`: Boolean flag controlling the display of 
#' messages.
#'
#' @return A named list by condition level. Each element contains:
#' \describe{
#'   \item{\code{alpha}}{`numeric(1)` The threshold used for that level.}
#'   \item{\code{pvc_adj_pvals}}{`matrix` Adjusted PVC p-values per feature and
#'   timepoint (after support filtering).}
#'   \item{\code{pvc_pattern_summary}}{`data.frame` Counts of excursion labels
#'   (\code{p}, \code{v}, \code{b}, \code{t}) per timepoint.}
#' }
#' Attributes on the returned list include \code{padjust_method},
#' \code{support}, and \code{batch_effects}.
#'
#' @details
#' A peak or valley is a timepoint whose value is significantly different
#' from both neighbors and deviates in the same direction: higher than both
#' (peak) or lower than both (valley). This is tested via a compound limma
#' contrast: \eqn{(T - T_{prev}) + (T - T_{next}) = 2T - T_{prev} - T_{next}}.
#' The resulting p-values are adjusted and compared to the per-level
#' \code{alpha}.
#'
#' @examples
#' set.seed(1)
#'
#' toy6 <- matrix(
#'   c(
#'     3, 5, 8, 12, 17, 23,
#'     23, 17, 13, 9, 6, 4,
#'     5, 3, 2, 2, 3, 5,
#'     1, 4, 9, 8, 4, 1,
#'     10, 10, 10, 10, 10, 10,
#'     2, 2, 2, 9, 12, 15,
#'     4, 5, 7, 10, 14, 19,
#'     12, 11, 9, 8, 9, 12
#'   ),
#'   nrow = 8,
#'   ncol = 6,
#'   byrow = TRUE,
#'   dimnames = list(paste0("f", 1:8), paste0("s", 1:6))
#' )
#'
#' wt0 <- rowMeans(toy6[, 1:3])
#' ko0 <- rowMeans(toy6[, 4:6])
#'
#' spike_tp <- 3
#' spike_amp <- 3
#'
#' flat4 <- function(base) cbind(base, base, base, base)
#' wt <- flat4(wt0)
#' wt[, spike_tp] <- wt[, spike_tp] + spike_amp
#' ko <- flat4(ko0)
#'
#' rep3 <- function(M, sd = 0.2) {
#'   do.call(cbind, lapply(1:3, function(i) {
#'     M + matrix(rnorm(length(M), sd = sd), nrow(M), ncol(M))
#'   }))
#' }
#'
#' toy_data <- cbind(rep3(wt), rep3(ko))
#' rownames(toy_data) <- rownames(toy6)
#' colnames(toy_data) <- paste0("s", seq_len(ncol(toy_data)))
#'
#' time <- 0:3
#' toy_meta <- data.frame(
#'   Time = rep(time, times = 2 * 3),
#'   condition = rep(c("WT", "KO"), each = 3 * length(time)),
#'   Replicate = rep(paste0("R", 1:3), each = length(time), times = 2),
#'   row.names = colnames(toy_data),
#'   stringsAsFactors = FALSE
#' )
#'
#' splineomics <- list(
#'   data = toy_data,
#'   meta = toy_meta,
#'   condition = "condition",
#'   meta_batch_column = "Replicate"
#' )
#'
#' res <- find_pvc(
#'   splineomics = splineomics,
#'   alphas = 0.05,
#'   padjust_method = "BH",
#'   support = 1
#' )
#'
#' @export
#'
find_pvc <- function(
    splineomics,
    alphas = 0.05,
    padjust_method = "BH",
    support = 1,
    verbose = FALSE
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
    meta_batch2_column <- splineomics[["meta_batch2_column"]]

    batch_effects <- c(meta_batch_column, meta_batch2_column)

    # Get all unique condition levels
    condition_levels <- unique(meta[[condition]])
    alphas <- normalize_alphas(
        alphas = alphas,
        condition_levels = condition_levels
    )

    # Initialize output list
    results <- list()

    for (level_index in seq_along(condition_levels)) {
        level <- condition_levels[[level_index]]

        # Select samples for this condition level
        selected_samples <- meta[[condition]] == level
        sub_data <- data[, selected_samples, drop = FALSE]
        sub_meta <- meta[selected_samples, , drop = FALSE]

        alpha <- alphas[[level]]

        # Run peak/valley detection on the subset
        pvc_pvals <- pvc_test(
            data = sub_data,
            meta = sub_meta,
            batch_effects = batch_effects,
            padjust_method = padjust_method
        )

        pvc_pvals <- filter_pvals_by_support(
            pvc_pvals,
            sub_data,
            sub_meta,
            support
        )

        labels <- classify_excursions(
            data = sub_data,
            meta = sub_meta,
            pvc_pvals = pvc_pvals,
            alpha = alpha
        )

        # Count each label type per column (i.e., per timepoint)
        pattern_counts <- apply(
            labels,
            2,
            function(col) table(factor(col, levels = c("p", "v", "b", "t")))
        )

        # Convert to data frame for easier aggregation
        pattern_df <- as.data.frame(pattern_counts)
        pattern_summary <- rowSums(pattern_df)
        total_hits <- sum(pattern_summary)
        
        if(isTRUE(verbose)) {
            message_pvc_pattern_hits(
                level = level,
                total_hits = total_hits,
                pattern_summary = pattern_summary,
                pattern_df = pattern_df
            )
        }
        results[[as.character(level)]][["pvc_adj_pvals"]] <- pvc_pvals
        results[[as.character(level)]][["pvc_pattern_summary"]] <- pattern_df
    }

    results
}


# Level 1 function definitions -------------------------------------------------


#' Normalize alpha thresholds for multiple condition levels
#'
#' @noRd
#'
#' @description
#' This helper function standardizes the `alphas` input to ensure it can be
#' consistently used across multiple condition levels. The input can either be
#' a single numeric value, which will be replicated and named for each
#'  condition
#' level, or a named list containing a numeric threshold for each condition
#' level.
#'
#' @param alphas Either a single numeric value (e.g., `0.05`) or a named list
#' of numeric values. If a single value is provided, it is applied to all
#' `condition_levels`. If a list is provided, it must be named, with each name
#' corresponding to a value in `condition_levels`.
#'
#' @param condition_levels A character vector of condition levels for which
#' alpha thresholds are required. These must match the names in the `alphas`
#' list if a list is provided.
#'
#' @return A named list of alpha thresholds, one per condition level.
#'
normalize_alphas <- function(
    alphas,
    condition_levels) {
    if (is.numeric(alphas) && length(alphas) == 1) {
        alphas <- setNames(
            as.list(rep(
                alphas,
                length(condition_levels)
            )),
            condition_levels
        )
    } else if (is.list(alphas)) {
        if (!all(condition_levels %in% names(alphas))) {
            stop_call_false(
                "Alpha list must be named for all condition levels."
                )
        }
    } else {
        stop_call_false(
            "Invalid alphas: must be a single numeric or a named list."
            )
    }
    return(alphas)
}


#' Detect peaks/valleys in time-series omics using compound contrasts in limma
#'
#' @noRd
#'
#' @description
#' This function identifies excursions in time-series omics data using a
#' rolling window Union-Intersection Test (UIT). A timepoint is flagged as
#' an excursion if it is significantly different from both its adjacent
#' neighbors in the same direction (higher or lower). The significance is
#' determined using limma's moderated t-test with FDR correction.
#'
#' @param data A numeric matrix, where rows correspond to features
#'  (e.g., genes,
#' proteins, metabolites) and columns correspond to samples.
#' @param meta A data frame containing metadata for the samples. Must include
#' a column named `"Time"` that specifies the timepoint for each sample.
#' @param batch_effects Vector of 2 strings describing the names of the meta
#' batch columns.
#' @param padjust_method Method to correct the p-values for multiple hypothesis
#'                       testing.
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
pvc_test <- function(
    data,
    meta,
    batch_effects,
    padjust_method = "BH") {
    # Extract unique timepoints in sorted order
    unique_timepoints <- sort(unique(meta$Time))
    num_timepoints <- length(unique_timepoints)

    # Ensure there are at least 3 timepoints for meaningful comparisons
    if (num_timepoints < 3) {
        stop("Not enough timepoints for compound-contrast-based comparisons.")
    }

    # Remove the batch effect(s)
    if (!identical(batch_effects, FALSE)) {
        if (length(batch_effects) == 1) {
            data <- limma::removeBatchEffect(
                data,
                batch = meta[[batch_effects[1]]]
            )
        } else if (length(batch_effects) == 2) {
            data <- limma::removeBatchEffect(
                data,
                batch = meta[[batch_effects[1]]],
                batch2 = meta[[batch_effects[2]]]
            )
        }
    }

    valid_timepoints <- make.names(as.character(unique_timepoints))

    time <- factor(meta$Time, levels = unique_timepoints)
    design <- stats::model.matrix(~ 0 + time)
    colnames(design) <- valid_timepoints

    fit <- limma::lmFit( # Fit base limma model
        data,
        design
    )

    # Initialize matrix to store compound contrast p-values
    peak_valley_pvals <- matrix(
        NA,
        nrow = nrow(data),
        ncol = num_timepoints - 2
    )
    rownames(peak_valley_pvals) <- rownames(data)
    colnames(peak_valley_pvals) <- paste0(
        unique_timepoints[-c(1, num_timepoints)],
        "_vs_neighbors"
    )

    # Loop over internal timepoints to compute compound contrast p-values
    for (t in 2:(num_timepoints - 1)) {
        compound_contrast <- paste0(
            "(", valid_timepoints[t], "-", valid_timepoints[t - 1], ") + ",
            "(", valid_timepoints[t], "-", valid_timepoints[t + 1], ")"
        )

        contrast_matrix <- limma::makeContrasts(
            contrasts = compound_contrast,
            levels = design
        )

        fit2 <- limma::contrasts.fit(fit, contrast_matrix)
        fit2 <- limma::eBayes(fit2)

        peak_valley_pvals[, t - 1] <- limma::topTable(
            fit2,
            coef = 1,
            number = Inf,
            adjust.method = padjust_method,
            sort.by = "none"
        )[, "adj.P.Val"]
    }

    return(peak_valley_pvals)
}


#' Filter p-values Based on Support Across Neighboring Timepoints
#'
#' @noRd
#'
#' @description
#' This function filters a p-value matrix by evaluating the data support
#' for each feature (row) and timepoint (column). A p-value is considered
#' valid only if the corresponding feature has at least a specified number
#' of non-NA values (`support`) across the current, previous, and next
#' timepoints.
#'
#' @details
#' Timepoint information is taken from the `Time` column in `sub_meta`,
#' which must be numeric. For each column in `pvc_pvals`, the function:
#'
#' * Extracts the timepoint from the column name (before the first underscore).
#' * Identifies the previous and next unique timepoints from `sub_meta$Time`.
#' * Finds all columns in `sub_data` associated with each of these timepoints.
#' * Checks if the feature has at least `support` non-NA values in each group.
#'
#' If any of the three groups (previous, current, next) fails the support
#' threshold, the corresponding p-value is replaced with `1`.
#'
#' Timepoints at the beginning or end of the time series (i.e., with no
#' previous or next timepoint) are treated as invalid and all their p-values
#' are set to `1`.
#'
#' @param pvc_pvals A numeric matrix of p-values. Each column corresponds to
#' a specific timepoint (with replicates), and each row to a feature.
#' @param sub_data A numeric matrix of feature data. Rows are features and
#' columns are samples. Column order must align with rows in `sub_meta`.
#' @param sub_meta A `data.frame` containing metadata for `sub_data` columns.
#' Must include a numeric `Time` column indicating sample timepoints.
#' @param support A non-negative integer specifying the minimum number of
#' non-NA values required for each of the three timepoint groups.
#'
#' @return A modified version of `pvc_pvals` where p-values are replaced with
#' `1` if support criteria are not met.
#'
filter_pvals_by_support <- function(
    pvc_pvals,
    sub_data,
    sub_meta,
    support) {
    # Ensure Time is numeric and sorted
    time_values <- sub_meta$Time
    if (!is.numeric(time_values)) {
        stop("sub_meta$Time must be numeric for sorting and comparison.")
    }

    # Get unique sorted timepoints
    unique_times <- sort(unique(time_values))

    # Create a lookup for each timepoint -> column indices in sub_data
    time_to_indices <- lapply(unique_times, function(t) which(time_values == t))
    names(time_to_indices) <- as.character(unique_times)

    # Initialize result
    filtered_pvals <- pvc_pvals

    # Loop over each column in pvc_pvals (each timepoint)
    for (col_idx in seq_len(ncol(pvc_pvals))) {
        # Extract timepoint t from column name (before first underscore)
        col_name <- colnames(pvc_pvals)[col_idx]
        t <- as.numeric(sub("_.*", "", col_name))
        if (is.na(t)) next # skip if invalid timepoint

        # Identify previous and next timepoints
        t_pos <- which(unique_times == t)
        if (length(t_pos) == 0) next # timepoint not found

        # Skip if t is the first or last (no prev/next)
        if (t_pos == 1 || t_pos == length(unique_times)) {
            filtered_pvals[, col_idx] <- 1
            next
        }

        t_prev <- unique_times[t_pos - 1]
        t_next <- unique_times[t_pos + 1]

        # Get column indices in sub_data for each group
        idx_prev <- time_to_indices[[as.character(t_prev)]]
        idx_curr <- time_to_indices[[as.character(t)]]
        idx_next <- time_to_indices[[as.character(t_next)]]

        # Combine data per feature (i.e., per row)
        for (row_idx in seq_len(nrow(pvc_pvals))) {
            vals_prev <- sub_data[row_idx, idx_prev]
            vals_curr <- sub_data[row_idx, idx_curr]
            vals_next <- sub_data[row_idx, idx_next]

            non_na_prev <- sum(!is.na(vals_prev))
            non_na_curr <- sum(!is.na(vals_curr))
            non_na_next <- sum(!is.na(vals_next))

            if (non_na_prev < support ||
                non_na_curr < support ||
                non_na_next < support) {
                filtered_pvals[row_idx, col_idx] <- 1
            }
        }
    }

    return(filtered_pvals)
}


#' Classify Peaks, Valleys, and Cliffs from Compound Contrast P-values
#'
#' @noRd
#'
#' @description
#' Assigns a label to each internal timepoint (T2 to Tn-1) indicating
#' whether it is a peak (\code{"p"}), valley (\code{"v"}), top of a
#' cliff (\code{"t"}), or bottom of a cliff (\code{"b"}). Timepoints
#' that are not statistically significant or cannot be classified are
#' labeled as \code{"."}.
#'
#' The classification is based on compound contrast p-values from the
#' \code{peak_valley_test()} function, as well as the relative symmetry
#' of changes between adjacent timepoints.
#'
#' @param data A numeric matrix or data frame of expression values, with
#' rows as features and columns matching \code{meta$Time}.
#' @param meta A data frame containing metadata with at least a column
#' \code{Time} indicating timepoint assignment.
#' @param peak_valley_pvals A matrix of adjusted p-values from a compound
#' contrast test, with rows matching \code{data} and columns corresponding
#' to internal timepoints (i.e., excluding first and last timepoints).
#' @param alpha Significance threshold for p-value filtering. Only
#' timepoints with \code{p < alpha} are considered for classification.
#' @param symmetry_ratio Minimum ratio of the smaller to larger slope
#' required to consider a pattern symmetric enough to be a peak or valley
#' (default: 0.3). Values below this threshold are classified as cliffs.
#'
#' @return A data frame with the same number of rows as \code{data} and
#' one column per timepoint. Each entry is a character label: \code{"p"},
#' \code{"v"}, \code{"t"}, \code{"b"}, or \code{"."} for no signal.
#'
classify_excursions <- function(
    data,
    meta,
    pvc_pvals,
    alpha = 0.05,
    symmetry_ratio = 0.3) {
    # Extract sorted timepoints
    unique_timepoints <- sort(unique(meta$Time))
    num_timepoints <- length(unique_timepoints)

    # Initialize label matrix with "."
    peak_valley_labels <- matrix(
        ".",
        nrow = nrow(data),
        ncol = num_timepoints
    )
    rownames(peak_valley_labels) <- rownames(data)
    colnames(peak_valley_labels) <- unique_timepoints

    # Classification logic for each triplet
    for (i in seq_len(nrow(data))) {
        for (t in 2:(num_timepoints - 1)) {
            p_val <- pvc_pvals[i, t - 1]

            # Skip if NA or not significant
            if (is.na(p_val) || p_val >= alpha) next

            # Extract timepoint means
            prev_mean <- mean(
                data[i, which(meta$Time == unique_timepoints[t - 1])],
                na.rm = TRUE
            )
            curr_mean <- mean(
                data[i, which(meta$Time == unique_timepoints[t])],
                na.rm = TRUE
            )
            next_mean <- mean(
                data[i, which(meta$Time == unique_timepoints[t + 1])],
                na.rm = TRUE
            )

            prev_change <- curr_mean - prev_mean
            next_change <- curr_mean - next_mean

            if (!is.na(prev_change) && !is.na(next_change)) {
                same_direction <- sign(prev_change) == sign(next_change)
                if (same_direction) {
                    abs_prev <- abs(prev_change)
                    abs_next <- abs(next_change)
                    smaller <- min(abs_prev, abs_next)
                    larger <- max(abs_prev, abs_next)

                    if (smaller / larger >= symmetry_ratio) {
                        if (prev_change > 0) {
                            peak_valley_labels[i, t] <- "p"
                        } else {
                            peak_valley_labels[i, t] <- "v"
                        }
                    } else {
                        peak_valley_labels[i, t] <- classify_cliff(
                            prev_change,
                            next_change
                        )
                    }
                } else {
                    peak_valley_labels[i, t] <- classify_cliff(
                        prev_change,
                        next_change
                    )
                }
            }
        }
    }

    # Return as data.frame for compatibility
    results_df <- data.frame(
        peak_valley_labels,
        check.names = FALSE
    )

    return(results_df)
}


#' Emit a concise PVC excursion summary message for one condition level
#'
#' @param level character(1)
#' @param total_hits integer(1)
#' @param pattern_summary numeric
#' @param pattern_df data.frame
#'
#' @return Invisibly returns NULL.
#' 
#' @noRd
#' 
message_pvc_pattern_hits <- function(
        level,
        total_hits,
        pattern_summary,
        pattern_df
) {
    pattern_type_line <- paste(
        names(pattern_summary),
        pattern_summary,
        sep = ": ",
        collapse = ", "
    )
    
    per_timepoint_lines <- paste(
        colnames(pattern_df),
        apply(
            pattern_df,
            2,
            function(x) {
                paste(
                    names(x),
                    x,
                    sep = "=",
                    collapse = "; "
                )
            }
        ),
        sep = ": ",
        collapse = "\n"
    )
    
    message(
        "\nDetected ",
        total_hits,
        " total pattern hits for condition level: ",
        level,
        "\n\n",
        "Summary by pattern type:\n",
        pattern_type_line,
        "\n\n",
        "Breakdown by timepoint:\n",
        per_timepoint_lines,
        "\n"
    )
    
    invisible(NULL)
}


# Level 2 function definitions -------------------------------------------------


#' Classify Cliff Pattern as Top or Bottom
#'
#' @noRd
#'
#' @description
#' Determines whether a significant excursion that is not a peak or valley
#' represents the top or bottom of a cliff, based on the direction and
#' magnitude of changes before and after a timepoint.
#'
#' @param prev_change Numeric value representing the change from the previous
#'                    timepoint to the current timepoint.
#' @param next_change Numeric value representing the change from the current
#'                    timepoint to the next timepoint.
#'
#' @return A single character: \code{"t"} for top of a cliff (upward into or
#' downward out of the timepoint),
#' or \code{"b"} for bottom of a cliff (downward into or upward out of the
#' timepoint).
#'
#' @examples
#' classify_cliff(-2, 0.5) # returns "b" (bottom of a cliff)
#' classify_cliff(1.5, -3) # returns "t" (top of a cliff)
#'
classify_cliff <- function(
    prev_change,
    next_change) {
    abs_prev <- abs(prev_change)
    abs_next <- abs(next_change)

    if (abs_prev > abs_next) {
        return(ifelse(prev_change < 0, "b", "t")) # drop into point -> bottom
    } else {
        return(ifelse(next_change < 0, "t", "b")) # drop after point -> top
    }
}