#' utils scripts contains shared functions that are used by at least two package
#' functions of the SplineOmics package.

# Level 1 internal functions ---------------------------------------------------


#' Create Progress Bar
#'
#' @noRd
#'
#' @description
#' Creates a progress bar for tracking the progress of an iterable task.
#'
#' @param iterable An iterable object (e.g., list or vector) whose length
#' determines the total number of steps.
#' @param message A message to display with the progress bar
#' (default is "Processing").
#'
#' @return A progress bar object from the 'progress' package.
#'
#' @importFrom progress progress_bar
#'
#' @seealso
#' \code{\link{progress_bar}}
#'
create_progress_bar <- function(
    iterable,
    message = "Processing") {
    # Create and return the progress bar
    pb <- progress::progress_bar$new(
        format = paste(" ", message, " [:bar] :percent :elapsed"),
        total = length(iterable),
        width = 60,
        clear = FALSE
    )

    return(pb)
}


#' Create Design Matrix for Splines
#'
#' @noRd
#'
#' @description
#' This function generates a design matrix using spline parameters and metadata.
#' It accommodates both B-splines and natural cubic splines based on the
#' provided spline type and parameters.
#'
#' @param meta A dataframe containing the metadata, including the Time column.
#' @param spline_params A list containing the spline parameters. This list can
#' include `dof` (degrees of freedom), `spline_type`, and `degree`.
#' @param level_index An integer representing the current level index for which
#' the design matrix is being generated.
#' @param design A character string representing the design formula to be used
#' for generating the model matrix.
#'
#' @return A design matrix constructed using the specified spline parameters and
#' design formula.
#'
#' @importFrom splines bs ns
#' @importFrom stats model.matrix as.formula
#'
design2design_matrix <- function(
    meta,
    spline_params,
    level_index,
    design
    ) {
    args <- list(
        x = meta$Time,
        intercept = FALSE
    ) 

    if (!is.null(spline_params$dof)) {
        args$df <- spline_params$dof[level_index]
    } 

    if (spline_params$spline_type[level_index] == "b") {
        args$degree <- spline_params$degree[level_index]
        meta$X <- do.call(splines::bs, args)
    } else { # natural cubic splines
        meta$X <- do.call(splines::ns, args)
    }

    design_matrix <- stats::model.matrix(
        stats::as.formula(design),
        data = meta
    )

    result <- list(
        design_matrix = design_matrix,
        meta = meta
    )
}


#' Merge Annotation with a Single Top Table
#'
#' @noRd
#'
#' @description
#' This function merges annotation information into a single `top_table`
#' dataframe based on the `feature_nr` column.
#'
#' @param top_table A dataframe containing the `top_table` with a `feature_nr`
#'                  column.
#' @param annotation A dataframe containing the annotation information.
#'
#' @return A dataframe with updated `top_table` containing merged annotation
#'         information.
#'
merge_top_table_with_annotation <- function(
    top_table,
    annotation) {
    if (!"feature_nr" %in% names(top_table)) {
        top_table$feature_nr <- seq_len(nrow(top_table))
    } else {
        top_table$feature_nr <- as.integer(as.character(top_table$feature_nr))
    }
    annotation_rows <- annotation[top_table$feature_nr, ]
    top_table <- cbind(top_table, annotation_rows)
}


#' Bind Data with Annotation
#'
#' @noRd
#'
#' @description
#' This function converts a matrix to a dataframe, adds row names as the first
#' column,
#' and binds it with annotation data.
#'
#' @param data A matrix containing the numeric data.
#' @param annotation A dataframe containing the annotation information.
#'
#' @return A dataframe with `data` and `annotation` combined, and the row names
#'  of `data`
#' as the first column named `feature_names`.
#'
bind_data_with_annotation <- function(
    data,
    annotation = NULL) {
    data_df <- as.data.frame(data)

    # Add row names as the first column named feature_names
    combined_df <- cbind(
        feature_name = rownames(data_df),
        data_df
    )

    # If annotation is not NULL, check row count and bind with annotation
    if (!is.null(annotation)) {
        if (nrow(data_df) != nrow(annotation)) {
            stop("The number of rows in data and annotation must be the same.",
                call. = FALSE
            )
        }

        # Bind the annotation with the data
        combined_df <- cbind(combined_df, annotation)
    }

    return(combined_df)
}


#' Print Informational Message
#'
#' @noRd
#'
#' @description
#' This function prints a nicely formatted informational message with a green
#'  "Info" label.
#'
#' @param message_prefix A custom message prefix to be displayed before the
#' success message.
#' @param report_dir The directory where the HTML reports are located.
#'
#' @return NULL
#'
print_info_message <- function(
    message_prefix,
    report_dir) {
    # Green color code for "Info"
    green_info <- "\033[32mInfo\033[0m"

    full_message <- paste(
        "\n",
        green_info,
        message_prefix,
        "completed successfully.\n",
        "Your HTML reports are located in the directory: ", report_dir, ".\n",
        "Please note that due to embedded files, the reports might be",
        "flagged as\n",
        "harmful by other software. Rest assured that they provide no harm.\n"
    )

    message(full_message)
}


#' Stop with custom message without call.
#'
#' @noRd
#'
#' @description
#' A helper function that triggers an error with the specified message and
#' suppresses the function call in the error output. This function behaves
#' similarly to the base `stop()` function but automatically concatenates
#' multiple message strings if provided.
#'
#' @param ... One or more character strings specifying the error message.
#'            If multiple strings are provided, they will be concatenated
#'            with a space between them.
#'
#' @return This function does not return a value; it stops execution and
#'         throws an error.
#'
stop_call_false <- function(...) {
    # Concatenate all arguments into a single string
    message_text <- paste(..., sep = " ")

    # Call stop with the concatenated message and call. = FALSE
    stop(message_text, call. = FALSE)
}


#' Extract fixed and random effects from a model formula string.
#'
#' @noRd
#'
#' @description
#' This function processes a model formula string by separating it into
#' fixed effects and random effects. Random effects are substrings in the
#' format "(...)", and fixed effects are the remaining part of the string.
#'
#' @param formula_string A character string representing the model formula.
#'
#' @return A list containing two components:
#'         - `fixed_effects`: The cleaned fixed effects string.
#'         - `random_effects`: The concatenated random effects string.
#'
#' @examples
#' extract_effects("~ 1 + Condition*Time + Plate + (1|Reactor)")
#' # Returns:
#' # $fixed_effects
#' # [1] "~ 1 + Condition*Time + Plate"
#' #
#' # $random_effects
#' # [1] "(1|Reactor)"
extract_effects <- function(formula_string) {
    # Match all substrings in the format "(...)"
    matches <- gregexpr("\\([^)]*\\)", formula_string, perl = TRUE)
    substrings <- regmatches(formula_string, matches)[[1]]

    # If there are random effects, clean them
    if (length(substrings) > 0) {
        # Remove one non-whitespace character before and after each match
        cleaned_random_effects <- vapply(substrings, function(substring) {
            substring <- sub("\\S\\(", "(", substring)
            substring <- sub("\\)\\S", ")", substring)
            return(substring)
        }, FUN.VALUE = character(1))

        # Concatenate random effects into a single string
        random_effects <- paste(cleaned_random_effects, collapse = " ")
    } else {
        random_effects <- ""
    }

    # Remove all random effects from the original string to get fixed effects
    fixed_effects <- formula_string
    if (length(substrings) > 0) {
        for (substring in substrings) {
            fixed_effects <- sub(substring, "", fixed_effects, fixed = TRUE)
        }
    }

    # Clean up any extra '+' or whitespace from fixed effects
    fixed_effects <- gsub("\\s*\\+\\s*$", "", fixed_effects)
    fixed_effects <- trimws(fixed_effects)

    # Return both fixed and random effects
    return(list(
        fixed_effects = fixed_effects,
        random_effects = random_effects
    ))
}


#' Check for Violation of Homoscedasticity in Linear (or Mixed) Model Inputs
#'
#' @noRd
#'
#' @description
#' This internal helper function tests whether the assumption of
#' homoscedasticity (constant variance of residuals) is violated across
#' experimental conditions for each feature (e.g., gene, protein, metabolite).
#'
#' For each feature independently, the function fits either the spline model
#' with limma or in case of random effects with dream. The model residuals are
#' then tested for heteroscedasticity using Levene's test.
#'
#' Levene's test evaluates whether the variance of residuals differs between
#' groups defined by the experimental condition (supplied via the `condition`
#' parameter). Under the null hypothesis, the residual variance is equal
#' across groups (homoscedasticity). A small p-value indicates that the
#' residual variance differs between groups (heteroscedasticity).
#'
#' This function applies Levene's test across all features independently and
#' counts how many features yield a p-value below the specified `p_threshold`.
#' If the fraction of violating features exceeds `fraction_threshold`, the
#' function concludes that the dataset likely violates the homoscedasticity
#' assumption and suggests switching to a more robust modeling strategy.
#'
#' @param data A numeric matrix of expression values (rows = features,
#'             columns = samples). Should already be log-transformed or
#'             otherwise variance-stabilized if required.
#' @param meta A data frame of sample metadata. Must have one row per column
#'             in `data`, and must contain a numeric `Time` column indicating
#'             the timepoint for each sample.
#' @param design A design formula or matrix for the limma analysis.
#' @param design2design_matrix_result The direct output of the internal function
#'                                    design2design_matrix (see source code
#'                                    of that function).
#' @param condition A character string of the column name of meta that contains
#'                  the levels of the experimental condition.
#' @param random_effects Boolean flag specifying if random effects are in the
#'                       design formula.
#' @param data_type String specifying the omics data type ("rna-seq" or
#'                  "other-omics"). Used to determine the recommendation
#'                  message in case of heteroscedasticity.
#' @param p_threshold Numeric. Significance threshold for the Breusch-Pagan
#'                    test per feature (default is 0.05).
#' @param fraction_threshold Numeric. Proportion of features that must violate
#'                           the homoscedasticity assumption (p < p_threshold)
#'                           to consider the whole dataset heteroscedastic
#'                           (default is 0.1, i.e., 10%).
#'
#' @return A list with two elements:
#' \describe{
#'   \item{violation}{A single logical value (`TRUE` or `FALSE`) indicating
#'                    whether a statistically significant
#'                    violation of homoscedasticity was detected across the
#'                    dataset
#'                    (based on the specified `p_threshold` and
#'                    `fraction_threshold`).}
#'   \item{bp_df}{A data frame with one row per feature, containing:
#'     \describe{
#'       \item{pval}{The p-value from Levene's test for each feature.}
#'       \item{max_var_group}{The group (condition) with the highest residual
#'                            variance.}
#'       \item{max_var}{The value of the maximum group variance.}
#'       \item{violation_flag}{Logical: `TRUE` if the feature's p-value is below
#'                             the threshold (i.e. it violates
#'                             homoscedasticity).}
#'     }
#'   }
#' }
#' A summary message about the detected violations and recommended next steps is
#' printed to the console.
#'
#' @importFrom limma lmFit eBayes
#' @importFrom stats residuals
#'
check_homoscedasticity_violation <- function(
    data,
    meta,
    design,
    design2design_matrix_result,
    condition,
    random_effects = FALSE,
    data_type = "other-omics",
    p_threshold = 0.05,
    fraction_threshold = 0.1) {
    message(paste(
        "Running Levene's test both feature wise and sample wise to implicitly",
        "decide whether to use the limma array weights or not."
    ))

    if (random_effects) {
        message(
            "Checking for heteroscedasticity violation with Levene's test in
      the presence of random effects is currently out of order. Hopefully,
      it will be available in future versions of SplineOmics"
        )

        return(list(
            violation = FALSE, # Single Boolean flag
            bp_df = NULL
        ))

        # message("Random effects detected: fitting model with dream()...")
        #
        # colnames(data) <- rownames(meta)   # dream requires this format
        #
        # # Fit dream model
        # fit <- variancePartition::dream(
        #   exprObj = data,
        #   formula = stats::as.formula(design),
        #   data = design2design_matrix_result[["meta"]],
        #   ddf = NULL
        # )
        # fit <- variancePartition::eBayes(fit = fit)
    } else {
        message("No random effects: fitting model with lmFit()...")

        # Fit limma model
        fit <- limma::lmFit(
            object = data,
            design = design2design_matrix_result[["design_matrix"]]
        )
        fit <- limma::eBayes(fit = fit)
    }

    # 2. Extract residuals   (# residuals: features x samples)
    residuals_matrix <- stats::residuals(
        fit,
        y = data
    )

    # Run Levene’s test feature-wise
    featurewise_levene_result <- check_featurewise_heteroscedasticity_levene(
        residuals_matrix = residuals_matrix,
        meta = meta,
        condition = condition,
        p_threshold = p_threshold,
        fraction_threshold = fraction_threshold
    )

    bp_df <- featurewise_levene_result$bp_df
    fraction_violated <- featurewise_levene_result$fraction_violated
    feature_wise_violation <- featurewise_levene_result$feature_wise_violation

    # Run sample-wise Levene's test
    groupwise_variation_flag <- check_samplewise_heteroscedasticity_levene(
        residuals_matrix = residuals_matrix,
        p_threshold = p_threshold
    )

    # Combine decision
    violation <- feature_wise_violation || groupwise_variation_flag

    # Report decision with detailed reasoning
    message("------------------------------------------------------------")
    if (violation) {
        message(
            "\u2757 Heteroscedasticity detected. Proceeding with",
            "robust strategy."
        )

        reason <- switch(paste(
            feature_wise_violation,
            groupwise_variation_flag,
            sep = "-"
        ),
        "TRUE-TRUE"  = "Both feature-wise & sample-wise tests show violation.",
        "TRUE-FALSE" = "Feature-wise test indicated violation.",
        "FALSE-TRUE" = "Sample-wise test indicated violation."
        )
        message(paste("Reason:", reason))

        if (data_type == "rna-seq") {
            message("\u27A1\uFE0F  Using robust RNA-seq strategy:")
            message("voomWithQualityWeights() to downweight noisy samples.")
        } else {
            message("\u27A1\uFE0F  Using robust modeling strategy:")
            message(
                "arrayWeights() and eBayes(robust = TRUE) to stabilize",
                "variance."
                )
        }
    } else {
        message("\u2705 No strong evidence for heteroscedasticity.")
        message("Proceeding WITHOUT using robust strategy.")
    }
    message("------------------------------------------------------------\n")

    return(list(
        violation = violation,
        bp_df = bp_df,
        percent_violated = 100 * fraction_violated
    ))
}


#' Sanitize metadata columns for modeling
#'
#' @noRd
#'
#' @description
#' Prepares a metadata `data.frame` or `DataFrame` for use in model design by:
#' - Applying `make.names()` to all **character and factor** columns.
#' - Leaving numeric and integer columns untouched.
#'
#' @param meta A `data.frame` or `S4Vectors::DataFrame` containing metadata.
#'
#' @return The same object with sanitized character/factor columns.
#'
sanitize_meta <- function(meta) {
    stopifnot(is.data.frame(meta) || inherits(meta, "DataFrame"))

    for (col in names(meta)) {
        col_data <- meta[[col]]

        if (is.character(col_data)) {
            meta[[col]] <- make.names(col_data)
        } else if (is.factor(col_data)) {
            levels(meta[[col]]) <- make.names(levels(col_data))
        }
        # numeric, integer, and logical columns are left unchanged
    }

    return(meta)
}


# Level 2 internal functions ---------------------------------------------------


#' Check feature-wise heteroscedasticity using Levene's test
#'
#' @noRd
#'
#' @description
#' This function performs Levene's test for each feature to assess whether
#' residual variances differ between groups defined by a condition.
#' It returns a detailed summary, a violation flag, and the annotated results.
#'
#' @param residuals_matrix Matrix of residuals (features x samples).
#' @param meta Sample metadata data frame.
#' @param condition The name of the column in `meta` defining the condition
#'  groups.
#' @param p_threshold Significance threshold (after adjustment).
#' @param fraction_threshold Minimum fraction of violating features to flag a
#'  violation.
#' @return A list containing:
#'   \item{bp_df}{Data frame with test results per feature}
#'   \item{fraction_violated}{Proportion of features violating homoscedasticity}
#'   \item{feature_wise_violation}{TRUE if the threshold is exceeded}
#'
#' @importFrom stats var
#' @importFrom car leveneTest
#' @importFrom pbapply pbapply
#'
check_featurewise_heteroscedasticity_levene <- function(
    residuals_matrix,
    meta,
    condition,
    p_threshold = 0.05,
    fraction_threshold = 0.1) {
    if (is.null(condition)) {
        message("No condition. Skipping feature-wise Levene's test.")
        message("-----------------------------------------------------------\n")
        return(list(
            bp_df = NULL,
            fraction_violated = NA_real_,
            feature_wise_violation = FALSE
        ))
    }

    message("Running feature wise Levene's test...")

    Phase <- as.factor(meta[[condition]])

    bp_results <- pbapply::pbapply(residuals_matrix, 1, function(res) {
        group_vars <- tapply(
            res,
            Phase,
            stats::var,
            na.rm = TRUE
        )

        if (any(is.na(group_vars)) ||
            any(group_vars == 0) ||
            length(group_vars) < 2
        ) {
            return(list(
                pval = NA,
                max_var_group = NA,
                max_var = NA
            ))
        }

        pval <- car::leveneTest(res ~ Phase)[["Pr(>F)"]][1]
        max_group <- names(which.max(group_vars))
        max_var <- max(group_vars)

        list(
            pval = pval,
            max_var_group = max_group,
            max_var = max_var
        )
    })

    # Assemble results
    bp_df <- as.data.frame(
        do.call(rbind, bp_results),
        stringsAsFactors = FALSE
    )
    rownames(bp_df) <- rownames(residuals_matrix)

    bp_df$adj_pval <- p.adjust(bp_df$pval, method = "BH")
    bp_df$violation_flag <- !is.na(
        bp_df$adj_pval
        ) & bp_df$adj_pval < p_threshold
    fraction_violated <- mean(bp_df$violation_flag)

    # Global summary
    n_total <- nrow(bp_df)
    n_violating <- sum(bp_df$violation_flag)
    fraction_violated <- n_violating / n_total

    message("\n------------------------------------------------------------")
    message(sprintf(
        "Fraction of features violating homoscedasticity
    (p < %.3f): %.2f%% (%d/%d features)",
        p_threshold,
        100 * fraction_violated,
        n_violating,
        n_total
    ))

    # Breakdown per group
    violating_groups <- bp_df$max_var_group[bp_df$violation_flag]

    if (length(violating_groups) > 0 && any(!is.na(violating_groups))) {
        group_contributions <- table(violating_groups)
        group_fractions <- prop.table(group_contributions)

        message("Distribution of max variance groups among violating features:")
        for (i in seq_along(group_contributions)) {
            group_name <- names(group_contributions)[i]
            n_group <- group_contributions[i]
            frac_group <- 100 * group_fractions[i]

            message(sprintf(
                "- %s: %.2f%% of violating features (%d/%d)",
                group_name,
                frac_group,
                n_group,
                n_violating
            ))
        }
    } else {
        message("No violating features found.")
    }

    message("------------------------------------------------------------\n")

    feature_wise_violation <- fraction_violated >= fraction_threshold

    return(list(
        bp_df = bp_df,
        fraction_violated = fraction_violated,
        feature_wise_violation = feature_wise_violation
    ))
}


#' Check for Heteroscedasticity Across Samples Using Levene's Test
#'
#' @noRd
#'
#' @description
#' This function tests whether the residual variances differ significantly
#' across samples
#' by applying Levene's test. Each sample is treated as a group in the test,
#' and the null hypothesis is that all samples have equal residual variance.
#'
#' The function is useful for detecting global heteroscedasticity that may not
#' be group-dependent, such as technical noise affecting specific samples. It
#' complements feature-wise group-level heteroscedasticity checks.
#'
#' @param residuals_matrix A numeric matrix of residuals with features in rows
#' and samples in columns.
#' @param p_threshold Significance threshold for rejecting the null hypothesis
#'  of equal variances. Defaults to 0.05.
#'
#' @return A logical scalar. \code{TRUE} if Levene's test detects significant
#'  inter-sample variance differences (p < \code{p_threshold}), otherwise
#'   \code{FALSE}.
#'
#' @details
#' The function reshapes the residual matrix into a long format where each
#' sample is treated
#' as a group, and applies \code{\link[car]{leveneTest}} to test for equality of
#'  variances
#' across samples. A message is printed summarizing the test outcome.
#'
#' @importFrom car leveneTest
#'
check_samplewise_heteroscedasticity_levene <- function(
    residuals_matrix,
    p_threshold = 0.05) {
    message(paste(
        "Running Levene's test across samples to detect inter-sample variance",
        "differences..."
    ))

    # Transpose: rows = samples, columns = features
    sample_residuals <- t(residuals_matrix)

    # Stack into long format
    long_data <- data.frame(
        residual = as.vector(sample_residuals),
        sample = factor(
            rep(
                rownames(sample_residuals),
                ncol(sample_residuals)
            )
        )
    )

    # Run Levene’s test across samples
    levene_res <- car::leveneTest(
        residual ~ sample,
        data = long_data
    )
    pval <- levene_res[["Pr(>F)"]][1]

    message(sprintf("Levene's test p-value (sample-level): %.4g", pval))

    if (!is.na(pval) && pval < p_threshold) {
        message(paste(
            "\u2757 Evidence of heteroscedasticity across samples",
            "(violation detected).\n"
        ))
        return(TRUE)
    } else {
        message(
            "\u2705 No strong evidence of inter-sample variance differences.\n"
            )
        return(FALSE)
    }
}