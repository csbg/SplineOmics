#' Perform a limma spline analysis to identify significant features
#'
#' @description
#' This is the core function, which performs a limma spline analysis to identify
#'  significant time-dependent changes in features (e.g., proteins) within an
#' omics time-series dataset. It evaluates features within each condition level
#' and between levels by comparing average differences and interactions
#' between time and condition.
#'
#' @param splineomics `SplineOmics`: An S3 object of class `SplineOmics` that 
#' contains the following elements:
#' \itemize{
#'   \item \code{data}: `matrix` The matrix of the omics dataset, with the
#'    feature
#'   names optionally as row headers.
#'   \item \code{rna_seq_data}: `list` An optional object containing the
#'    preprocessed
#'   RNA-seq data, such as the output from `limma::voom` or a similar
#'    preprocessing pipeline. This must only be provided when the input is
#'    RNA-seq data.
#'   \item \code{meta}: `data.frame` A dataframe containing metadata
#'    corresponding to the
#'   \code{data}, must include a 'Time' column and the column specified by
#'   \code{condition}. The columns of meta describe the meta info, such as
#'   the time and condition, and each row corresponds to a column in data, and
#'   therefore, contains the meta info for that data column. It is important
#'   that meta and data are matched in this way.
#'   \item \code{padjust_method}: `character(1)` Statistical method that is
#'    used for multiple
#'   hypothesis correction. Supported methods include all that are included in
#'   the p.adjust() function in R: "holm", "hochberg", "hommel", "bonferroni",
#'   "BH" or "fdr", "BY", or "none" for no correction. Default of this package
#'   is "BH".
#'   \item \code{design}: `character(1)` A character string representing the
#'    limma design
#'   formula, such as "~ 1 + Phase*Time + Reactor" for an integrated design, or
#'   "~ 1 + Time + Reactor" for an isolated design.
#'   \item\code{dream_params}: `list | NULL` A named list or NULL. When not
#'    NULL, it can contain the following named
#'      elements:
#'      - `dof`: `integer(1)` An integer greater than 1, specifying the
#'       degrees of freedom
#'      for  the dream topTable. If set to 0, then the best dof is automatically
#'      found with the help of leave-one-out-crossvalidation (loocv). The dof
#'      with the lowest error on the loocv is chosen.
#'      - `KenwardRoger`: `logical(1)` A boolean indicating whether to use the
#'       Kenward-Roger
#'      approximation for mixed models.
#'      Note that random effects are now directly specified in the design
#'      formula and not in `dream_params`.
#'   \item \code{mode}: `character(1)` Specifies how the design formula is
#'    constructed: either `"isolated"` or `"integrated"`.
#'
#'   - `"isolated"`: Each level is analyzed independently, using only the
#'     subset of data corresponding to that level. The design formula does
#'     not include the condition variable, since only one condition is
#'     present in each subset.
#'
#'   - `"integrated"`: All levels are analyzed together in a single model,
#'     using the full dataset. The design formula includes the condition
#'     variable (and optionally interaction terms with it) so that results
#'     are estimated jointly across all levels.
#'   \item \code{condition}: `character(1)` A character string specifying the
#'    column name
#'   in \code{meta} used to define groups for analysis. The condition column
#'   contains the levels of the experiment (such as control and treatment).
#'   \item \code{spline_params}: `list` A list of spline parameters used in the
#'   analysis, including:
#'     \itemize{
#'       \item \code{spline_type}: `character(1)`The type of spline (e.g., `"n"`
#'        for natural
#'       splines or `"b"` for B-splines).
#'       \item \code{dof}: `integer(1)` Degrees of freedom for the spline. If
#'        set to `0`,
#'       `SplineOmics` automatically determines the optimal number of degrees
#'       of freedom using leave-one-out cross-validation and selects the value
#'       that yields the best predictive performance.
#'       \item \code{degree}: `integer(1)` Polynomial degree of the spline basis
#'       (for B-splines only).
#'     }
#'   \item \code{use_array_weights}: `logical(1)` Boolean flag indicating if 
#'    the robust fit
#'    strategy to deal with heteroscedasticity should be used or not. If set to
#'    NULL, then this is handeled implicitly based on the result of the Levene
#'    test. If this test is significant for at least 10% of the features,
#'    then the robust strategy is used. The robust strategy uses the function
#'    voomWithQualityWeights for RNA-seq data instead of the normal voom
#'    function. For other, non-count-based data, the function
#'    limma::arrayWeights is used instead, combined with setting the robust
#'    argument to TRUE for the limma::eBayes function. In summary, the strategy
#'    employed by those functions is to downweights samples with higher
#'    variance. This can be neccessary, because linear models have the
#'    assumption of homoscedasticity, which means that the variance is
#'    (approx.) the same across all datapoints where the linear model is fitted.
#'     If this is violated, then the resulting p-values cannot be trusted
#'    (common statistical wisdom).
#'   \item \code{bp_cfg}: `numeric()` A named numeric vector specifying the
#'    parallelization
#'   configuration, with expected names `"n_cores"`: `integer(1)` and
#'    `"blas_threads"`: `integer(1)`.
#'
#'   This controls how many R worker processes (`n_cores`) and how many
#'   BLAS/OpenBLAS threads per process (`blas_threads`) should be used
#'   during parallel computation.
#'
#'   If `bp_cfg` is `NULL`, missing, or any of its required fields is
#'   `NA`, both `n_cores` and `blas_threads` default to `1`. This effectively
#'   disables parallelization and avoids oversubscription of CPU threads.
#' }
#'
#' @param verbose `logical(1)`: Boolean flag controlling the display of
#'  messages.
#'
#' @return The SplineOmics object, updated with a list with three elements:
#'         - `time_effect`: `list` A list of top tables for each level with the
#'          time effect.
#'         - `avrg_diff_conditions`: `list` A list of top tables for each 
#'         comparison between the levels. The comparison is the average
#'          difference of the values.
#'         - `interaction_condition_time`: `list` A list of top tables for each
#'                                         comparison between levels. The
#'                                         comparison is the interaction between
#'                                         the condition and the time.
#'
#' @importFrom purrr partial map map_chr map2
#' @importFrom stats setNames
#' @importFrom utils combn
#'
#' @examples
#' # Toy data: 4 features x 6 samples (two conditions, three time points)
#' toy_data <- matrix(
#'     c(
#'         3, 5, 8, 12, 17, 23, # f1
#'         23, 17, 13, 9, 6, 4, # f2
#'         5, 3, 2, 2, 3, 5, # f3
#'         1, 4, 9, 8, 4, 1, # f4
#'         10, 10, 10, 10, 10, 10, # f5
#'         2, 2, 2, 9, 12, 15, # f6
#'         4, 5, 7, 10, 14, 19, # f7
#'         12, 11, 9, 8, 9, 12 # f8
#'     ),
#'     nrow = 8, ncol = 6, byrow = TRUE,
#'     dimnames = list(paste0("f", 1:8), paste0("s", 1:6))
#' )
#'
#' toy_meta <- data.frame(
#'     Time = c(0, 1, 2, 0, 1, 2),
#'     condition = rep(c("WT", "KO"), each = 3),
#'     Replicate = rep(c("R1", "R2"), each = 3),
#'     row.names = colnames(toy_data),
#'     stringsAsFactors = FALSE
#' )
#'
#' toy_annot <- data.frame(
#'     feature_nr = 1:8,
#'     gene = c("G1", "G2", "G3", "G4"),
#'     stringsAsFactors = FALSE
#' )
#'
#' # Stub limma "top tables" with minimal required fields
#' # (feature_nr + adj.P.Val)
#' tt_wt <- data.frame(feature_nr = 1:4, adj.P.Val = c(0.01, 0.20, 0.04, 0.60))
#' tt_ko <- data.frame(feature_nr = 1:4, adj.P.Val = c(0.50, 0.03, 0.70, 0.02))
#' tt_c2 <- data.frame(feature_nr = 1:4, adj.P.Val = c(0.04, 0.70, 0.80, 0.90))
#' tt_c3 <- data.frame(feature_nr = 1:4, adj.P.Val = c(0.20, 0.90, 0.03, 0.80))
#'
#' design_str <- "~ 1 + Time*condition"
#'
#' # Minimal spline parameters required by spline machinery
#' spline_params <- list(
#'     spline_type = "n", # natural cubic splines
#'     dof = 1L # degrees of freedom for the spline basis
#' )
#'
#' toy_splineomics <- list(
#'     data = toy_data,
#'     meta = toy_meta,
#'     annotation = toy_annot,
#'     report_info = list(
#'         omics_data_type = "RNA-seq",
#'         data_description = "toy example",
#'         data_collection_date = "2025-01-01",
#'         analyst_name = "Example",
#'         contact_info = "example@example.org",
#'         project_name = "ToyProject"
#'     ),
#'     design = design_str,
#'     mode = "integrated",
#'     condition = "condition",
#'     spline_params = spline_params,
#'     meta_batch_column = NULL,
#'     meta_batch2_column = NULL,
#'     limma_splines_result = list(
#'         time_effect                  = list(WT = tt_wt, KO = tt_ko),
#'         avrg_diff_conditions         = tt_c2,
#'         interaction_condition_time   = tt_c3
#'     ),
#'     feature_name_columns = "gene"
#' )
#' class(toy_splineomics) <- "SplineOmics"
#'
#' toy_splineomics <- run_limma_splines(toy_splineomics)
#'
#' @export
#'
run_limma_splines <- function(
    splineomics,
    verbose = TRUE
    ) {
    start_time <- Sys.time()
    check_splineomics_elements(
        splineomics = splineomics,
        func_type = "run_limma_splines"
    )

    args <- lapply(
        as.list(match.call()[-1]),
        eval,
        parent.frame()
    )
    args[["verbose"]] <- verbose
    check_null_elements(args)
    input_control <- InputControl$new(args)
    input_control$auto_validate()

    data <- splineomics[["data"]]
    rna_seq_data <- splineomics[["rna_seq_data"]]
    meta <- sanitize_meta(splineomics[["meta"]])
    spline_params <- splineomics[["spline_params"]]
    padjust_method <- splineomics[["padjust_method"]]
    design <- splineomics[["design"]]
    dream_params <- splineomics[["dream_params"]]
    mode <- splineomics[["mode"]]
    condition <- splineomics[["condition"]]
    if (is.null(rna_seq_data)) {     
        use_array_weights <- splineomics[["use_array_weights"]]
    } else {     # For RNA-seq data, array weights were handled beforehand.
        use_array_weights <- FALSE
        message(
            "Array weights are ignored at this stage for RNA-seq data, ",
            "as they were already incorporated during the preprocessing step."
        )
    }
    
    bp_cfg <- splineomics[["bp_cfg"]]

    # Because at first I enforced that X in the design formula stands for the
    # time and I heavily oriented my code towards that. But then I realised that
    # it is nonsense to encode the time as X, and now it is explicitly "Time"
    # (because meta must contain the exact name "Time" for this respective
    # column).
    design <- gsub(
        "Time",
        "X",
        design
    )

    feature_names <- rownames(data)
    meta[[condition]] <- factor(meta[[condition]])

    if (mode == "isolated") {
        levels <- levels(meta[[condition]])

        # Get hits for level (within level analysis)
        process_level_with_params <- purrr::partial(
            fit_within_condition_isolated, # the function that is called
            spline_params = spline_params,
            data = data,
            rna_seq_data = rna_seq_data,
            meta = meta,
            design = design,
            dream_params = dream_params,
            condition = condition,
            feature_names = feature_names,
            padjust_method = padjust_method,
            mode = mode,
            use_array_weights = use_array_weights,
            bp_cfg = bp_cfg,
            verbose = verbose
        )

        results_nested <- purrr::imap(
            levels,
            process_level_with_params
        )
        results_list <- purrr::flatten(results_nested)

        merged_dof <- Reduce(
            pmax.int,
            lapply(
                results_list,
                \(x) x$spline_params[["dof"]]
            )
        )
        spline_params[["dof"]] <- merged_dof

        isolated_within_level_top_tables <- stats::setNames(
            purrr::map(results_list, "top_table"),
            names(results_list)
        )

        isolated_fits <- stats::setNames(
            purrr::map(results_list, "fit"),
            names(results_list)
        )

        limma_splines_result <- list(
            time_effect = isolated_within_level_top_tables
        )
    } else {   # mode == integrated
        effects <- extract_effects(design)

        if (spline_params[["dof"]][1] == 0) { # auto-dof.
            best_dof <- select_spline_dof_loocv(
                data = data,
                meta = meta,
                spline_params = spline_params,
                level_index = 1,
                fixed_effects = effects[["fixed_effects"]]
            )
            spline_params$dof[1] <- best_dof
        }
        
        design2design_matrix_result <- design2design_matrix(
            meta = meta,
            spline_params = spline_params,
            level_index = 1,
            design = effects[["fixed_effects"]]
        )

        contrast_plan <- build_contrast_plan(
            condition = condition,
            meta = design2design_matrix_result[["meta"]],
            cols = colnames(design2design_matrix_result[["design_matrix"]]),
            dof = spline_params[["dof"]]
        )
        
        # Step 1: Fit the global model once
        fit_obj <- fit_global_model(
            data = data,
            rna_seq_data = rna_seq_data,
            meta = meta,
            design = design,
            design2design_matrix_result = design2design_matrix_result,
            contrast_plan = contrast_plan,
            effects = effects,
            dream_params = dream_params,
            spline_params = spline_params,
            condition = condition,
            padjust_method = padjust_method,
            use_array_weights = use_array_weights,
            bp_cfg = bp_cfg,
            verbose = verbose
        )
        
        limma_splines_result <- create_all_toptables(
            fit_obj = fit_obj,
            condition = condition,
            feature_names = feature_names
        )
    }

    args <- list(
        splineomics = splineomics,
        limma_splines_result = limma_splines_result,
        fit = if (exists("fit_obj")) fit_obj[["fit"]] else isolated_fits,
        spline_params = spline_params,
        meta = meta
    )

    if (!"use_array_weights" %in% names(args) && exists("results_list")) {
        args$use_array_weights <- purrr::map_lgl(
            results_list,
            "use_array_weights"
        ) |>
            purrr::set_names(names(results_list))
    } else if (!"homosc_violation_result" %in% names(splineomics)) {
        args$homosc_violation_result <- if (exists("fit_obj")) {
            fit_obj[["homosc_violation_result"]]
        } else {
            NULL
        }
    }

    if (verbose) {
        end_time <- Sys.time()
        elapsed <- difftime(end_time, start_time, units = "min")
        message(
            sprintf(
                "\033[32mInfo\033[0m Finished limma spline analysis in %.1f ",
                as.numeric(elapsed)
            ),
            "min"
        )
    }

    splineomics <- do.call(
        update_splineomics,
        args
    )
}


# Level 1 internal functions ---------------------------------------------------


#' Within level analysis for the mode: isolated.
#'
#' @noRd
#'
#' @description
#' Processes a single level within a condition, performing limma analysis
#' and generating the top table of results.
#'
#' @param level The level within the condition to process.
#' @param level_index The index of the level within the condition.
#' @param spline_params A list of spline parameters for the analysis.
#' @param data A matrix of data values.
#' @param rna_seq_data An object containing the preprocessed RNA-seq data,
#' such as the output from `limma::voom` or a similar preprocessing pipeline.
#' @param meta A dataframe containing the metadata for data.
#' @param design A design formula or matrix for the limma analysis.
#' @param dream_params A named list or NULL. When not NULL, it must at least
#' contain the named element 'random_effects', which must contain a string that
#' is a formula for the random effects of the mixed models by dream.
#' Additionally, it can contain the named elements dof, which must be a int
#' bigger than 1, which is the degree of freedom for the dream topTable, and
#' the named element KenwardRoger, which must be a bool, specifying whether
#' to use that method or not.
#' @param condition A character string specifying the condition.
#' @param feature_names A non-empty character vector of feature names.
#' @param padjust_method A character string specifying the p-adjustment method.
#' @param mode A character string specifying the mode
#'            ('isolated' or 'integrated').
#' @param use_array_weights Boolean value specifying whether to use array
#' weights for limma or variancePartition.
#' @param bp_cfg A named numeric vector specifying the parallelization
#'   configuration, with expected names `"n_cores"` and `"blas_threads"`.
#'
#'   This controls how many **R worker processes** (`n_cores`) and how many
#'   **BLAS/OpenBLAS threads per process** (`blas_threads`) should be used
#'   during parallel computation.
#'
#'   If `bp_cfg` is `NULL`, missing, or any of its required fields is
#'   `NA`, both `n_cores` and `blas_threads` default to `1`. This effectively
#'   disables parallelization and avoids oversubscription of CPU threads.
#' @param verbose Boolean flag controlling the display of messages.
#'
#' @return A list containing the name of the results and the top table of
#'          results.
#'
#' @importFrom stats relevel
#'
fit_within_condition_isolated <- function(
    level,
    level_index,
    spline_params,
    data,
    rna_seq_data,
    meta,
    design,
    dream_params,
    condition,
    feature_names,
    padjust_method,
    mode,
    use_array_weights,
    bp_cfg,
    verbose
    ) {
    samples <- which(meta[[condition]] == level)
    data_copy <- data[, samples]
    meta_copy <- meta[meta[[condition]] == level, , drop = FALSE]

    if (!is.null(rna_seq_data) && ncol(rna_seq_data$E) != nrow(meta_copy)) {
        stop_call_false(
            "Mismatch detected: rna_seq_data$E has ",
            ncol(rna_seq_data$E),
            " columns, but meta_copy has ",
            nrow(meta_copy),
            " rows. The most likely cause for this is that you selected ",
            "mode == 'isolated', but passed the full data in rna_seq_data. ",
            "For RNA-seq data, you must pass the respective datasets of the ",
            "different conditions individually (for all other omics datasets, ",
            "this is handled implicitly by SplineOmics: it splits up the data ",
            "and meta into the different conditions. However, that is not ",
            "possible with the RNA-seq data objects)."
        )
    }

    result <- process_within_level(
        data = data_copy,
        rna_seq_data = rna_seq_data,
        meta = meta_copy,
        design = design,
        dream_params = dream_params,
        spline_params = spline_params,
        level_index = level_index,
        padjust_method = padjust_method,
        use_array_weights = use_array_weights,
        bp_cfg = bp_cfg,
        verbose = verbose
    )

    top_table <- process_top_table(
        result,
        feature_names
    )

    results_name <- paste(
        condition,
        level,
        sep = "_"
    )

    out <- list(
        top_table = top_table,
        fit = result[["fit"]],
        spline_params = result[["spline_params"]],
        use_array_weights = result[["use_array_weights"]]
    )

    setNames(list(out), results_name)
}


#' Build a unified contrast plan (joint L matrix + mapping)
#'
#' @noRd
#'
#' @description
#' Builds a single joint contrast matrix \code{L_all} (rows = model
#' coefficients, columns = contrasts) and a mapping object
#' \code{map_for_L_all} that groups contrast column names into reporting
#' categories:
#' \itemize{
#'   \item within-level spline time effects (per level)
#'   \item between-level condition main effect differences (per pair)
#'   \item between-level condition-by-time differences (per pair)
#' }
#'
#' The output supports a unified downstream workflow:
#' \itemize{
#'   \item limma: \code{contrasts.fit(fit, L_all)} post hoc
#'   \item dream: pass \code{L_all} as \code{L=} into
#'     \code{variancePartition::dream()}
#' }
#'
#' @param condition Character scalar. Condition or factor name.
#' @param meta data.frame. Sample metadata containing
#'   \code{meta[[condition]]}.
#' @param cols Character vector. Coefficient names used by the fitted
#'   model. Must match \code{colnames(fit$coefficients)} exactly.
#' @param dof Integer. Number of spline basis terms.
#' @param include_within_time Logical. Include within-level time effects.
#' @param include_between_pairs Logical. Include all pairwise contrasts.
#'
#' @return A list with elements:
#' \describe{
#'   \item{L_all}{Numeric matrix. Rows correspond to \code{cols}.}
#'   \item{map_for_L_all}{Nested list mapping categories to contrast names.}
#'   \item{spline_cols}{Character vector of spline basis column names.}
#'   \item{levels_in_condition}{Levels of the condition factor.}
#' }
#'
build_contrast_plan <- function(
        condition,
        meta,
        cols,
        dof
) {
    dof <- as.integer(dof)
    
    levels_in_condition <- levels(factor(meta[[condition]]))
    if (length(levels_in_condition) < 2L) {
        stop("Condition must have at least two levels.")
    }
    
    spline_cols <- grep("^X$|^X\\d+$", cols, value = TRUE)
    if (length(spline_cols) < 1L) {
        stop("No spline basis columns found in 'cols'.")
    }
    
    if (length(spline_cols) != dof) {
        stop(
            "Mismatch between dof (", dof,
            ") and detected spline columns (",
            length(spline_cols), ")."
        )
    }
    
    main_terms <- paste0(condition, levels_in_condition)
    idx_main <- match(main_terms, cols)
    missing_main <- is.na(idx_main)
    
    if (sum(missing_main) == 1L) {
        baseline <- levels_in_condition[missing_main]
    } else {
        stop(
            "Could not infer baseline of design matrix."
        )
    }
    
    dummy_design_matrix <- matrix(
        0,
        nrow = 1L,
        ncol = length(cols),
        dimnames = list(NULL, cols)
    )
    
    L_blocks <- list()
    
    map_for_L_all <- list(
        within_time = list(),
        between_condition_only = list(),
        between_main_and_time = list(),
        between_condition_time = list()
    )
    
    for (lev in levels_in_condition) {
        L_lev <- build_spline_contrast(
            lev = lev,
            baseline = baseline,
            condition = condition,
            spline_terms = spline_cols,
            design_cols = cols,
            dof = dof
        )
        
        if (!identical(rownames(L_lev), cols)) {
            stop(
                "Spline contrast rownames do not match 'cols' ",
                "for level '", lev, "'."
            )
        }
        
        lev_key <- as.character(lev)
        
        new_names <- paste0(
            "within_",
            as.character(condition),
            "_",
            lev_key,
            "__",
            spline_cols
        )
        
        colnames(L_lev) <- new_names
        
        L_blocks[[paste0("within_", lev_key)]] <- L_lev
        map_for_L_all$within_time[[lev_key]] <- new_names
    }
    
    pairs <- utils::combn(levels_in_condition, 2, simplify = FALSE)
    
    for (pair in pairs) {
        lev1 <- pair[1]
        lev2 <- pair[2]
        
        cm <- get_condition_contrast_coefs(
            condition = condition,
            level_pair = c(lev1, lev2),
            design_matrix = dummy_design_matrix
        )
        
        L_pair <- cm$L
        
        if (!identical(rownames(L_pair), cols)) {
            stop(
                "Pairwise contrast rownames do not match 'cols' ",
                "for pair '", lev1, "' vs '", lev2, "'."
            )
        }
        
        key <- paste0(as.character(lev1), "_vs_", as.character(lev2))
        
        prefix <- paste0(
            "pair_",
            as.character(condition),
            "_",
            key,
            "__"
        )
        
        colnames(L_pair) <- paste0(prefix, colnames(L_pair))
        
        L_blocks[[paste0("pair_", key)]] <- L_pair
        
        map_for_L_all$between_condition_only[[key]] <-
            paste0(prefix, "diff_main")
        
        map_for_L_all$between_main_and_time[[key]] <-
            c(
                paste0(prefix, "diff_main"),
                paste0(prefix, paste0("diff_", spline_cols))
            )
        
        map_for_L_all$between_condition_time[[key]] <-
            paste0(prefix, paste0("diff_", spline_cols))
    }
    
    if (length(L_blocks) < 1L) {
        stop("No contrast blocks were constructed.")
    }
    
    L_all <- do.call(cbind, L_blocks)
    L_all <- as.matrix(L_all)
    storage.mode(L_all) <- "double"
    
    if (!identical(rownames(L_all), cols)) {
        stop("Internal error: L_all rownames do not match 'cols'.")
    }
    
    if (anyDuplicated(colnames(L_all)) > 0L) {
        stop("Internal error: L_all has duplicated column names.")
    }
    
    mapped <- unique(unlist(map_for_L_all, use.names = FALSE))
    missing <- setdiff(mapped, colnames(L_all))
    
    if (length(missing) > 0L) {
        stop(
            "Internal error: map_for_L_all contains names not in L_all: ",
            paste(missing, collapse = ", ")
        )
    }
    
    list(
        L_all = L_all,
        map_for_L_all = map_for_L_all,
        spline_cols = spline_cols,
        levels_in_condition = levels_in_condition,
        baseline = baseline
    )
}


#' Between Level Analysis
#'
#' @noRd
#'
#' @description
#' Performs a between-level analysis using limma to compare specified levels
#' within a condition.
#'
#' @param data A matrix of data values.
#' 
#' @param rna_seq_data An object containing the preprocessed RNA-seq data,
#' such as the output from `limma::voom` or a similar preprocessing pipeline.
#' 
#' @param meta A dataframe containing metadata, including a 'Time' column.
#' 
#' @param design A design formula or matrix for the limma analysis.
#' 
#' @param design2design_matrix_result A named list returned by
#'   `design2design_matrix()` containing at least the elements
#'   `design_matrix`, which is the fixed-effects design matrix used for
#'   limma fits, and `meta`, which is the spline-transformed metadata
#'   passed to `variancePartition::dream()` when random effects are used.
#'   
#' @param contrast_plan A list describing the unified contrast structure
#'   for the analysis. It must contain a joint contrast matrix `L_all`,
#'   a mapping object `map_for_L_all` that defines how contrast columns
#'   correspond to reported effect categories, and `spline_cols`. For
#'   variancePartition fits, this contrast plan is applied at model
#'   fitting time via the `L` argument to `dream()`. For limma fits, the
#'   same plan is used downstream to apply contrasts post hoc.
#'   
#' @param effects Strings of the fixed and potentially random effects of the 
#' design.
#' 
#' @param dream_params A named list or NULL. When not NULL, it must at least
#' contain the named element 'random_effects', which must contain a string that
#' is a formula for the random effects of the mixed models by dream.
#' Additionally, it can contain the named elements dof, which must be a int
#' bigger than 1, which is the degree of freedom for the dream topTable, and
#' the named element KenwardRoger, which must be a bool, specifying whether
#' to use that method or not.
#' 
#' @param spline_params A list of spline parameters for the analysis.
#' 
#' @param condition A character string of the column name of meta that contains
#'                  the levels of the experimental condition.
#'                  
#' @param padjust_method A character string specifying the p-adjustment method.
#' 
#' @param use_array_weights Boolean value specifying whether to use array
#' weights for limma or variancePartition.
#' 
#' @param bp_cfg A named numeric vector specifying the parallelization
#'   configuration, with expected names `"n_cores"` and `"blas_threads"`.
#'
#'   This controls how many **R worker processes** (`n_cores`) and how many
#'   **BLAS/OpenBLAS threads per process** (`blas_threads`) should be used
#'   during parallel computation.
#'
#'   If `bp_cfg` is `NULL`, missing, or any of its required fields is
#'   `NA`, both `n_cores` and `blas_threads` default to `1`. This effectively
#'   disables parallelization and avoids oversubscription of CPU threads.
#'   
#' @param verbose Boolean flag controlling the display of messages.
#'
#' @return A list containing top tables for the factor only and factor-time
#' contrast.
#'
#' @seealso
#' \code{\link[splines]{bs}}, \code{\link[splines]{ns}},
#' \code{\link[limma]{lmFit}}, \code{\link[limma]{eBayes}},
#' \code{\link[limma]{topTable}}, \code{\link{modify_limma_top_table}}
#'
#' @importFrom splines bs ns
#' @importFrom stats as.formula model.matrix
#' @importFrom limma lmFit eBayes topTable
#' @importFrom variancePartition dream
#' @importFrom BiocParallel bpstop
#'
fit_global_model <- function(
    data,
    rna_seq_data,
    meta,
    design,
    design2design_matrix_result,
    contrast_plan,
    effects,
    dream_params,
    spline_params,
    condition,
    padjust_method,
    use_array_weights,
    bp_cfg,
    verbose
    ) {
    design_matrix <- design2design_matrix_result[["design_matrix"]]

    if (!is.null(rna_seq_data)) {
        data <- rna_seq_data # Just having one variable makes the code easier
    }

    aw_result <- resolve_array_weights(
        data = data,
        rna_seq_data = rna_seq_data,
        meta = meta,
        design = design,
        design2design_matrix_result = design2design_matrix_result,
        condition = condition,
        use_array_weights = use_array_weights,
        use_random_effects = effects[["random_effects"]] != ""
    )

    if (verbose) {
        message("\nFitting global model...")
    }

    if (effects[["random_effects"]] != "") { # variancePartition approach
        colnames(data) <- rownames(meta) # dream requires this format

        # Apply the Kenward-Roger method if specified
        if (isTRUE(dream_params[["KenwardRoger"]])) {
            method <- "Kenward-Roger"
        } else {
            method <- NULL
        }

        param <- bp_setup(
            bp_cfg = bp_cfg,
            verbose = verbose
        )
        
        contrast_plan <- filter_contrast_plan_for_dream(contrast_plan)
        fit <- variancePartition::dream(
            exprObj = data,
            formula = stats::as.formula(design),
            # Spline transformed meta.
            data = design2design_matrix_result[["meta"]], 
            ddf = method,
            useWeights = aw_result[["use_weights"]],
            weightsMatrix = aw_result[["weights"]],
            L = contrast_plan[["L_all"]],
            BPPARAM = param # parallelization
        )

        fit <- variancePartition::eBayes(
            fit = fit,
            robust = aw_result[["use_weights"]]
        )
        # includes SOCK/FORK clusters on Windows
        if (inherits(param, "SnowParam")) { 
            BiocParallel::bpstop(param) # cleanly shuts down workers
        }
    } else { # limma approach
        fit <- limma::lmFit(
            object = data,
            design = design_matrix,
            weights = aw_result[["weights"]]
        )
    }

    list(
        fit = fit,
        contrast_plan = contrast_plan,
        design_matrix = design_matrix,
        meta = design2design_matrix_result[["meta"]],
        condition = condition,
        feature_names = rownames(data),
        padjust_method = padjust_method,
        homosc_violation_result = aw_result[["homosc_violation_result"]],
        spline_params = spline_params
    )
}


#' Create All Top Tables for Spline-Based Analyses
#'
#' @noRd
#'
#' @description
#' Creates all required top tables from a single global spline-based model
#' fit using a unified contrast plan. The function automatically dispatches
#' between limma and variancePartition (dream) based on the fitted model
#' object and applies the appropriate contrast and moderation strategy.
#'
#' For limma fits, contrasts are applied post hoc using
#' `limma::contrasts.fit()` followed by moderation with
#' `limma::eBayes()`. For variancePartition fits, contrasts are evaluated
#' during model fitting via `dream(..., L = ...)`, and statistics are
#' extracted directly from the fitted object.
#'
#' The function produces three categories of results:
#' \itemize{
#'   \item within-level spline time effects for each condition level
#'   \item average differences between condition levels (main and time terms)
#'   \item condition-time interaction differences (time terms only)
#' }
#'
#' For the between-level categories, the coefficient columns in the returned
#' top tables are renamed to match the historical convention:
#' `diff_main` followed by `diff_X1`, `diff_X2`, ... (or `diff_X` when
#' dof = 1) for average differences, and `diff_X*` for interaction-only
#' tests. Within-level time effect tables use `spline_c1`, `spline_c2`,
#' ... as coefficient column names.
#'
#' @param fit_obj A list containing the fitted model object in `fit`, the
#'   unified contrast plan in `contrast_plan`, and analysis metadata. For
#'   limma fits, `fit` must be the output of `limma::lmFit()` (not already
#'   moderated). For variancePartition fits, `fit` must be the output of
#'   `variancePartition::dream()` with contrasts supplied at fit time.
#' @param condition A character string of the column name of `meta` that
#'   contains the experimental condition levels. This is used to generate
#'   stable names for the within-level time effect results.
#' @param feature_names Optional character vector of feature identifiers.
#'   When `NULL`, defaults to `fit_obj[["feature_names"]]`.
#'
#' @return A named list containing three result categories:
#' \describe{
#'   \item{time_effect}{A list of within-level spline time-effect results,
#'     named by condition and level.}
#'   \item{avrg_diff_conditions}{A list of between-level average-difference
#'     results named with the prefix `avrg_diff_`.}
#'   \item{interaction_condition_time}{A list of between-level interaction
#'     results named with the prefix `time_interaction_`.}
#' }
#'
#' @seealso
#' \code{\link[limma]{lmFit}}, \code{\link[limma]{contrasts.fit}},
#' \code{\link[limma]{eBayes}}, \code{\link[limma]{topTable}},
#' \code{\link[variancePartition]{dream}},
#' \code{\link[variancePartition]{topTable}}
#'
#' @importFrom limma contrasts.fit eBayes topTable
#'
create_all_toptables <- function(
        fit_obj,
        condition,
        feature_names = NULL
) {
    if (is.null(feature_names)) {
        feature_names <- fit_obj[["feature_names"]]
    }
    
    fit <- fit_obj[["fit"]]
    if (is.null(fit)) {
        stop("fit_obj$fit is required.")
    }
    
    contrast_plan <- fit_obj[["contrast_plan"]]
    if (is.null(contrast_plan)) {
        stop("fit_obj$contrast_plan is required.")
    }
    
    map <- contrast_plan[["map_for_L_all"]]
    spline_cols <- contrast_plan[["spline_cols"]]
    
    if (is.null(map) || is.null(spline_cols)) {
        stop("contrast_plan must contain map_for_L_all and spline_cols.")
    }
    
    padj <- fit_obj[["padjust_method"]]
    if (is.null(padj)) {
        padj <- "BH"
    }
    
    is_dream_fit <- inherits(fit, "MArrayLM2") || inherits(fit, "dream")
    use_dream <- isTRUE(is_dream_fit)
    
    if (!use_dream) {
        L_all <- contrast_plan[["L_all"]]
        if (is.null(L_all)) {
            stop("contrast_plan$L_all is required for limma fits.")
        }
        
        fit_c <- limma::contrasts.fit(fit, L_all)
        
        fit_c <- limma::eBayes(
            fit_c,
            robust = isTRUE(
                fit_obj[["homosc_violation_result"]][["use_weights"]]
            )
        )
        
        top_fun <- function(coef_names) {
            limma::topTable(
                fit_c,
                coef = coef_names,
                adjust.method = padj,
                number = Inf,
                sort.by = if (length(coef_names) > 1L) "F" else "t"
            )
        }
        
        fit_for_proc <- fit_c
        intercept_fit <- fit
    } else {
        top_fun <- function(coef_names) {
            variancePartition::topTable(
                fit,
                coef = coef_names,
                adjust.method = padj,
                number = Inf,
                sort.by = if (length(coef_names) > 1L) "F" else "t"
            )
        }
        
        fit_for_proc <- fit
        intercept_fit <- fit
    }
    
    rename_spline_c <- function(tt) {
        coef_cols <- grep("^Coef[0-9]+$", colnames(tt), value = TRUE)
        if (length(coef_cols) < 1L) {
            return(tt)
        }
        
        colnames(tt)[match(coef_cols, colnames(tt))] <-
            paste0("spline_c", seq_along(coef_cols))
        tt
    }
    
    rename_diff_cols <- function(tt, has_main) {
        coef_cols <- grep("^Coef[0-9]+$", colnames(tt), value = TRUE)
        
        if (length(coef_cols) < 1L) {
            return(tt)
        }
        
        new_names <- character(length(coef_cols))
        
        if (isTRUE(has_main)) {
            new_names[1] <- "diff_main"
            if (length(coef_cols) > 1L) {
                new_names[-1] <- paste0(
                    "diff_",
                    spline_cols[seq_len(length(coef_cols) - 1L)]
                )
            }
        } else {
            new_names <- paste0(
                "diff_",
                spline_cols[seq_len(length(coef_cols))]
            )
        }
        
        colnames(tt)[match(coef_cols, colnames(tt))] <- new_names
        tt
    }
    
    lev_keys <- names(map$within_time)
    time_effect <- lapply(lev_keys, function(lev_key) {
        coef_names <- map$within_time[[lev_key]]
        top <- top_fun(coef_names)
        top <- rename_spline_c(top)
        
        process_top_table(
            list(top_table = top, fit = fit_for_proc),
            feature_names = feature_names,
            intercept_fit = intercept_fit
        )
    })
    names(time_effect) <- paste(condition, lev_keys, sep = "_")
    
    if (!is.null(map$between_main_and_time)) {
        pair_keys <- names(map$between_main_and_time)
        cond_only <- lapply(pair_keys, function(pair_key) {
            coef_names <- map$between_main_and_time[[pair_key]]
            top <- top_fun(coef_names)
            top <- rename_diff_cols(top, has_main = TRUE)
            
            process_top_table(
                list(top_table = top, fit = fit_for_proc),
                feature_names = feature_names,
                intercept_fit = intercept_fit
            )
        })
        names(cond_only) <- paste0("avrg_diff_", pair_keys)
    } else {
        pair_keys <- names(map$between_condition_only)
        cond_only <- lapply(pair_keys, function(pair_key) {
            coef_name <- map$between_condition_only[[pair_key]]
            top <- top_fun(coef_name)
            top <- rename_diff_cols(top, has_main = TRUE)
            
            process_top_table(
                list(top_table = top, fit = fit_for_proc),
                feature_names = feature_names,
                intercept_fit = intercept_fit
            )
        })
        names(cond_only) <- paste0("avrg_diff_", pair_keys)
    }
    
    pair_keys2 <- names(map$between_condition_time)
    cond_time <- lapply(pair_keys2, function(pair_key) {
        coef_names <- map$between_condition_time[[pair_key]]
        top <- top_fun(coef_names)
        top <- rename_diff_cols(top, has_main = FALSE)
        
        process_top_table(
            list(top_table = top, fit = fit_for_proc),
            feature_names = feature_names,
            intercept_fit = intercept_fit
        )
    })
    names(cond_time) <- paste0("time_interaction_", pair_keys2)
    
    list(
        time_effect = time_effect,
        avrg_diff_conditions = cond_only,
        interaction_condition_time = cond_time
    )
}


# Level 2 internal functions ---------------------------------------------------


#' Process Top Table
#'
#' @noRd
#'
#' @description
#' Processes the top table from a LIMMA analysis, adding feature names and
#' intercepts.
#'
#' @param process_within_level_result List of lists containing the limma
#'                                    topTable, and fit. All of this is from
#'                                    one specific level.
#' @param feature_names A non-empty character vector of feature names.
#'
#' @return A dataframe containing the processed top table with added intercepts.
#'
#' @seealso
#' \link{modify_limma_top_table}, \link[limma]{lmFit}
#'
#' @importFrom stats coef
#'
process_top_table <- function(
        process_within_level_result,
        feature_names,
        intercept_fit = NULL
) {
    top_table <- process_within_level_result$top_table
    fit <- process_within_level_result$fit
    
    top_table <- modify_limma_top_table(
        top_table,
        feature_names
    )
    
    fit_i <- intercept_fit
    if (is.null(fit_i)) {
        fit_i <- fit
    }
    
    coefs_i <- stats::coef(fit_i)
    if (!"(Intercept)" %in% colnames(coefs_i)) {
        stop(
            "No (Intercept) in fit used for intercept extraction."
        )
    }
    
    intercepts <- as.data.frame(coefs_i[, "(Intercept)", drop = FALSE])
    intercepts_ordered <- intercepts[top_table$feature_nr, , drop = FALSE]
    top_table$intercept <- intercepts_ordered[, 1]
    
    top_table
}


#' Process Within Level
#'
#' @noRd
#'
#' @description
#' Performs a within-level analysis using limma to generate top tables and fit
#' objects based on the specified spline parameters. Performs the limma spline
#' analysis for a selected level of a factor
#'
#' @param data A matrix of data values.
#' @param rna_seq_data An object containing the preprocessed RNA-seq data,
#' such as the output from `limma::voom` or a similar preprocessing pipeline.
#' @param meta A dataframe containing metadata, including a 'Time' column.
#' @param design A design formula or matrix for the limma analysis.
#' @param dream_params A named list or NULL. When not NULL, it it can contain
#' the named elements dof, which must be a int
#' bigger than 1, which is the degree of freedom for the dream topTable, and
#' the named element KenwardRoger, which must be a bool, specifying whether
#' to use that method or not.
#' @param spline_params A list of spline parameters for the analysis.
#' @param level_index The index of the level within the factor.
#' @param padjust_method A character string specifying the p-adjustment method.
#' @param use_array_weights Boolean value specifying whether to use array
#' weights for limma or variancePartition.
#' @param bp_cfg A named numeric vector specifying the parallelization
#'   configuration, with expected names `"n_cores"` and `"blas_threads"`.
#'
#'   This controls how many **R worker processes** (`n_cores`) and how many
#'   **BLAS/OpenBLAS threads per process** (`blas_threads`) should be used
#'   during parallel computation.
#'
#'   If `bp_cfg` is `NULL`, missing, or any of its required fields is
#'   `NA`, both `n_cores` and `blas_threads` default to `1`. This effectively
#'   disables parallelization and avoids oversubscription of CPU threads.
#' @param verbose Boolean flag controlling the display of messages.
#'
#' @return A list containing the top table and the fit object from the limma
#' analysis, and the potentially updated spline_params.
#'
#' @seealso
#' \link[splines]{bs}, \link[splines]{ns}, \link[limma]{lmFit},
#' \link[limma]{eBayes}, \link[limma]{topTable}
#'
#' @importFrom splines bs ns
#' @importFrom stats as.formula model.matrix
#' @importFrom limma lmFit eBayes topTable
#' @importFrom variancePartition dream
#'
process_within_level <- function(
    data,
    rna_seq_data,
    meta,
    design,
    dream_params,
    spline_params,
    level_index,
    padjust_method,
    use_array_weights,
    bp_cfg,
    verbose
    ) {
    effects <- extract_effects(design)

    if (spline_params[["dof"]][level_index] == 0) {
        best_dof <- select_spline_dof_loocv(
            data = data,
            meta = meta,
            spline_params = spline_params,
            level_index = level_index,
            fixed_effects = effects[["fixed_effects"]]
        )
        spline_params$dof[level_index] <- best_dof
    }

    design2design_matrix_result <- design2design_matrix(
        meta = meta,
        spline_params = spline_params,
        level_index = level_index,
        design = effects[["fixed_effects"]]
    )

    design_matrix <- design2design_matrix_result[["design_matrix"]]

    aw_result <- resolve_array_weights(
        data = data,
        rna_seq_data = rna_seq_data,
        meta = meta,
        design = design,
        design2design_matrix_result = design2design_matrix_result,
        use_array_weights = use_array_weights,
        use_random_effects = effects[["random_effects"]] != "",
        verbose = verbose
    )

    if (!is.null(rna_seq_data)) {
        data <- rna_seq_data
    }

    message(paste("\nFitting model for level ", level_index))

    if (effects[["random_effects"]] != "") {
        colnames(data) <- rownames(meta) # dream wants it like this.

        if (isTRUE(dream_params[["KenwardRoger"]])) {
            method <- "Kenward-Roger"
        } else {
            # Kenward-Roger for < 20 samples, else Satterthwaite
            method <- "adaptive" 
        }

        param <- bp_setup(
            bp_cfg = bp_cfg,
            verbose = verbose
        )

        fit <- variancePartition::dream(
            exprObj = data,
            formula = stats::as.formula(design),
            data = design2design_matrix_result[["meta"]],
            ddf = method,
            useWeights = aw_result[["use_weights"]],
            weightsMatrix = aw_result[["weights"]],
            BPPARAM = param
        )
        fit <- variancePartition::eBayes(
            fit = fit,
            robust = aw_result[["use_weights"]]
        )
        # includes SOCK/FORK clusters on Windows
        if (inherits(param, "SnowParam")) { 
            BiocParallel::bpstop(param) # cleanly shuts down workers
        }

        num_matching_columns <- sum(grepl("^X\\d+$", colnames(design_matrix)))
        coeffs <- paste0("X", seq_len(num_matching_columns))

        if (!is.null(dream_params[["dof"]])) {
            dof <- dream_params[["dof"]]
        } else {
            dof <- Inf
        }

        top_table <- variancePartition::topTable(
            fit,
            adjust.method = padjust_method,
            number = dof,
            coef = coeffs,
            sort.by = if (length(coeffs) > 1) "F" else "t"
        )
    } else {
        fit <- limma::lmFit(
            object = data,
            design = design_matrix,
            weights = aw_result[["weights"]]
        )
        fit <- limma::eBayes(
            fit = fit,
            robust = aw_result[["use_weights"]]
        )

        # Extract the top table based on coefficients
        num_matching_columns <- sum(grepl("^X\\d+$", colnames(design_matrix)))
        coeffs <- paste0("X", seq_len(num_matching_columns))

        top_table <- limma::topTable(
            fit,
            adjust.method = padjust_method,
            number = Inf,
            coef = coeffs
        )
    }

    attr(top_table, "adjust.method") <- padjust_method

    list(
        top_table = top_table,
        fit = fit,
        spline_params = spline_params,
        use_array_weights = aw_result[["use_weights"]]
    )
}


#' Remove intercept from a formula
#'
#' @noRd
#'
#' @description
#' This function modifies a given formula by replacing the first occurrence
#' of a standalone intercept (`1`) with `0`. It works even if the `1` is
#' preceded by a tilde (`~`), ensuring that the intercept is removed while
#' leaving other parts of the formula intact.
#'
#' @param formula A formula object. The formula can include an intercept (`1`)
#'   and other terms. If a `1` is found, it is replaced with `0`.
#'
#' @return A modified formula with the intercept removed. The first standalone
#'   occurrence of `1` will be replaced by `0`.
#'
remove_intercept <- function(formula) {
    # Deparse the formula into a single string
    formula_str <- paste(deparse(formula), collapse = " ")

    # Regular expression to match the first standalone 1
    # (with optional spaces around it)
    pattern <- "(~\\s*|\\s+)(1)(\\s|$)"

    # Replace the first occurrence of "1" with "0"
    formula_str <- sub(pattern, "\\10\\3", formula_str)

    new_formula <- as.formula(formula_str)
}


#' Select Optimal Spline Degrees of Freedom via Analytical LOOCV
#'
#' @noRd
#'
#' @description
#' This internal helper function selects the optimal number of degrees of
#' freedom (DoF) for a spline-based time model by minimizing the average
#' analytical leave-one-out cross-validation (LOOCV) error across all features
#' (e.g., genes).
#'
#' The method assumes that spline basis functions are included in a linear model
#' and uses \code{lm()} to fit each feature individually. The DoF that minimizes
#' the average LOOCV error across all features is returned.
#'
#' @param data A numeric matrix with features (e.g., genes) in rows and samples
#'   in columns.
#' @param meta A data frame containing sample metadata. Must include a
#'   \code{Time} column.
#' @param spline_params A list of spline parameters. Must include
#'   \code{spline_type} and \code{dof}. The entry at \code{level_index} must be
#'   set to \code{"auto"} to trigger optimization.
#' @param level_index Integer index specifying which element in
#'   \code{spline_params$dof} to update.
#' @param fixed_effects A formula or character string specifying the fixed
#'   effects for model construction.
#'
#' @return An integer giving the optimal degrees of freedom selected via LOOCV.
#'
select_spline_dof_loocv <- function(
    data,
    meta,
    spline_params,
    level_index,
    fixed_effects) {
    if (spline_params$spline_type[level_index] == "n") {
        candidate_dofs <- 2:min(10, length(unique(meta$Time)))
    } else if (spline_params$spline_type[level_index] == "b") {
        candidate_dofs <- 4:min(10, length(unique(meta$Time)))
    }

    avg_loocv_errors <- numeric(length(candidate_dofs))

    pb <- progress_bar$new(
        format = "  Optimizing DoF [:bar] :current/:total (:percent) ETA: :eta",
        total = length(candidate_dofs),
        clear = FALSE,
        width = 60
    )

    for (i in seq_along(candidate_dofs)) {
        dof_i <- candidate_dofs[i]
        spline_params_i <- spline_params
        spline_params_i$dof <- numeric(length(spline_params$dof))
        spline_params_i$dof[level_index] <- as.numeric(dof_i)

        design_result <- design2design_matrix(
            meta = meta,
            spline_params = spline_params_i,
            level_index = level_index,
            design = fixed_effects
        )
        design_matrix_i <- design_result$design_matrix

        loocv_errors <- apply(
            data,
            1,
            analytic_loocv,
            design_matrix = design_matrix_i
        )
        avg_loocv_errors[i] <- mean(
            loocv_errors,
            na.rm = TRUE
        )
        pb$tick()
    }

    best_dof <- candidate_dofs[which.min(avg_loocv_errors)]
    message(sprintf("Selected optimal spline DoF = %d", best_dof))

    return(best_dof)
}


#' Set up a BiocParallel backend in a cross-platform, RAM-safe way
#'
#' @noRd
#'
#' @description
#' Creates and registers a `BiocParallelParam` that is appropriate for the
#' current platform **and** safe on memory-constrained machines:
#'
#' * **Linux / macOS** – uses `MulticoreParam()` (fork).
#' * **Windows**       – uses `SnowParam()`  (SOCK cluster).
#' * **n_cores <= 1**  – falls back to `SerialParam()`.
#'
#' Before spawning workers it throttles the number of BLAS/OpenBLAS threads
#' **per worker** to avoid *workers × threads* oversubscription, which can slow
#' jobs down or exhaust RAM.
#'
#' @param bp_cfg A named numeric vector specifying the parallelization
#'   configuration, with expected names `"n_cores"` and `"blas_threads"`.
#'
#'   This controls how many **R worker processes** (`n_cores`) and how many
#'   **BLAS/OpenBLAS threads per process** (`blas_threads`) should be used
#'   during parallel computation.
#'
#'   If `bp_cfg` is `NULL`, missing, or any of its required fields is
#'   `NA`, both `n_cores` and `blas_threads` default to `1`. This effectively
#'   disables parallelization and avoids oversubscription of CPU threads.
#' @param verbose Boolean flag controlling the display of messages.
#'
#' @return A `BiocParallelParam` object (already registered via
#'   `BiocParallel::register()`).  Use it directly or rely on the global
#'   registration inside BiocParallel-aware functions.
#'
#' @importFrom RhpcBLASctl blas_set_num_threads
#' @importFrom BiocParallel SerialParam SnowParam MulticoreParam register
#'                          bpstart
#'
bp_setup <- function(
    bp_cfg,
    verbose) {
    # Defaults
    default <- list(n_cores = 1L, blas_threads = 1L)

    # NULL case — fallback
    if (is.null(bp_cfg)) {
        n_cores <- default$n_cores
        blas_threads <- default$blas_threads
    } else {
        # One or both may be missing — check existence before extracting
        n_cores <- if ("n_cores" %in% names(bp_cfg)) {
            as.integer(bp_cfg[["n_cores"]])
        } else {
            default$n_cores
        }

        blas_threads <- if ("blas_threads" %in% names(bp_cfg)) {
            as.integer(bp_cfg[["blas_threads"]])
        } else {
            default$blas_threads
        }
    }

    # 1.  Throttle BLAS threads in *this* session (forks inherit)
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
        RhpcBLASctl::blas_set_num_threads(blas_threads)
    }

    # 2.  Choose backend
    is_windows <- identical(.Platform$OS.type, "windows")
    if (n_cores <= 1) {
        param <- BiocParallel::SerialParam()
    } else if (is_windows) {
        param <- BiocParallel::SnowParam(
            workers = n_cores,
            type = "PSOCK",
            exportglobals = TRUE
        )
    } else {
        param <- BiocParallel::MulticoreParam(workers = n_cores)
    }

    if (verbose) {
        message(paste(
            "\nNOTE: If you manually stop run_limma_splines() in RStudio and ",
            "used parallelization for variancePartition::dream(), then those",
            "parallelized",
            "processes may continue running. Use your system's",
            "process manager to terminate them manually!\n"
        ))
    }

    # Stop any existing cluster before switching backends
    try(BiocParallel::bpstop(BiocParallel::bpparam()), silent = TRUE)

    # 3) Start, register, and auto-stop on exit
    BiocParallel::bpstart(param)
    on.exit(
        {
            try(BiocParallel::bpstop(param), silent = TRUE)
        },
        add = TRUE
    )
    BiocParallel::register(param)

    invisible(param)
}


#' Build Contrast Matrix for Spline-Based Time Models
#'
#' @noRd
#'
#' @description
#' Constructs a contrast matrix for testing spline-based time effects within a
#' linear model, optionally including interaction terms with a condition factor.
#' The contrast matrix is structured to test each spline basis coefficient
#' separately, comparing a specific condition level against the baseline.
#'
#' For the baseline level, only the spline terms are tested. For other levels,
#' both the spline terms and the corresponding interaction terms are included.
#' Interaction terms are identified dynamically from the design matrix column
#' names to ensure robustness against formula term ordering.
#'
#' @param lev Character. The condition level for which the contrast is built.
#' @param baseline Character. The baseline condition level.
#' @param condition Character. The condition factor name used in the model.
#' @param spline_terms Character vector. The names of the spline terms
#' (e.g., "X1", "X2").
#' @param design_cols Character vector. Column names of the design matrix
#' (typically from \code{colnames(fit$coefficients)}).
#' @param dof Integer. Degrees of freedom, i.e., the number of spline basis
#'  terms.
#'
#' @return A contrast matrix with rows corresponding to model coefficients and
#' columns corresponding to pseudo-contrasts for each spline term.
#' This matrix can be used with \code{limma::contrasts.fit()}.
#'
build_spline_contrast <- function(
    lev,
    baseline,
    condition,
    spline_terms,
    design_cols,
    dof
    ) {
    contrast_matrix <- matrix(
        0,
        nrow = dof,
        ncol = length(design_cols),
        dimnames = list(
            paste0("spline_c", seq_len(dof)), # pseudo coef names
            design_cols # must match colnames(fit$coefficients)
        )
    )

    if (lev == baseline) {
        # Only spline coefficients are tested
        for (k in seq_len(dof)) {
            contrast_matrix[k, spline_terms[k]] <- 1
        }
    } else {
        for (k in seq_len(dof)) {
            contrast_matrix[k, spline_terms[k]] <- 1
            # Try both possible interaction name orders:
            # ConditionLev:Xk  and  Xk:ConditionLev
            pattern1 <- paste0(
                "^",
                condition,
                lev,
                ":",
                spline_terms[k],
                "$"
            )
            pattern2 <- paste0(
                "^",
                spline_terms[k],
                ":",
                condition,
                lev,
                "$"
            )
            col_idx <- grep(pattern1, design_cols)
            if (length(col_idx) == 0L) {
                col_idx <- grep(pattern2, design_cols)
            }
            if (length(col_idx) != 1L) {
                stop_call_false(
                    "Could not uniquely identify interaction term for ",
                    spline_terms[k], " and ", condition, lev
                )
            }
            contrast_matrix[k, col_idx] <- 1
        }
    }

    # Transpose → rows = coefficients in the model,
    # columns = pseudo-contrasts (one per spline term)
    t(contrast_matrix)
}


#' Filter Contrast Plan for variancePartition (dream)
#'
#' @noRd
#'
#' @description
#' Adapts a unified contrast plan for use with variancePartition (dream) by
#' removing redundant contrasts that consist of a single non-zero coefficient.
#' Such contrasts are already evaluated by default in dream and would
#' otherwise trigger warnings when supplied via the `L` argument.
#'
#' The function removes these redundant contrast columns from the joint
#' contrast matrix and updates the corresponding contrast mappings so that
#' downstream coefficient extraction remains consistent. Dropped contrasts
#' are replaced by aliases pointing to the underlying model coefficient
#' names, ensuring that all requested tests can still be resolved by
#' `variancePartition::topTable()`.
#'
#' This transformation is intended to be applied only in the dream fitting
#' branch and allows the same high-level contrast definitions to be reused
#' across both limma and variancePartition backends.
#'
#' @param contrast_plan A list containing a joint contrast matrix `L_all` and
#'   a mapping object `map_for_L_all`, typically produced by
#'   `build_contrast_plan()`.
#'
#' @return A modified contrast plan with redundant single-term contrasts
#'   removed from `L_all`, updated mappings in `map_for_L_all`, and an
#'   additional element `alias` that records the correspondence between
#'   dropped contrast names and the underlying model coefficients.
#'
#' @seealso
#' \code{\link[variancePartition]{dream}},
#' \code{\link[variancePartition]{topTable}}
#'
#' @importFrom variancePartition dream
#'
filter_contrast_plan_for_dream <- function(contrast_plan) {
    L_all <- contrast_plan[["L_all"]]
    map <- contrast_plan[["map_for_L_all"]]
    
    if (is.null(L_all) || is.null(map)) {
        stop("contrast_plan must contain L_all and map_for_L_all.")
    }
    
    nz <- colSums(L_all != 0)
    keep <- nz != 1L
    drop <- nz == 1L
    
    alias <- list()
    
    if (any(drop)) {
        drop_names <- colnames(L_all)[drop]
        
        for (cn in drop_names) {
            v <- L_all[, cn]
            i <- which(v != 0)
            if (length(i) != 1L) {
                next
            }
            alias[[cn]] <- rownames(L_all)[i]
        }
    }
    
    L_new <- L_all[, keep, drop = FALSE]
    
    replace_alias <- function(x) {
        if (length(x) < 1L) {
            return(x)
        }
        
        unname(vapply(x, function(nm) {
            if (!is.null(alias[[nm]])) {
                alias[[nm]]
            } else {
                nm
            }
        }, character(1)))
    }
    
    map_new <- map
    map_new$within_time <- lapply(map$within_time, replace_alias)
    
    if (!is.null(map$between_condition_only)) {
        map_new$between_condition_only <-
            lapply(map$between_condition_only, replace_alias)
    }
    
    if (!is.null(map$between_condition_time)) {
        map_new$between_condition_time <-
            lapply(map$between_condition_time, replace_alias)
    }
    
    if (!is.null(map$between_main_and_time)) {
        map_new$between_main_and_time <-
            lapply(map$between_main_and_time, replace_alias)
    }
    
    contrast_plan[["L_all"]] <- L_new
    contrast_plan[["map_for_L_all"]] <- map_new
    contrast_plan[["alias"]] <- alias
    
    contrast_plan
}


# Level 3 internal functions ---------------------------------------------------


#' Modify limma Top Table
#'
#' @noRd
#'
#' @description
#' Modifies the limma top table to include feature indices and names.
#'
#' @param top_table A dataframe containing the top table results from limma
#' @param feature_names A character vector of feature names.
#'
#' @return A tibble with feature indices and names included.
#'
#' @importFrom tibble as_tibble
#' @importFrom dplyr relocate last_col mutate
#' @importFrom rlang .data
#'
modify_limma_top_table <- function(
    top_table,
    feature_names
    ) {
    is_integer_string <- function(x) {
        return(grepl("^[0-9]+$", x))
    }

    # Because the row headers of a potential rna_seq_data object were not
    # converted to ints (written as strings) beforehand. This is run only when
    # the row headers are still "real" strings.
    if (!all(vapply(
        rownames(top_table),
        is_integer_string,
        logical(1)
    ))
    ) {
        rownames(top_table) <- vapply(
            rownames(top_table),
            function(id) {
                # Find the index of the current row name in feature_names
                index <- which(feature_names == id)
                # Return the index as a string
                return(as.character(index))
            },
            character(1)
        )
    }

    top_table <- tibble::as_tibble(
        top_table,
        rownames = "feature_nr"
    )

    # Convert feature_nr to integer
    top_table <- top_table |>
        dplyr::mutate(
            feature_nr = as.integer(.data$feature_nr)
        ) |>
        dplyr::relocate(
            "feature_nr",
            .after = dplyr::last_col()
        )

    # Sort and add feature names based on the feature_nr
    sorted_feature_names <- feature_names[top_table$feature_nr]
    top_table <- top_table |> dplyr::mutate(
        feature_names = sorted_feature_names
        )

    return(top_table)
}


#' Construct contrast matrix for condition differences over time
#'
#' @noRd
#'
#' @description
#' Builds a contrast matrix \code{L} for comparing two condition levels in a
#' spline-based time model. The first contrast (\code{diff_main}) captures the
#' average condition difference, while the remaining contrasts (\code{diff_X*})
#' capture differences in temporal shape via condition-by-spline interactions.
#' Assumes a design matrix in which at most one level is represented explicitly
#' as main and interaction terms, with the other treated as the reference.
#'
#' @param condition A string giving the name of the condition factor.
#'
#' @param level_pair A character vector of length 2 specifying the two levels
#'   to be compared.
#'
#' @param design_matrix The design matrix used for model fitting, whose column
#'   names determine which condition level is modeled explicitly and which is
#'   treated as the reference.
#'
#' @return A named list with elements:
#'   - \code{L}: contrast matrix with rows aligned to
#'     \code{colnames(design_matrix)} and columns \code{diff_main} and
#'     \code{diff_<spline>} for each spline basis term
#'   - \code{main_and_time_names}: column names of \code{L}, including the main
#'     difference and all shape terms
#'   - \code{shape_names}: column names of \code{L} excluding
#'     \code{diff_main}
#'   - \code{spline_cols}: detected spline basis column names (e.g. \code{X},
#'     \code{X1}, \code{X2}, ...)
#'   - \code{main_cols}: named character vector mapping \code{level1} and
#'     \code{level2} to their expected main-effect column names in the design
#'     matrix
#'
get_condition_contrast_coefs <- function(
        condition,
        level_pair,
        design_matrix
) {
    level1 <- level_pair[1]
    level2 <- level_pair[2]
    
    cols <- colnames(design_matrix)
    if (is.null(cols) || anyNA(cols) || any(cols == "")) {
        stop("design_matrix must have non-empty column names.")
    }
    
    spline_cols <- grep("^X$|^X\\d+$", cols, value = TRUE)
    if (length(spline_cols) < 1L) {
        stop("No spline basis columns found (expected X or X1, X2, ...).")
    }
    
    main1 <- paste0(condition, level1)
    main2 <- paste0(condition, level2)
    
    idx_main1 <- match(main1, cols)
    idx_main2 <- match(main2, cols)
    
    blk1 <- match_int_either(
        main = main1,
        spline_cols = spline_cols,
        cols = cols
    )
    blk2 <- match_int_either(
        main = main2,
        spline_cols = spline_cols,
        cols = cols
    )
    
    idx_int1 <- blk1$idx
    idx_int2 <- blk2$idx
    
    if (all(is.na(c(idx_main1, idx_main2, idx_int1, idx_int2)))) {
        stop(
            "Neither level appears in design_matrix columns. ",
            "Check 'condition' prefix and factor levels."
        )
    }
    
    check_block_complete(
        level = level1,
        idx_main = idx_main1,
        idx_int = idx_int1,
        main = main1,
        spline_cols = spline_cols
    )
    check_block_complete(
        level = level2,
        idx_main = idx_main2,
        idx_int = idx_int2,
        main = main2,
        spline_cols = spline_cols
    )
    
    p <- ncol(design_matrix)
    k <- length(spline_cols)
    
    L <- matrix(
        0,
        nrow = p,
        ncol = 1L + k
    )
    rownames(L) <- cols
    colnames(L) <- c(
        "diff_main",
        paste0("diff_", spline_cols)
    )
    
    if (!is.na(idx_main2)) {
        L[idx_main2, "diff_main"] <- 1
    }
    if (!is.na(idx_main1)) {
        L[idx_main1, "diff_main"] <- -1
    }
    
    for (j in seq_len(k)) {
        cn <- paste0("diff_", spline_cols[j])
        
        if (!is.na(idx_int2[j])) {
            L[idx_int2[j], cn] <- 1
        }
        if (!is.na(idx_int1[j])) {
            L[idx_int1[j], cn] <- -1
        }
    }
    
    list(
        L = L,
        main_and_time_names = colnames(L),
        shape_names = colnames(L)[-1L],
        spline_cols = spline_cols,
        main_cols = c(
            level1 = main1,
            level2 = main2
        )
    )
}


#' Analytical Leave-One-Out Cross-Validation computation
#'
#' @noRd
#'
#' @description
#' Computes the analytical Leave-One-Out Cross-Validation (LOOCV) error
#' for a given expression vector and design matrix using a linear model.
#' This function uses the hat matrix and residuals to avoid repeated model
#' fitting and provides a fast estimate of predictive error.
#'
#' @param expr_vec A numeric vector of expression values for a single feature.
#' @param design_matrix A numeric design matrix (no intercept column) used to
#'   model the time trend or other covariates.
#'
#' @return A numeric scalar representing the LOOCV error.
#'
#' @importFrom stats lm hatvalues residuals
#'
analytic_loocv <- function(
    expr_vec,
    design_matrix) {
    fit <- stats::lm(expr_vec ~ design_matrix - 1) # no intercept, already in
    h <- stats::hatvalues(fit)
    r <- stats::residuals(fit)
    mean((r / (1 - h))^2)
}


#' Resolve whether array weights should be used, and return them if needed
#'
#' @noRd
#'
#' @description
#' This function encapsulates the logic for deciding whether to use array
#' weights
#' in a linear modeling pipeline based on the presence of heteroscedasticity in
#' the data. It supports both automatic detection (based on Levene's test) and
#' manual override via the `use_array_weights` argument. If `rna_seq_data` is
#' provided, the use of weights is suppressed since RNA-seq workflows are
#' expected
#' to have applied the voom transformation already. If heteroscedasticity is
#' detected or explicitly requested, the function computes array weights using
#' limma's `arrayWeights()` function. If random effects are present, the
#' weights are reshaped into a matrix suitable for use with
#' `variancePartition::dream()`.
#' The function returns a structured list containing the final weight matrix
#'  (or NULL),
#' the resolved boolean flag for weight usage, and the heteroscedasticity test
#'  result (if available) for diagnostic purposes.
#'
#' @param data Expression matrix (features x samples)
#' @param rna_seq_data Optional voom-transformed data. If provided, weights are
#'  skipped.
#' @param meta Metadata data frame
#' @param design Full design formula (string)
#' @param design2design_matrix_result Result from design2design_matrix()
#' @param condition Condition column name used for group testing
#' @param use_array_weights User-specified flag (TRUE/FALSE/NULL)
#' @param use_random_effects Boolean, whether model includes random effects
#' @param data_type Type of omics data (used for logging only)
#' @param p_threshold p-value threshold for Levene's test
#' @param fraction_threshold Fraction of features needing to fail to trigger
#'  use of weights
#' @param verbose Boolean flag controlling the display of messages.
#'
#' @return A list with `weights`, `use_weights`, and `homosc_violation_result`
#'
resolve_array_weights <- function(
    data,
    rna_seq_data,
    meta,
    design,
    design2design_matrix_result,
    condition = NULL,
    use_array_weights = TRUE,
    use_random_effects = FALSE,   
    verbose = FALSE
    ) {
    # Auto-detect violation if user has not explicitly set the flag
    if (is.null(use_array_weights)) {
        homosc_violation_result <- check_homoscedasticity_violation(
            data = data,
            meta = meta,
            design = design,
            design2design_matrix_result = design2design_matrix_result,
            condition = condition,
            use_random_effects = use_random_effects
        )
        use_array_weights <- homosc_violation_result$violation
    } else {
        homosc_violation_result <- NULL
    }

    if (isTRUE(use_array_weights)) {
        if (verbose) {
            message(
                "Using arrayWeights strategy for heteroscedasticity adjustment."
            )
        }
        weights <- limma::arrayWeights(
            object = data,
            design = design2design_matrix_result$design_matrix
        )
        pretty_weights <- data.frame(
            Sample = colnames(data),
            Weight = round(weights, 4)
        )
        header <- sprintf("%-25s %10s", "Sample", "Weight")
        rows <- sprintf(
            "%-25s %10.4f",
            pretty_weights$Sample,
            pretty_weights$Weight
            )
        weights_msg <- paste(
            c("Array weights:", header, rows),
            collapse = "\n"
        )
        if (verbose) {
            message(weights_msg)
        }
        
        if (use_random_effects) { # no random effects present
            weights <- matrix(
                weights,
                nrow = nrow(data),
                ncol = ncol(data),
                byrow = TRUE
            )
        }
    } else {                 # use_array_weights is FALSE
        weights <- NULL
    }

    list(
        weights = weights,
        use_weights = use_array_weights,
        homosc_violation_result = homosc_violation_result
    )
}


# Level 4 internal functions ---------------------------------------------------


#' Check completeness of an interaction block
#'
#' @noRd
#'
#' @description
#' Validate that, for a given level, all required interaction terms between a
#' main effect and its spline terms are present whenever any part of the block
#' is used. If some interactions are missing, an error is raised.
#'
#' @param level Character scalar identifying the level being checked.
#' @param idx_main Integer index of the main-effect column, or \code{NA} if
#'   absent.
#' @param idx_int Integer vector of indices for interaction terms corresponding
#'   to \code{spline_cols}, with \code{NA} indicating missing interactions.
#' @param main Character scalar giving the name of the main effect.
#' @param spline_cols Character vector of spline term names.
#'
#' @return Invisibly returns \code{NULL}. Called for its side effect of throwing
#'   an error when the interaction block is incomplete.
#'   
check_block_complete <- function(
        level,
        idx_main,
        idx_int,
        main,
        spline_cols
) {
    has_any <- !is.na(idx_main) || any(!is.na(idx_int))
    
    if (has_any && any(is.na(idx_int))) {
        missing_spl <- spline_cols[is.na(idx_int)]
        missing_terms <- c(
            paste0(missing_spl, ":", main),
            paste0(main, ":", missing_spl)
        )
        
        stop(
            "Interaction block incomplete for level '", level,
            "' (", main, "). Missing (either order): ",
            paste(missing_terms, collapse = ", ")
        )
    }
}


#' Match interaction columns irrespective of order
#'
#' @noRd
#'
#' @description
#' Identify interaction terms between a main effect and spline terms in a
#' vector of column names, allowing for either ordering of the interaction
#' (\code{spline:main} or \code{main:spline}). Ambiguous cases where both
#' orderings are present are rejected.
#'
#' @param main Character scalar giving the name of the main effect.
#' @param spline_cols Character vector of spline term names.
#' @param cols Character vector of column names in which to search.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{idx}{Integer vector of indices into \code{cols} for the matched
#'       interaction terms, or \code{NA} if no match is found.}
#'     \item{terms_a}{Character vector of interaction terms of the form
#'       \code{spline:main}.}
#'     \item{terms_b}{Character vector of interaction terms of the form
#'       \code{main:spline}.}
#'   }
#'   
match_int_either <- function(
        main,
        spline_cols,
        cols
) {
    a <- paste0(spline_cols, ":", main)
    b <- paste0(main, ":", spline_cols)
    
    idx_a <- match(a, cols)
    idx_b <- match(b, cols)
    
    idx <- ifelse(!is.na(idx_a), idx_a, idx_b)
    
    both <- !is.na(idx_a) & !is.na(idx_b)
    if (any(both)) {
        stop(
            "Ambiguous interaction columns for '", main,
            "': both orders present for ",
            paste(spline_cols[both], collapse = ", "),
            ". Remove one convention or disambiguate before fitting."
        )
    }
    
    list(
        idx = idx,
        terms_a = a,
        terms_b = b
    )
}