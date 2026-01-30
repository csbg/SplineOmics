testthat::test_that("find_pvc phosphoproteomics pipeline snapshot", {
    data <- readRDS(
        xzfile(
            system.file(
                "extdata",
                "phosphoproteomics_data.rds.xz",
                package = "SplineOmics"
            )
        )
    )
    
    meta <- read.csv(
        system.file(
            "extdata",
            "phosphoproteomics_meta.csv",
            package = "SplineOmics"
        ),
        stringsAsFactors = FALSE
    )
    
    first_na_col <- which(is.na(data[1, ]))[1]
    annotation <- data |>
        dplyr::select((first_na_col + 1):ncol(data))
    
    data <- SplineOmics::extract_data(
        data = data,
        feature_name_columns = c("T..Gene"),
        use_row_index = TRUE,
        top_row = 1,
        bottom_row = 1000,
        right_col = 36,
        left_col = 1
    )
    
    splineomics <- SplineOmics::create_splineomics(
        data = data,
        meta = meta,
        annotation = annotation,
        condition = "Phase",
        meta_batch_column = "Reactor"
    )
    
    warning_message <- testthat::capture_warnings({
        pvc_results <- SplineOmics::find_pvc(
            splineomics = splineomics,
            alphas = 0.025,
            padjust_method = "BH",
            support = 1
        )
    })
    
    testthat::expect_type(
        pvc_results,
        "list"
    )
    testthat::expect_setequal(
        names(pvc_results),
        c("Exponential", "Stationary")
    )
    
    exp <- pvc_results$Exponential
    stat <- pvc_results$Stationary
    
    testthat::expect_type(
        exp,
        "list"
    )
    testthat::expect_type(
        stat,
        "list"
    )
    
    round6 <- function(x) {
        if (is.numeric(x)) {
            return(round(x, 6))
        }
        if (is.matrix(x)) {
            return(apply(x, 2, function(col) round(col, 6)))
        }
        x
    }
    
    snapshot_payload <- list(
        warnings = warning_message,
        structure = list(
            exp_names = names(exp),
            stat_names = names(stat)
        ),
        exponential = list(
            alpha = exp$alpha,
            pvc_adj_pvals_dim = dim(exp$pvc_adj_pvals),
            pvc_adj_pvals_head5 = round6(
                exp$pvc_adj_pvals[
                    seq_len(min(5, nrow(exp$pvc_adj_pvals))),
                    ,
                    drop = FALSE
                ]
            ),
            pvc_pattern_summary = exp$pvc_pattern_summary
        ),
        stationary = list(
            alpha = stat$alpha,
            pvc_adj_pvals_dim = dim(stat$pvc_adj_pvals),
            pvc_adj_pvals_head5 = round6(
                stat$pvc_adj_pvals[
                    seq_len(min(5, nrow(stat$pvc_adj_pvals))),
                    ,
                    drop = FALSE
                ]
            ),
            pvc_pattern_summary = stat$pvc_pattern_summary
        ),
        attributes = list(
            padjust_method = attr(pvc_results, "padjust_method", exact = TRUE),
            support = attr(pvc_results, "support", exact = TRUE),
            batch_effects = attr(pvc_results, "batch_effects", exact = TRUE)
        )
    )
    
    testthat::expect_snapshot_value(
        snapshot_payload,
        style = "json2"
    )
    
    testthat::expect_silent(
        SplineOmics::create_pvc_report(
            splineomics = splineomics,
            pvc_results = pvc_results,
            verbose = FALSE,
            report_dir = withr::local_tempdir()
        )
    )
})
