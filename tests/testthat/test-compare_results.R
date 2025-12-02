test_that("compare_results snapshot with run_limma_splines", {
    skip_on_cran()
    set.seed(123)
    
    # --- Toy data -------------------------------------------------------------
    toy_data <- matrix(
        c(
            3, 5, 8, 12, 17, 23,   # f1
            23, 17, 13, 9, 6, 4,   # f2
            5, 3, 2, 2, 3, 5,      # f3
            1, 4, 9, 8, 4, 1,      # f4
            10, 10, 10, 10, 10, 10, # f5
            2, 2, 2, 9, 12, 15,    # f6
            4, 5, 7, 10, 14, 19,   # f7
            12, 11, 9, 8, 9, 12    # f8
        ),
        nrow = 8, ncol = 6, byrow = TRUE,
        dimnames = list(paste0("f", 1:8), paste0("s", 1:6))
    )
    
    toy_meta <- data.frame(
        Time      = c(0, 1, 2, 0, 1, 2),
        condition = rep(c("WT", "KO"), each = 3),
        Replicate = rep(c("R1", "R2"), each = 3),
        row.names = colnames(toy_data),
        stringsAsFactors = FALSE
    )
    
    toy_annot <- data.frame(
        feature_nr = 1:8,
        gene       = paste0("G", 1:8),
        stringsAsFactors = FALSE
    )
    
    design_str <- "~ 1 + Time*condition"
    
    spline_params <- list(
        spline_type = "n",
        dof         = 1L
    )
    
    # Stub limma "top tables" (will typically be overwritten by 
    # run_limma_splines)
    tt_stub <- data.frame(
        feature_nr = 1:4,
        adj.P.Val  = c(0.05, 0.5, 0.5, 0.5)
    )
    
    toy_splineomics <- list(
        data               = toy_data,
        meta               = toy_meta,
        annotation         = toy_annot,
        report_info        = list(
            omics_data_type      = "RNA-seq",
            data_description     = "toy example",
            data_collection_date = "2025-01-01",
            analyst_name         = "Example",
            contact_info         = "example@example.org",
            project_name         = "ToyProject1"
        ),
        design             = design_str,
        mode               = "integrated",
        condition          = "condition",
        spline_params      = spline_params,
        meta_batch_column  = NULL,
        meta_batch2_column = NULL,
        limma_splines_result = list(
            time_effect                = list(WT = tt_stub, KO = tt_stub),
            # NOTE: using list() here to match compare_results' expectations
            avrg_diff_conditions       = list(KO_vs_WT = tt_stub),
            interaction_condition_time = list(KO_vs_WT = tt_stub)
        ),
        feature_name_columns = "gene"
    )
    class(toy_splineomics) <- "SplineOmics"
    
    # --- First run ------------------------------------------------------------
    splineomics1 <- run_limma_splines(toy_splineomics)
    
    # --- Second run with slight perturbation ----------------------------------
    toy_data2 <- toy_data
    toy_data2["f1", ] <- toy_data2["f1", ] + 1  # change p-values a bit
    
    toy_splineomics2 <- toy_splineomics
    toy_splineomics2$data <- toy_data2
    toy_splineomics2$report_info$project_name <- "ToyProject2"
    
    splineomics2 <- run_limma_splines(toy_splineomics2)
    
    # --- Compare and snapshot -------------------------------------------------
    res <- compare_results(
        splineomics1,
        splineomics2,
        splineomics1_description = "run1",
        splineomics2_description = "run2",
        adj_p_tresh1 = 0.05,
        adj_p_tresh2 = 0.05
    )
    
    # Reduce to a stable, text-friendly object for snapshotting
    snapshot_obj <- list(
        correlation_summary = res$correlation_summary,
        hits_summary        = res$hits_summary,
        plot_names          = names(res$plots)
    )
    
    expect_snapshot_value(snapshot_obj, style = "json2")
})
