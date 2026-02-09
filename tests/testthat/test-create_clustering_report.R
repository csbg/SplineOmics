testthat::test_that(
    "create_clustering_report runs on cluster_hits payload and writes HTML",
    {
        testthat::skip_on_cran()
        
        ## Setup (mirrors cluster_hits pipeline test)
        set.seed(42)
        
        data <- readRDS(
            xzfile(
                system.file(
                    "extdata",
                    "proteomics_data.rds.xz",
                    package = "SplineOmics"
                )
            )
        )
        
        meta <- read.csv(
            system.file(
                "extdata",
                "proteomics_meta.csv",
                package = "SplineOmics"
            ),
            stringsAsFactors = FALSE
        )
        
        first_na_col <- which(is.na(data[1, ]))[1]
        
        annotation <- data |>
            dplyr::select((first_na_col + 1):ncol(data)) |>
            dplyr::slice(-c(1:3))
        
        extracted_data <- SplineOmics::extract_data(
            data = data,
            feature_name_columns = c("Gene_name"),
            top_row = 4,
            bottom_row = 1165,
            right_col = 37,
            left_col = 2
        )
        
        report_info <- list(
            omics_data_type = "PTX",
            data_description = "Test data",
            data_collection_date = "2024",
            analyst_name = "Thomas Rauter",
            contact_info = "thomas.rauter@plus.ac.at",
            project_name = "DGTX"
        )
        
        splineomics <- SplineOmics::create_splineomics(
            data = extracted_data,
            meta = meta,
            annotation = annotation,
            report_info = report_info,
            condition = "Phase",
            meta_batch_column = "Reactor"
        )
        
        splineomics <- SplineOmics::update_splineomics(
            splineomics = splineomics,
            design = "~ 1 + Time + Reactor",
            mode = "isolated",
            use_array_weights = FALSE,
            spline_params = list(
                spline_type = c("n", "n"),
                dof = c(2L, 2L)
            )
        )
        
        splineomics <- SplineOmics::run_limma_splines(splineomics)
        
        nr_clusters <- list(
            Exponential = 6,
            Stationary = 2:3
        )
        
        gene_column_name <- "Gene_symbol"
        genes <- annotation[[gene_column_name]]
        
        clustering_results <- SplineOmics::cluster_hits(
            splineomics = splineomics,
            nr_clusters = nr_clusters,
            genes = genes
        )
        
        testthat::expect_named(
            clustering_results,
            c("cluster_table", "spline_results", "report_payload"),
            ignore.order = FALSE
        )
        
        report_payload <- clustering_results[["report_payload"]]
        
        report_dir <- withr::local_tempdir()
        
        testthat::expect_silent(
            {
                plots <- SplineOmics::create_clustering_report(
                    report_payload = report_payload,
                    report_dir = report_dir,
                    verbose = FALSE,
                    max_hit_number = 10
                )
            }
        )
        
        testthat::expect_type(plots, "list")
        
        html_files <- list.files(
            report_dir,
            pattern = "\\.html$",
            recursive = TRUE,
            full.names = TRUE
        )
        
        testthat::expect_true(length(html_files) >= 1)
    }
)

testthat::test_that(
    "create_clustering_report rejects invalid max_hit_number",
    {
        testthat::skip_on_cran()
        
        splineomics <- list(
            data = matrix(0, nrow = 2, ncol = 2),
            meta = data.frame(
                condition = c("A", "B"),
                time = c(0, 1),
                batch = c("b1", "b1"),
                stringsAsFactors = FALSE
            ),
            condition = "condition",
            mode = "timecourse",
            spline_params = list(),
            report_info = list(),
            design = ~ condition,
            meta_batch_column = "batch",
            meta_batch2_column = NULL,
            feature_name_columns = NULL,
            annotation = NULL,
            homosc_violation_result = list(
                violation = FALSE,
                percent_violated = 0
            )
        )
        
        report_payload <- list(
            splineomics = splineomics,
            design = splineomics$design,
            all_levels_clustering = list(),
            genes = c("g1", "g2"),
            predicted_timecurves = array(0, dim = c(2, 2, 2)),
            category_2_and_3_hits = list(),
            adj_pthresh_time_effect = 0.05,
            adj_pthresh_avrg_diff_conditions = 0.05,
            adj_pthresh_interaction_condition_time = 0.05,
            min_effect_size = 0
        )
        
        testthat::expect_error(
            SplineOmics::create_clustering_report(
                report_payload = report_payload,
                report_dir = withr::local_tempdir(),
                max_hit_number = -1
            )
        )
    }
)
