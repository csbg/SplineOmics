test_that("create_limma_report runs without error, writes to tempdir, and produces stable output", {
    data <- readRDS(xzfile(system.file(
        "extdata", "proteomics_data.rds.xz",
        package = "SplineOmics"
    )))
    meta <- read.csv(system.file(
        "extdata", "proteomics_meta.csv",
        package = "SplineOmics"
    ), stringsAsFactors = FALSE)

    summarise_plots <- function(x) {
        # x is whatever create_limma_report() returns (you call it "plots")
        # If your function returns a list that includes non-plot strings (headers),
        # we keep them too so the structure is visible in the snapshot.
        lapply(x, function(el) {
            if (inherits(el, "ggplot")) {
                list(
                    type   = "ggplot",
                    layers = length(el$layers),
                    title  = el$labels$title %||% NA_character_
                )
            } else if (inherits(el, "recordedplot")) {
                list(
                    type = "recordedplot"
                    # recordedplot has no easy stable textual representation
                )
            } else if (is.list(el) && !is.null(el$plots)) {
                # If you ever pass through nested sections that have $plots
                list(
                    type   = "section",
                    nplots = length(el$plots)
                )
            } else if (is.character(el) && length(el) == 1) {
                list(
                    type  = "header",
                    text  = el
                )
            } else {
                list(type = paste(class(el), collapse = "/"))
            }
        })
    }

    extracted_data <- SplineOmics::extract_data(
        data = data,
        feature_name_columns = c("Gene_name"),
        top_row = 4,
        bottom_row = 1165,
        right_col = 37,
        left_col = 2
    )

    first_na_col <- which(is.na(data[1, ]))[1]
    annotation <- data |>
        dplyr::select((first_na_col + 1):ncol(data)) |>
        dplyr::slice(-c(1:3))

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

    # --- First run ---
    splineomics <- SplineOmics::update_splineomics(
        splineomics = splineomics,
        design = "~ 1 + Time + Reactor",
        mode = "isolated",
        spline_params = list(
            spline_type = c("n", "n"),
            dof = c(2L, 2L)
        )
    )
    splineomics <- SplineOmics::run_limma_splines(splineomics)

    temp_report_dir <- tempdir()

    # Run report and capture output
    expect_error(
        plots_first <- SplineOmics::create_limma_report(
            splineomics = splineomics,
            report_dir = temp_report_dir
        ),
        regexp = NA
    )
    expect_snapshot(
        summarise_plots(plots_first)
    )

    splineomics <- SplineOmics::create_splineomics(
        data = extracted_data,
        meta = meta,
        annotation = annotation,
        report_info = report_info,
        condition = "Phase",
        meta_batch_column = "Reactor"
    )

    # --- Second run with changed params ---
    splineomics <- SplineOmics::update_splineomics(
        splineomics = splineomics,
        design = "~ 1 + Phase*Time + Reactor",
        mode = "integrated",
        spline_params = list(
            spline_type = c("n"),
            dof = c(2L)
        )
    )
    splineomics <- SplineOmics::run_limma_splines(splineomics)

    expect_error(
        plots_second <- SplineOmics::create_limma_report(
            splineomics = splineomics,
            report_dir = temp_report_dir
        ),
        regexp = NA
    )
    expect_snapshot(
        summarise_plots(plots_second)
    )
})
