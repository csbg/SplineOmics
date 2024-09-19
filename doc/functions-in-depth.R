## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----extract_data() usage, eval = FALSE---------------------------------------
#  clean_data <- extract_data(data,
#                             feature_name_columns = NA)

## ----explore_data() usage, eval = FALSE---------------------------------------
#  plots <- explore_data(data,
#                        meta,
#                        condition,
#                        report_info,
#                        meta_batch_column = NA,
#                        meta_batch2_column = NA,
#                        report_dir = here::here(),
#                        report = TRUE)

## ----screen_limma_hyperparams() usage, eval = FALSE---------------------------
#  screen_limma_hyperparams(datas,
#                           datas_descr,
#                           metas,
#                           designs,
#                           condition,
#                           spline_test_configs,
#                           report_info,
#                           report_dir = here::here(),
#                           adj_pthresholds = c(0.05),
#                           meta_batch_column = NA,
#                           meta_batch2_column = NA,
#                           time_unit = "m",
#                           padjust_method = "BH")

## ----run_limma_splines() usage, eval = FALSE----------------------------------
#  top_tables <- run_limma_splines(data,
#                                  meta,
#                                  design,
#                                  condition,
#                                  spline_params =
#                                    list(spline_type = c("n"),
#                                         dof = c(2L)),
#                                  padjust_method = "BH")

## ----create_limma_report() usage, eval = FALSE--------------------------------
#  create_limma_report(run_limma_splines_result,
#                      report_info,
#                      adj_pthresh = 0.05,
#                      report_dir = here::here())

## ----cluster_hits() usage, eval = FALSE---------------------------------------
#  clustering_results <- cluster_hits(top_tables,
#                                     data,
#                                     meta,
#                                     design,
#                                     condition,
#                                     report_info,
#                                     spline_params =
#                                      list(spline_type = c("n"),
#                                           dof = c(2L)),
#                                     adj_pthresholds = c(0.05),
#                                     clusters = c("auto"),
#                                     meta_batch_column = NA,
#                                     meta_batch2_column = NA,
#                                     time_unit = "min",
#                                     report_dir = here::here(),
#                                     report = TRUE)

## ----download_enrichr_databases() usage, eval = FALSE-------------------------
#  download_enrichr_databases(gene_set_lib,
#                             output_dir = here::here())

## ----run_gsea() usage, eval = FALSE-------------------------------------------
#  plots <- run_gsea(levels_clustered_hits,
#                    genes,
#                    databases,
#                    report_info,
#                    params = NA,
#                    plot_titles = NA,
#                    background = NULL,
#                    report_dir = here::here())

