#' Package Name: SplineOmics
#'
#' @description
#' The R package SplineOmics finds the significant features (hits) of
#' time-series -omics data by using splines and limma for hypothesis testing.
#' It then clusters the hits based on the spline shape while showing all
#' results in summary HTML reports.
#'
#' For detailed documentation, vignettes, and examples, please visit the
#' [SplineOmics GitHub page](https://github.com/csbg/SplineOmics.git).
#'
#' @section Key Functions and Classes:
#' - extract_data: Extracts data matrix from Excel file.
#' - create_splineomics: Creates the SplineOmics object, which contains arguments
#'                       used by several package functions.
#' - explore_data: Performs exploratory data analysis with the data, and outputs
#'                 an HTML report containg various plots, such as density plots
#'                 and correlation heatmaps.
#' - screen_limma_hyperparams: Allows the specify lists of different hyperparameters
#'                             to test, such as a degree of freedom of 2, 3, 4,
#'                             and adj.p-val thresholds, such as 0.1 and 0.05,
#'                             and tests all specified different values for all
#'                             limma spline hyperparameters in a semi-combinatorial
#'                             way.
#' - update_splineomics: Allows to change values of the SplineOmics object, for
#'                       example after observing that outliers should be removed
#'                       from the data (update the data parameter).
#' - run_limma_splines: Central function of the script, is called by the
#'                      screen_limma_hyperparams function and can be called to
#'                      get the limma spline analysis results (p-values for all
#'                      features (e.g. proteins)) with the hyperparameters, that
#'                       were selected finally.
#' - create_limma_report: Creates an HTML report showing the run_limma_splines
#'                        results
#' - cluster_hits: Clusters the splines of the hits (significant features) based
#'                  on their shape and shows all results as plots in an HTML
#'                  report.
#' - download_enrichr_databases: Allows to download the Enrichr databases for
#'                               runnin clusterProfiler in the run_gsea function
#'                               with them.
#' - run_gsea: Runs clusterProfiler with the clustered hits by using the Enrichr
#'             databases.
#'
#' @section Package Options:
#' None
#'
#' @section Dependencies:
#' -   **ComplexHeatmap**: For creating complex heatmaps with advanced features.
#' -   **base64enc**: For encoding/decoding base64.
#' -   **dendextend**: For extending `dendrogram` objects in R, allowing for easier manipulation of dendrograms.
#' -   **dplyr**: For data manipulation.
#' -   **ggplot2**: For creating elegant data visualizations using the grammar of graphics.
#' -   **ggrepel**: For better label placement in ggplot2.
#' -   **here**: For constructing paths to your project’s files.
#' -   **limma**: For linear models for microarray data.
#' -   **openxlsx**: For reading, writing, and editing xlsx files.
#' -   **patchwork**: For combining multiple ggplot objects into a single plot.
#' -   **pheatmap**: For creating pretty heatmaps.
#' -   **progress**: For adding progress bars to your loops and apply functions.
#' -   **purrr**: For functional programming tools.
#' -   **rlang**: For tools to work with core language features of R and R’s base types.
#' -   **scales**: For scale functions for visualization.
#' -   **tibble**: For creating tidy data frames that are easy to work with.
#' -   **tidyr**: For tidying your data.
#' -   **zip**: For combining files into a zip file.
#'
#' Optional dependencies
#'
#' These dependencies are only necessary for some functions:
#'
#' -   **edgeR**: For preprocessing RNA-seq data in the run_limma_splines() fun.
#' -   **clusterProfiler**: For the run_gsea() function (gene set enrichment).
#' -   **rstudioapi**: For the open_tutorial() and open_template() functions.
#'
#' @section Authors:
#' - [Thomas-Rauter](https://github.com/Thomas-Rauter) - Wrote the package and
#'   developed the approach with VSchaepertoens under guidance from nfortelny
#'   and skafdasschaf.
#' - [nfortelny](https://github.com/nfortelny) - Principal Investigator,
#'   provided guidance and support.
#' - [skafdasschaf](https://github.com/skafdasschaf) - Helped review code and
#'   provided improvement suggestions.
#' - [VSchaepertoens](https://github.com/VSchaepertoens) - Developed an internal
#'   plotting function and contributed to exploratory data analysis and the
#'   overall approach.
#'
#' @section Maintainer:
#' - Name: Thomas Rauter
#' - Email: thomas.rauter@plus.ac.at
#'
#' @section License:
#' - License: MIT
#'
#' @section Useful URLs:
#' - [GitHub repo of the package](https://github.com/csbg/SplineOmics.git)
#'
#' @section Additional Information:
#' None
#'
#' @keywords omics, time-series, splines, limma, clustering, GSEA, HTML reports
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
