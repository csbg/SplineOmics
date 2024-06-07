
# SplineOmics <img src="man/figures/hex_logo.png" style="float: right; width: 150px; margin-left: 300px; vertical-align: middle;"/>

The R package SplineOmics gets the significant features (hits) of
time-series omics data by using splines and limma for hypothesis testing
and clusters the hits based on the spline shape, showing all results in
HTML reports.

## Table of Contents

- [Introduction](#introduction)
- [Installation](#installation)
- [Usage](#usage)
  - [extract_data](#extract_data)
  - [explore_data](#explore_data)
  - [run_limma_splines](#run_limma_splines)
  - [cluster_hits](#cluster_hits)
  - [run_gsea](#run_gsea)
- [Docker Container](#docker-container)
- [Dependencies](#dependencies)
- [Getting Help](#getting-help)
- [Contributing](#contributing)
- [License](#license)
- [Citation](#citation)

## Introduction

Welcome to `SplineOmics`, an R package designed to streamline the
analysis of omics time-series data, followed by automated HTML report
generation.

### Is the SplineOmics package of use for me?

If you have omics data over time, the package will help you to run limma
with splines, decide on which parameters to use and how to cluster and
plot results.

You need to have a data matrix and a meta table from a time-series omics
experiment (e.g. proteomics metabolomics, transcriptomics, etc.). In the
data matrix, every row is a feature (protein, metabolite, gene) and
every column is a sample, taken at a specific time, and the values are
the measurement values (e.g. intensities, etc.).

In the meta table, every row describes a sample (column of the data
matrix), and the columns contain the different meta information, such as
as the time when the sample was taken.

If you want to achieve one or more of the following things:

- Find out which hyperparameters, such as degree of freedom, are needed
  for the limma spline analysis

- Find out which of the features (rows of the data matrix) (e.g. which
  of the proteins, metabolites, etc.) are changed significantly over the
  time and get a p-value for them.

- Cluster the hits (the significant features) based on their temporal
  pattern.

- Run GSEA (gene set enrichment analysis) with the clustered hits.

- Automatically generate a reports, showing all results

Then the `SplineOmics` package could be of interest to you. This package
finds the hits by applying splines (piece wise polynomial curves) to the
time-series omics data and using the R package `limma` to find out, in a
nutshell, which of the features have splines that are significantly
different from horizontal (which would mean that nothing happens over
the time, therefore no significant change).

### Summary

With `SplineOmics`, you can:

- **Explore Various Hyperparameters:** The `limma_hyperparams_screen()`
  function offers a comprehensive way to test combinations of
  hyperparameters, such as different datasets, `limma` design formulas,
  degrees of freedom, and p-value thresholds. This enables users to
  evaluate the impact of various settings on their analysis results
  systematically.

- **Perform Limma Spline Analysis:** Once the optimal hyperparameters
  are identified, `run_limma_splines()` performs the limma analysis
  using splines.

- **Cluster Significant Features:** The `cluster_hits()` function goes a
  step further by clustering the significant features (hits) identified
  in the spline analysis. It organizes these features into meaningful
  groups (clusters) and generates a comprehensive report, facilitating
  the interpretation and communication of the results.

- **Run GSEA with the clustered hits:** The `gsea_report()` function
  takes the clustered hits and performs a gene set enrichment analysis,
  using clusterProfiler, with them.

## Installation

Follow these steps to install the `SplineOmics` package from its GitHub
repository into your R environment.

#### Prerequisites

- Ensure **R** is installed on your system. If not, download and install
  it from [CRAN](https://cran.r-project.org/).
- **RStudio** is recommended for a more user-friendly experience with R.
  Download and install RStudio from
  [rstudio.com](https://www.rstudio.com/products/rstudio/download/).

#### Installation Steps

1.  **Open RStudio** or your R console.

2.  **Install and load `devtools`, and install `SplineOmics` from
    GitHub**:

    Copy and paste the following code block into your R console. This
    will check if the `devtools` package is installed, install it if it
    is not, then use it to install the `SplineOmics` package from
    GitHub.

    ``` r
    # Check if devtools is installed; if not, install it
    if (!requireNamespace("devtools", quietly = TRUE)) {
      install.packages("devtools")
    }

    # Load devtools package
    library(devtools)

    Sys.setenv(GITHUB_PAT = "your_GitHub_PAT")

    # Install the SplineOmics package from GitHub
    devtools::install_github("csbg/SplineOmics")
    ```

3.  **Load the `SplineOmics` package**:

    Once the installation is complete, load the `SplineOmics` package
    into your R session or script to start using it:

    ``` r
    library(SplineOmics)
    ```

#### Troubleshooting

If you encounter errors related to dependencies or package versions
during installation, try updating your R and RStudio to the latest
versions and repeat the installation steps.

For issues specifically related to the `SplineOmics` package, check the
Issues section of the GitHub repository for similar problems or to post
a new issue.

## Usage

The `SplineOmics` has following functions available (for a more detailed
description, see the usage vignette):

### extract_data()

Extract the data section and the feature description out of an initial
Excel file:

``` r
data <- extract_data(excel_table,
                     "Unique identifier")
```

### explore_data()

Generate an HTML with the exploratory data analysis (and additionally
return the plots in R):

``` r
plots <- explore_data(data,
                      meta,
                      condition,
                      report_info,
                      meta_batch_column,
                      report_dir)
```

### limma_hyperparams_screen()

Allows to systematically test different hyperparameters for the
subsequent timeseries analysis with limma and splines. Test different
spline parameters, limma design formulas, verions of data (outliers
removed vs. not removed), etc. in a semi-combinatorial fashion.

``` r
result <- limma_hyperparams_screen(datas,
                                   datas_descr,
                                   metas,
                                   designs,
                                   condition,
                                   spline_test_configs,
                                   report_info,
                                   report_dir,
                                   pthresholds,
                                   meta_batch_column)
```

### run_limma_splines()

Get the limma topTable, with the adjusted p-values for each feature, by
analysing the time-series data with splines and limma:

``` r
top_tables <- run_limma_splines(data, 
                            meta, 
                            design, 
                            spline_params,
                            condition)
```

### cluster_hits()

Cluster the significant features (hits) via hierarchical clustering:

``` r
clustering_results <- cluster_hits(top_tables, 
                                   data, 
                                   meta, 
                                   design,
                                   condition, 
                                   spline_params,
                                   adj_pthresholds,
                                   clusters,
                                   report_info,
                                   meta_batch_column,
                                   report_dir)
```

### run_gsea()

Perform a gene set enrichment analysis (gsea), using clusterProfiler,
with the clustered hits:

``` r
run_gsea <- function(levels_clustered_hits,
                     genes,
                     databases,
                     report_info,
                     params = NA,
                     plot_titles = NA,
                     background = NULL,
                     report_dir = here::here()) {
```

## Docker Container

To facilitate reproducible analysis, we provide a Docker container that
encapsulates the necessary environment and dependencies. Follow the
instructions below to pull the Docker container and run your analysis.

### Pulling the Docker Container

You can pull the Docker container of the desired version of the
`SplineOmics` package from Docker Hub using the following command (here
it downloads version 0.1.0):

``` sh
docker pull thomasrauter/splineomics:0.1.0
```

### Running the Analysis Script

To run the analysis script, you need to mount your data files and the
parameter file into the container. Ensure that you have the following
files prepared:

data.xlsx: Your data file. meta.xlsx: Your metadata file. params.json: A
JSON file with the analysis parameters.

Here is the bash command to run the analysis script:

``` sh
docker run -it --rm \
    -v /path/to/data.xlsx:/workspace/data.xlsx \
    -v /path/to/meta.xlsx:/workspace/meta.xlsx \
    -v /path/to/params.json:/workspace/params.json \
    -v /path/to/output:/output \
    thomasrauter/splineomics Rscript /workspace/run_analysis.R
```

### Example params.json File

Below is an example of the params.json file that you need to provide:

``` json
{
  "report_info": {
    "omics_data_type": "PPTX",
    "data_description": "Old phosphoproteomics data with the missing two samples",
    "data_collection_date": "February 2024",
    "analyst_name": "Thomas Rauter",
    "contact_info": "thomas.rauter@plus.ac.at",
    "project_name": "DGTX"
  },
  "condition": "Phase",
  "meta_batch_column": "Reactor",
  "design": "~ 1 + Phase*X + Reactor",
  "spline_params": {
    "spline_type": ["n"],
    "dof": [2]
  },
  "adj_pthresholds": [0.05, 0.05],
  "clusters": [2, 2]
}
```

### Running the Analysis

Ensure you replace /path/to/ with the actual paths to your files. The
analysis results will be saved in the /path/to/output directory.

## Dependencies

The `SplineOmics` package relies on several other R packages for its
functionality. Below is a list of dependencies that will automatically
be installed along with `SplineOmics`. If you already have these
packages installed, ensure they are up to date to avoid any
compatibility issues.

- **ComplexHeatmap**: For creating complex heatmaps with advanced
  features.
- **RColorBrewer**: For providing color palettes for data visualization.
- **base64enc**: For encoding/decoding base64.
- **circlize**: For circular visualization of data.
- **cluster**: For clustering algorithms such as hierarchical
  clustering.
- **clusterProfiler**: For functional enrichment analysis.
- **data.table**: For high-performance data manipulation.
- **dendextend**: For extending `dendrogram` objects in R, allowing for
  easier manipulation of dendrograms.
- **dplyr**: For data manipulation.
- **fs**: For file system operations.
- **ggplot2**: For creating elegant data visualizations using the
  grammar of graphics.
- **ggrepel**: For better label placement in ggplot2.
- **here**: For constructing paths to your project’s files.
- **htmltools**: For HTML rendering and output.
- **kableExtra**: For producing beautiful tables in R Markdown
  documents.
- **knitr**: For dynamic report generation in R.
- **limma**: For linear models for microarray data.
- **magrittr**: For piping operators.
- **openxlsx**: For reading, writing, and editing xlsx files.
- **patchwork**: For combining multiple ggplot objects into a single
  plot.
- **pheatmap**: For creating pretty heatmaps.
- **progress**: For adding progress bars to your loops and apply
  functions.
- **purrr**: For functional programming tools.
- **ragg**: For creating high-quality images for graphics devices.
- **readr**: For reading rectangular data.
- **reshape2**: For flexible reshaping of data.
- **rlang**: For tools to work with core language features of R and R’s
  base types.
- **scales**: For scale functions for visualization.
- **stringr**: For simplifying the manipulation of strings.
- **tibble**: For creating tidy data frames that are easy to work with.
- **tidyr**: For tidying your data.
- **viridis**: For colorblind-friendly color maps.

### R Version

- Recommended: R 4.3.3 or higher
  - Note: This project was developed using R 4.3.3. While it should be
    compatible with newer versions, this is the version guaranteed to
    work as tested.

## Getting Help

If you encounter a bug or have a suggestion for improving the
`SplineOmics` package, we encourage you to [open an
issue](https://github.com/csbg/SplineOmics/issues) on our GitHub
repository. Before opening a new issue, please check to see if your
question or bug has already been reported by another user. This helps
avoid duplicate reports and ensures that we can address problems
efficiently.

For more detailed questions, discussions, or contributions regarding the
package’s use and development, please refer to the [GitHub
Discussions](https://github.com/csbg/SplineOmics/discussions) page for
`SplineOmics`. This forum is a great place to ask for help, share your
experiences, and connect with the community.

Thank you for using and contributing to the development of
`SplineOmics`!

## Contributing

We welcome contributions to the `SplineOmics` package! Whether you’re
interested in fixing bugs, adding new features, or improving
documentation, your help is greatly appreciated.

Here’s how you can contribute:

1.  **Report a Bug or Request a Feature:** If you encounter a bug or
    have an idea for a new feature, please [open an
    issue](https://github.com/csbg/SplineOmics/issues) on our GitHub
    repository. Before opening a new issue, check to see if the issue
    has already been reported or the feature requested by another user.

2.  **Submit a Pull Request:** If you’ve developed a bug fix or a new
    feature that you’d like to share, submit a pull request. Here are
    the steps:

    - Fork the repository.
    - Create a new branch in your fork for your contributions.
    - Commit your changes to the branch.
    - Push the branch to your fork.
    - Submit a pull request to the `SplineOmics` repository from your
      fork and branch.
    - Please describe your changes clearly in the pull request
      description and reference any related issues.

3.  **Improve Documentation:** Good documentation is crucial for any
    project. If you notice missing or incorrect documentation, please
    feel free to contribute.

Please adhere to our [Code of Conduct](CODE_OF_CONDUCT.md) in all your
interactions with the project.

Thank you for considering contributing to `SplineOmics`. Your efforts
are what make the open-source community a fantastic place to learn,
inspire, and create.

## License

This package is licensed under the MIT License with additional terms.
Please see the [LICENSE](./docs/LICENSE) file for full terms and
conditions.

© 2024 Thomas Rauter. All rights reserved.

## Citation

The `SplineOmics` package is currently not published in a peer-reviewed
scientific journal or similar outlet, and as such, there is no formal
citation requirement associated with its use. You are free to use
`SplineOmics` without citing it. However, we appreciate acknowledgements
in your projects or publications that benefit from this package. Also,
if you like the package, consider giving the GitHub Repository a star.
Your support helps us in the continued development and improvement of
`SplineOmics`. Thank you for using our package!

The `SplineOmics` relies substantially on following R libraries in
particular, which should be cited in your publication:

- limma –\> Ritchie, M.E., Phipson, B., Wu, D., Hu, Y., Law, C.W., Shi,
  W., and Smyth, G.K. (2015). limma powers differential expression
  analyses for RNA-sequencing and microarray studies. Nucleic Acids
  Research 43(7), e47. (DOI:
  [10.1093/nar/gkv007](https://doi.org/10.1093/nar/gkv007))

- ggplot2 –\> H. Wickham. ggplot2: Elegant Graphics for Data Analysis.
  Springer-Verlag New York, 2016.

- dplyr –\> Wickham H, François R, Henry L, Müller K, Vaughan D (2023).
  *dplyr: A Grammar of Data Manipulation*. R package version 1.1.4,
  <https://CRAN.R-project.org/package=dplyr>.

- tidyr –\> Wickham H, Vaughan D, Girlich M (2024). *tidyr: Tidy Messy
  Data*. R package version 1.3.1,
  <https://CRAN.R-project.org/package=tidyr>.

- ComplexHeatmap –\> Gu, Z. (2016) Complex heatmaps reveal patterns and
  correlations in multidimensional genomic data. Bioinformatics.

- pheatmap –\> Kolde R (2019). *pheatmap: Pretty Heatmaps*. R package
  version 1.0.12, <https://CRAN.R-project.org/package=pheatmap>.

- purrr –\> Wickham H, Henry L (2023). *purrr: Functional Programming
  Tools*. R package version 1.0.2,
  <https://CRAN.R-project.org/package=purrr>.

- tibble –\> Müller K, Wickham H (2023). *tibble: Simple Data Frames*. R
  package version 3.2.1, <https://CRAN.R-project.org/package=tibble>.
