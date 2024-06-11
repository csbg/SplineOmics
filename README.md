
# SplineOmics <img src="man/figures/hex_logo.png" style="float: right; width: 150px; margin-left: 300px; vertical-align: middle;"/>

![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
![Maintained?
Yes](https://img.shields.io/badge/Maintained%3F-yes-brightgreen.svg)

The R package `SplineOmics` gets the significant features (hits) of
time-series omics data by using splines and limma for hypothesis testing
and clusters the hits based on the spline shape, showing all results in
HTML reports.

## Table of Contents

- [📘 Introduction](#introduction)
- [🔧 Installation](#installation)
- [📂 Usage](#usage)
  - [📖 Tutorial](#tutorial)
  - [📊 Overview](#overview)
    - [📄 extract_data](#extract_data)
    - [🔍 explore_data](#explore_data)
    - [🧮️ screen_limma_hyperparams](#screen_limma_hyperparams)
    - [⚙️ run_limma_splines](#run_limma_splines)
    - [📈 create_limma_report](#create_limma_report)
    - [🌐 cluster_hits](#cluster_hits)
    - [📥 download_enrichr_databases](#download_enrichr_databases)
    - [🔬 create_gsea_report](#create_gsea_report)
  - [🛠️ Functions in Depth](#functions-in-depth)
- [🐳 Docker Container](#docker-container)
- [📦 Dependencies](#dependencies)
- [❓ Getting Help](#getting-help)
- [🤝 Contributing](#contributing)
- [📜 License](#license)
- [🎓 Citation](#citation)

## 📘 Introduction

Welcome to `SplineOmics`, an R package designed to streamline the
analysis of omics time-series data, followed by automated HTML report
generation.

### Is the SplineOmics package of use for me?

If you have omics data over time, the package will help you to run
`limma` with splines, decide on which parameters to use, perform the
clustering, run GSEA and show result plots in HTML reports.

### What do I need precisely?

1.  A data matrix where each row is a feature (e.g., protein,
    metabolite, etc.) and each column is a sample taken at a specific
    time.

2.  A table with metadata on the rows/features (e.g., gene and protein
    name)

3.  A table with metadata on the columns/samples (e.g., reactor, time
    point, etc.)

### Capabilities

With `SplineOmics`, you can:

- **Automatically perform exploratory data analysis (EDA):**

  The `explore_data()` function generates an HTML report, containing
  various EDA plots, such as densitiy, PCA, and correlation heatmap
  plots.

- **Explore various limma splines hyperparameters:**

  Test combinations of hyperparameters, such as different datasets,
  `limma` design formulas, degrees of freedom, and p-value thresholds
  using the `screen_limma_hyperparams()` function.

- **Perform limma spline analysis:**

  Use the `run_limma_splines()` function to perform the `limma` analysis
  with splines once the optimal hyperparameters are identified.

- **Cluster significant features:**

  Cluster the significant features (hits) identified in the spline
  analysis with the `cluster_hits()` function.

- **Run GSEA with clustered hits:**

  Perform gene set enrichment analysis using the clustered hits with the
  `create_gsea_report()` function.

- **Generate reports:**

  Automatically generate reports to showcase all results.

## 🔧 Installation

Follow these steps to install the `SplineOmics` package from its GitHub
repository into your R environment.

#### Prerequisites

- Ensure **R** is installed on your system. If not, download and install
  it from [CRAN](https://cran.r-project.org/).
- **RStudio** is recommended for a more user-friendly experience with R.
  Download and install RStudio from
  [posit.co](https://posit.co/download/rstudio-desktop/).

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
[Issues section](https://github.com/%3Cuser%3E/%3Crepo%3E/issues) of the
GitHub repository for similar problems or to post a new issue.

## 📂 Usage

### 📖 Tutorial

[This
tutorial](https://raw.githubusercontent.com/csbg/SplineOmics/main/doc/get-started.html)
covers a real CHO cell time-series proteomics example from start to the
end.

### 🛠 Functions in Depth

A detailed description of all arguments and outputs of all the available
package functions can be found
[here](https://raw.githubusercontent.com/csbg/SplineOmics/main/doc/functions-in-depth.html)

## 🐳 Docker Container

To facilitate reproducible analysis, we provide a Docker container that
encapsulates the necessary environment and dependencies. Follow the
instructions below to pull the Docker container and run your analysis.

### Pulling the Docker Container

You can pull the Docker container of the desired version of the
`SplineOmics` package from GitHub using the following command (here it
downloads version 0.1.0):

``` sh
docker pull ghcr.io/thomas-rauter/splineomics:0.1.0
```

### Running the Docker Container in Interactive Mode

To run the Docker container in interactive mode, you can use the
following command. This will start a container and open a bash shell,
allowing you to run your analysis interactively within the container.
This command needs to be run in a dir where the subdirs input and output
exist. Place your data and meta (and annotation) files in input, and
receive your output from the package in the output dir.

``` sh
docker run -it -d \
    -v $(pwd)/input:/home/rstudio/input \
    -v $(pwd)/output:/home/rstudio/output \
    -p 8888:8787 \
    -e PASSWORD=password \
    --name splineomics_container \
    thomasrauter/splineomics:0.1.0
```

Once the container is running, open a web browser and navigate to
<http://localhost:8888>. Use rstudio as the username and the password
you set with the -e PASSWORD=password option. As long as the container
is running, you can work on that localhost page with RStudio, where also
the `SplineOmics` package is installed.

Stop the container:

``` sh
docker stop splineomics_container
```

Start the container again:

``` sh
docker start splineomics_container
```

### Inspect Docker container installations

To see all the R packages and system installations that make up the
Docker container, you can run the following command in the terminal of
RStudio on your localhost page (<http://localhost:8888>). Because with
the above command the ‘/home/rstudio/output’ dir is mounted to your
local filesystem, this will make the installation log files available.

``` sh
cp -r /log home/rstudio/output
```

## 📦 Dependencies

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

## ❓ Getting Help

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

## 🤝 Contributing

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
    feature that you’d like to share, submit a pull request.

3.  **Improve Documentation:** Good documentation is crucial for any
    project. If you notice missing or incorrect documentation, please
    feel free to contribute.

Please adhere to our [Code of Conduct](CODE_OF_CONDUCT.md) in all your
interactions with the project.

Thank you for considering contributing to `SplineOmics`. Your efforts
are what make the open-source community a fantastic place to learn,
inspire, and create.

## 📜 License

This package is licensed under the MIT License: [LICENSE](./LICENSE)

© 2024 Thomas Rauter. All rights reserved.

## 🎓 Citation

The `SplineOmics` package is currently not published in a peer-reviewed
scientific journal or similar outlet. However, if this package helped
you in your work, consider citing this GitHub repository.

Also, if you like the package, consider giving the GitHub repository a
star. Your support helps us in the continued development and improvement
of `SplineOmics`. Thank you for using our package!
