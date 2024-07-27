
# SplineOmics <img src="man/figures/hex_logo.png" style="float: right; width: 150px; margin-left: 300px; vertical-align: middle;"/>

![Version](https://img.shields.io/badge/version-0.1.0-blue) [![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE)
![Maintained?
Yes](https://img.shields.io/badge/Maintained%3F-yes-brightgreen.svg) ![R
CMD
Check](https://img.shields.io/badge/R%20CMD%20check-passed-brightgreen)
[![Docker](https://img.shields.io/badge/docker-pull-blue)](https://hub.docker.com/r/thomasrauter/splineomics)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1234567.svg)](https://doi.org/10.5281/zenodo.1234567)

The R package `SplineOmics` finds the significant features (hits) of
time-series omics data by using splines and `limma` for hypothesis
testing and clusters the hits based on the spline shape while plotting
showing all results in summary HTML reports.

## Table of Contents

- [📘 Introduction](#-introduction)
- [🔧 Installation](#-installation)
  - [🐳 Docker Container](#-docker-container)
- [▶️ Usage](#-usage)
  - [Tutorial](#-tutorial)
  - [Functions in Depth](#-functions-in-depth)
- [📦 Dependencies](#-dependencies)
- [❓ Getting Help](#-getting-help)
- [🤝 Contributing](#-contributing)
- [💬 Feedback](#-feedback)
- [📜 License](#-license)
- [🎓 Citation](#-citation)
- [🌟 Contributors](#-contributors)

## 📘 Introduction

Welcome to `SplineOmics`, an R package designed to streamline the
analysis of -omics time-series data, followed by automated HTML report
generation.

### Is the SplineOmics package of use for me?

If you have -omics data over time, the package will help you to run
`limma` with splines, decide on which parameters to use, perform the
clustering, run GSEA and show result plots in HTML reports. Any data that is a valid input to the `limma` package is also a valid input to the `SplineOmics` package.

### What do I need precisely?

1.  A data matrix where each row is a feature (e.g., protein,
    metabolite, etc.) and each column is a sample taken at a specific
    time.

2.  A table with identifiers on the rows/features of the data matrix
    (e.g., gene and protein name)

3.  A table with metadata on the columns/samples of the data matrix
    (e.g., batch, time point, etc.)

### Capabilities

With `SplineOmics`, you can:

- **Automatically perform exploratory data analysis:**

  The `explore_data()` function generates an HTML report, containing
  various plots, such as density, PCA, and correlation heatmap plots.

- **Explore various limma splines hyperparameters:**

  Test combinations of hyperparameters, such as different datasets,
  `limma` design formulas, degrees of freedom, p-value thresholds, etc.,
  using the `screen_limma_hyperparams()` function.

- **Perform limma spline analysis:**

  Use the `run_limma_splines()` function to perform the `limma` analysis
  with splines once the optimal hyperparameters are identified.

- **Cluster significant features:**

  Cluster the significant features (hits) identified in the spline
  analysis with the `cluster_hits()` function.

- **Run GSEA with clustered hits:**

  Perform gene set enrichment analysis (GSEA) using the clustered hits
  with the `create_gsea_report()` function.

- **Generate reports:**

  Automatically generate reports to showcase all results.

## 🔧 Installation

Follow the steps below to install the `SplineOmics` package from the GitHub
repository into your R environment. 

#### Prerequisites

- Ensure **R** is installed on your system. If not, download and install
  it from [CRAN](https://cran.r-project.org/) (to open the links in this README section in a new tab, press the scroll wheel (middle button) of your mouse).
- **RStudio** is recommended for a more user-friendly experience with R.
  Download and install RStudio from
  [posit.co](https://posit.co/download/rstudio-desktop/).

#### Installation Steps

1.  **Open RStudio** or your R console.

2.  **Install** `SplineOmics` **from
    GitHub** with all dependencies.

Copy and paste the following code block into your R console or run it as
a script. 

> **Note for Windows Users:**  
> Please read the text below this code block before running it!

``` r
# Function to ensure a package is installed
ensure_installed <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Install packages if not already available
ensure_installed("BiocManager")
ensure_installed("remotes")

# Load packages
library(BiocManager)
library(remotes)

# Install Bioconductor dependencies
BiocManager::install(c(
  "ComplexHeatmap",
   "clusterProfiler",
   "limma"
  ), force = TRUE)

# Delete once repo is public
Sys.setenv(GITHUB_PAT = "your_GitHub_PAT")

# Install SplineOmics from GitHub
remotes::install_github(
  "csbg/SplineOmics@ad35d9aef2e8a8b19572c83ec771f4d92b343a4e",
  dependencies = TRUE,  # Install all dependencies
  force = TRUE,         # Force reinstallation
  upgrade = "always",   # Always upgrade dependencies
)

# Verify the installation
if ("SplineOmics" %in% rownames(installed.packages())) {
  message("SplineOmics installed successfully.")
} else {
  message("SplineOmics installation failed.")
}
```

Note that when some installation paths are not writable on **Windows**, it is necessary running `RStudio` as administrator once for the installation. Otherwise, set up a library path (code block below) for the installation and (re)run the code block above.

```r
# Create a directory for R libraries
dir.create("~/Rlibs", showWarnings = FALSE)

# Set the library path to include the new directory
.libPaths(c("~/Rlibs", .libPaths()))
```

3.  **Load the `SplineOmics` package**:

Once the installation is complete, load the `SplineOmics` package into
your R session or script to start using it:

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

### 🐳 Docker Container

Alternatively, you can run your analysis in a `Docker` container. The underlying `Docker` image
encapsulates the `SplineOmics` package together with the necessary
environment and dependencies. This ensures higher levels of reproducibility because the analysis is carried out in a consistent environment, independent of the operating system and its custom configurations.

More information about `Docker` containers can be found on the [official Docker page](https://www.docker.com/resources/what-container/).

For instructions on downloading the image of the `SplineOmics` package and running the container, please refer to the
[Docker instructions](https://raw.githubusercontent.com/csbg/SplineOmics/main/doc/Docker_instructions.html).

## ▶️ Usage

### Tutorial

[This
tutorial](https://raw.githubusercontent.com/csbg/SplineOmics/main/doc/get-started.html)
covers a real CHO cell time-series proteomics example from start to
end.

When you have the `SplineOmics` package installed, you can also run the
following commands in `RStudio` to start the tutorial as an interactive
demo:

``` r
library(SplineOmics)
i_demo()
```

This opens the R Markdown file of the demo in `RStudio`. Note that this
Markdown requires the readxls and the conflicted package in addition to
the dependencies of the `SplineOmics`\`package.

To get a template for your analysis, run:

``` r
library(SplineOmics)  # Once the package is loaded, this is not needed anymore.
template()
```

This opens an R Markdown file in `RStudio` that contains the "skeleton" of the whole workflow.

### Functions in Depth

A detailed description of all arguments and outputs of all the available
package functions can be found
[here](https://raw.githubusercontent.com/csbg/SplineOmics/main/doc/functions-in-depth.html).

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

## 📚 Further Reading

For those interested in gaining a deeper understanding of the methodologies used in the `SplineOmics` package, here are some recommended publications:

- **Splines**: To learn more about splines, you can refer to this [review](https://doi.org/10.1186/s12874-019-0666-3).

- **limma**: To read about the limma R package, you can refer to this [publication](10.1093/nar/gkv007).

- **Hierarchical clustering**: To get information about hierarchical clustering, you can refer to this [web article](https://towardsdatascience.com/understanding-the-concept-of-hierarchical-clustering-technique-c6e8243758ec).

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

## 💬 Feedback

We appreciate your feedback! Besides raising issues, you can provide
feedback in the following ways:

- **Direct Email**: Send your feedback directly to [Thomas
  Rauter](mailto:thomas.rauter@plus.ac.at).

- **Anonymous Feedback**: Use [this Google
  Form](https://forms.gle/jocMXSxLf3GrGBdT9) to provide anonymous
  feedback by answering questions.

Your feedback helps us improve the project and address any issues you
may encounter.

## 📜 License

This package is licensed under the MIT License: [LICENSE](./LICENSE)

© 2024 Thomas Rauter. All rights reserved.

## 🎓 Citation

The `SplineOmics` package is currently not published in a peer-reviewed
scientific journal or similar outlet. However, if this package helped
you in your work, consider citing this GitHub repository.

To cite this package, you can use the citation information provided in
the [`CITATION.cff`](./CITATION.cff) file.

You can also generate a citation in various formats using the
`CITATION.cff` file by visiting the top right of this repo and clicking
on the “Cite this repository” button.

Also, if you like the package, consider giving the GitHub repository a
star. Your support helps us in the continued development and improvement
of `SplineOmics`. Thank you for using our package!

## 🌟 Contributors

- [Thomas-Rauter](https://github.com/Thomas-Rauter) - 🚀 Wrote the
  package, developed the approach together with VSchaepertoens under
  guidance from nfortelny and skafdasschaf.
- [nfortelny](https://github.com/nfortelny) - 🧠 Principal Investigator,
  provided guidance and support for the overall approach.
- [skafdasschaf](https://github.com/skafdasschaf) - 🔧 Helped reviewing
  code, delivered improvement suggestions and scientific guidance to
  develop the approach.
- [VSchaepertoens](https://github.com/VSchaepertoens) - ✨ Developed one
  internal plotting function, as well as some code for the exploratory
  data analysis plots, and the overall approach together with
  Thomas-Rauter.
