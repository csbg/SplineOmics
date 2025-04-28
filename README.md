
# SplineOmics

![Version](https://img.shields.io/badge/version-0.1.2-blue) [![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE)
![Maintained?
Yes](https://img.shields.io/badge/Maintained%3F-yes-brightgreen.svg) ![R
CMD
Check](https://img.shields.io/badge/R%20CMD%20check-passed-brightgreen)
[![Docker](https://img.shields.io/badge/docker-pull-blue)](https://hub.docker.com/r/thomasrauter/splineomics)
![Dependencies](https://img.shields.io/badge/dependencies-19-blue)
![Platforms](https://img.shields.io/badge/platforms-all-brightgreen)

<img src="man/figures/hex_logo.png" align="left" width="65"/>

The R package `SplineOmics` finds the significant features (hits) of
time-series -omics data by using splines and `limma` for hypothesis
testing. It then clusters the hits based on the spline shape while
showing all results in summary HTML reports.

The graphical abstract below shows the full workflow streamlined by
`SplineOmics`:

<figure>
<img src="man/figures/SplineOmics_graphical_abstract.png"
alt="Graphical Abstract of SplineOmics Workflow" />
<figcaption aria-hidden="true">Graphical Abstract of SplineOmics
Workflow</figcaption>
</figure>

## Table of Contents

- [📘 Introduction](#-introduction)
- [🔧 Installation](#-installation)
  - [🐳 Docker Container](#-docker-container)
- [▶️ Usage](#-usage)
  - [🎓 Tutorial](#-tutorial)
  - [📋 Details](#-details)
  - [🧬 RNA-seq and Glycan Data](#-rna-seq-and-glycan-data)
- [📦 Dependencies](#-dependencies)
- [📚 Further Reading](#-further-reading)
- [❓ Getting Help](#-getting-help)
- [🤝 Contributing](#-contributing)
- [💬 Feedback](#-feedback)
- [📜 License](#-license)
- [🎓 Citation](#-citation)
- [🌟 Contributors](#-contributors)
- [🙏 Acknowledgements](#-acknowledgements)

## 📘 Introduction

Welcome to `SplineOmics`, an R package designed to streamline the
analysis of -omics time-series data, followed by automated HTML report
generation.

### Is the SplineOmics package of use for me?

If you have -omics data over time, the package will help you to run
`limma` with splines, decide on which parameters to use, perform the
clustering, run GSEA and show result plots in HTML reports. Any
time-series data that is a valid input to the `limma` package is also a
valid input to the `SplineOmics` package (such as transcriptomics,
proteomics, phosphoproteomics, metabolomics, glycan fractional
abundances, etc.).

### What do I need precisely?

1.  **Data**: A data matrix where each row is a feature (e.g., protein,
    metabolite, etc.) and each column is a sample taken at a specific
    time. The data must have no NA values, should have normally
    distributed features and no dependence between the samples.

2.  **Meta**: A table with metadata on the columns/samples of the data
    matrix (e.g., batch, time point, etc.)

3.  **Annotation** (optional): A table with identifiers on the
    rows/features of the data matrix (e.g., gene and protein name).

### Capabilities

With `SplineOmics`, you can:

- **Automatically perform exploratory data analysis:**

  The `explore_data()` function generates an HTML report, containing
  various plots, such as density, PCA, and correlation heatmap plots
  ([example
  report](https://csbg.github.io/SplineOmics_html_reports/explore_data_PTX.html)).

- **Explore various limma splines hyperparameters:**

  Test combinations of hyperparameters, such as different datasets,
  `limma` design formulas, degrees of freedom, p-value thresholds, etc.,
  using the `screen_limma_hyperparams()` function ([example
  report](https://csbg.github.io/SplineOmics_html_reports/Data_1_Design_1_vs_Data_1_Design_2_PTX.html)
  (along with the
  [encoding](https://csbg.github.io/SplineOmics_html_reports/hyperparams_screen_meta_table.html))).

- **Perform limma spline analysis:**

  Use the `run_limma_splines()` function to perform the `limma` analysis
  with splines once the optimal hyperparameters are identified ([example
  report](https://csbg.github.io/SplineOmics_html_reports/create_limma_report_PTX.html)).

- **Cluster significant features:**

  Cluster the significant features (hits) identified in the spline
  analysis with the `cluster_hits()` function ([example
  report](https://csbg.github.io/SplineOmics_html_reports/report_clustered_hits_PTX.html)).

- **Run GSEA with clustered hits:**

  Perform gene set enrichment analysis (GSEA) using the clustered hits
  with the `create_gsea_report()` function ([example
  report](https://csbg.github.io/SplineOmics_html_reports/create_gsea_report_PTX.html)).

## 🔧 Installation

Follow the steps below to install the `SplineOmics` package from the
GitHub repository into your R environment.

> **Note** Carefully read the terminal messages of the installations. It
> can happen that installations fail due to missing dependencies, which
> you then must resolve using other commands not necessarily written
> down here.

#### Prerequisites

- Ensure **R** is installed on your system. If not, download and install
  it from [CRAN](https://cran.r-project.org/).
- **RStudio** is recommended for a more user-friendly experience with R.
  Download and install RStudio from
  [posit.co](https://posit.co/download/rstudio-desktop/).

#### Installation Steps

> **Note for Windows Users:**

Note that some installation paths potentially are not writable on
**Windows**. Therefore, it can be necessary to set up a library path and
use that path for the installations:

``` r
custom_lib_path <- "C:/Rlibs"  # Replace with your desired path

# Create the directory if it doesn't exist
if (!dir.exists(custom_lib_path)) {
    dir.create(
      custom_lib_path,
      showWarnings = FALSE,
      recursive = TRUE
      )
}

# Set the library path to include the new directory
.libPaths(c(custom_lib_path, .libPaths()))

# Check if the new library path is added successfully
if (custom_lib_path %in% .libPaths()) {
    message("Library path set to: ", custom_lib_path)
} else {
    stop("Failed to set library path.")
}
```

Alternatively, you can run `RStudio` as administrator once for the
installation (which is however generally not recommended, because it is
a security risk).

> **Additional note for Windows Users:**

During the installation on Windows, you might see a message indicating
that Rtools is not installed, which is typically required for compiling
R packages from source. However, for this installation, Rtools is not
necessary, and you can safely ignore this message.

1.  **Open RStudio** or your R console.

2.  **Install `BiocManager`** for Bioconductor dependencies (if not
    already installed)

``` r
install.packages(
  "BiocManager"
  # lib = custom_lib_path 
)
```

3.  **Install Bioconductor dependencies** separately using `BiocManager`

``` r
BiocManager::install(
  c("ComplexHeatmap", "limma", "variancePartition")
  # force = TRUE   # when encountering issues
  # lib = custom_lib_path 
)
```

4.  **Install** the **`remotes`** package for GitHub downloads (if not
    already installed)

``` r
install.packages(
  "remotes"
  # lib = custom_lib_path
)
```

5.  **Install** the **`SplineOmics`** package from GitHub and all its
    non-Bioconductor dependencies, using `remotes`

``` r
remotes::install_github(
  "csbg/SplineOmics",   # GitHub repository
  ref = "0.1.2",        # Specify the tag to install
  dependencies = TRUE,  # Install all dependencies
  upgrade = "always"    # Always upgrade dependencies
  # force = TRUE        # when encountering issues
  # lib = custom_lib_path 
)
```

6.  **Verify** the **installation** of the `SplineOmics` package

``` r
# Verify the installation of the SplineOmics package
if ("SplineOmics" %in% rownames(installed.packages())) {
  message("SplineOmics was installed successfully.")
} else {
  message("SplineOmics installation failed. Please check for errors during installation.")
}
```

#### Troubleshooting

If you encounter errors related to dependencies or package versions
during installation, try updating your R and RStudio to the latest
versions and repeat the installation steps.

For issues specifically related to the `SplineOmics` package, check the
[Issues section](https://github.com/csbg/SplineOmics/issues) of the
GitHub repository for similar problems or to post a new issue.

### 🐳 Docker Container

Alternatively, you can run your analysis in a `Docker` container. The
underlying `Docker` image encapsulates the `SplineOmics` package
together with the necessary environment and dependencies. This ensures
higher levels of reproducibility because the analysis is carried out in
a consistent environment, independent of the operating system and its
custom configurations.

Please note that you must have the `Docker Engine` installed on your
machine. For instructions on how to install it, consult the official
[Docker Engine installation
guide](https://docs.docker.com/engine/install/).

More information about `Docker` containers can be found on the [official
Docker page](https://www.docker.com/resources/what-container/).

For instructions on downloading the image of the `SplineOmics` package
and running the container, please refer to the [Docker
instructions](https://csbg.github.io/SplineOmics/articles/Docker-instructions.html).

#### Troubleshooting

If you face “permission denied” issues on Linux distributions, check
this
[vignette](https://csbg.github.io/SplineOmics/articles/Docker_permission_denied.html).

## ▶️ Usage

### 🎓 Tutorial

[This
tutorial](https://csbg.github.io/SplineOmics/articles/get-started.html)
covers a real CHO cell time-series proteomics example from start to end.

To open an R Markdown file of the **tutorial** in `RStudio`, run:

``` r
library(SplineOmics)
open_tutorial()  
```

### 📋 Details

A detailed description of all arguments and outputs of all the functions
in the package (exported and internal functions) can be found
[here](https://csbg.github.io/SplineOmics/reference/).

#### Design `limma` design formula

A quick guide on how to design a `limma` design formula can be found
[here](https://csbg.github.io/SplineOmics/articles/design_limma_design_formula.html).

An explanation of the three different `limma` results is
[here](https://csbg.github.io/SplineOmics/articles/limma_result_categories.html).

### 🧬 RNA-seq and Glycan Data

#### RNA-seq data

Transcriptomics data must be preprocessed for `limma`. You need to
provide an appropriate object, such as a `voom` object, in the
`rna_seq_data` argument of the `SplineOmics` object (see
[documentation](https://csbg.github.io/SplineOmics/reference/create_splineomics.html)).
Along with this, the normalized matrix (e.g., the `$E` slot of the
`voom` object) must be passed to the `data` argument. This allows
flexibility in preprocessing; you can use any method you prefer as long
as the final object and matrix are compatible with limma. One way to
preprocess your RNA-seq data is by using the `preprocess_rna_seq_data()`
function included in the `SplineOmics` package (see
[documentation](https://csbg.github.io/SplineOmics/reference/preprocess_rna_seq_data.html)).

[Here](https://csbg.github.io/SplineOmics/articles/RNA-seq%20analysis.html)
you can find an example analysis of RNA-seq data with the SplineOmics
package.

#### Glycan fractional abundance data

The glycan fractional abundance data matrix, where each row represents a
type of glycan and the columns correspond to timepoints, must be
transformed before analysis. This preprocessing step is essential due to
the compositional nature of the data. In compositional data, an increase
in the abundance of one component (glycan) necessarily results in a
decrease in others, introducing a dependency among the variables that
can bias the analysis. One way to address this issue is by applying the
Centered Log Ratio (CLR) transformation to the data with the clr
function from the compositions package:

``` r
library(compositions)
clr_transformed_data <- clr(data_matrix)  # use as SplineOmics input
```

The results from clr transformed data can be harder to understand and
interpret however. If you prefer ease of interpretation and are fine
that the results contain some artifacts due to the compositional nature
of the data, log2 transform your data instead and use that as input for
the `SplineOmics` package.

``` r
log2_transformed_data <- log2(data_matrix)  # use as SplineOmics input
```

## 📦 Dependencies

The `SplineOmics` package relies on several other R packages for its
functionality. Below is a list of dependencies that will automatically
be installed along with `SplineOmics`. If you already have these
packages installed, ensure they are up to date to avoid any
compatibility issues.

- **ComplexHeatmap** (\>= 2.18.0): For creating complex heatmaps with
  advanced features.
- **base64enc** (\>= 0.1-3): For encoding/decoding base64.
- **dendextend** (\>= 1.17.1): For extending `dendrogram` objects,
  allowing for easier manipulation of dendrograms.
- **dplyr** (\>= 1.1.4): For data manipulation.
- **ggplot2** (\>= 3.5.1): For creating elegant data visualizations
  using the grammar of graphics.
- **ggrepel** (\>= 0.9.5): For better label placement in `ggplot2`.
- **here** (\>= 1.0.1): For constructing paths to your project’s files.
- **limma** (\>= 3.58.1): For linear models in microarray and RNA-seq
  analysis.
- **openxlsx** (\>= 4.2.6.1): For reading, writing, and editing `.xlsx`
  files.
- **patchwork** (\>= 1.2.0): For combining multiple `ggplot` objects
  into a single plot.
- **pheatmap** (\>= 1.0.12): For creating aesthetically pleasing
  heatmaps.
- **progress** (\>= 1.2.3): For adding progress bars to loops and apply
  functions.
- **purrr** (\>= 1.0.2): For functional programming tools.
- **rlang** (\>= 1.1.4): For working with core language features of R.
- **scales** (\>= 1.3.0): For scale functions in visualizations.
- **svglite** (\>= 2.1.3): For creating high-quality vector graphics
  (SVG).
- **tibble** (\>= 3.2.1): For creating tidy data frames.
- **tidyr** (\>= 1.3.1): For tidying data.
- **zip** (\>= 2.3.1): For compressing and combining files into zip
  archives.

### Optional Dependencies

These packages are optional and are only needed for specific
functionality:

- **edgeR** (\>= 4.0.16): For preprocessing RNA-seq data in the
  `preprocess_rna_seq_data()` function.
- **clusterProfiler** (\>= 4.10.1): For the `run_gsea()` function (gene
  set enrichment analysis).
- **rstudioapi** (\>= 0.16.0): For the `open_tutorial()` function.

### R Version

- Recommended: R 4.3.3 or higher

## 📚 Further Reading

For those interested in gaining a deeper understanding of the
methodologies used in the `SplineOmics` package, here are some
recommended publications:

- **Splines**: To learn more about splines, you can refer to this
  [review](https://doi.org/10.1186/s12874-019-0666-3).

- **limma**: To read about the `limma` R package, you can refer to this
  [publication](https://doi.org/10.1093/nar/gkv007).

- **PCA**: To learn more about PCA, download and read this
  [document](https://github.com/csbg/SplineOmics/raw/main/docs/Points_of_Significance_PCA.pdf).

- **Hierarchical clustering**: To get information about hierarchical
  clustering, you can refer to this [web
  article](https://towardsdatascience.com/understanding-the-concept-of-hierarchical-clustering-technique-c6e8243758ec).

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
`SplineOmics`.

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

## 🙏 Acknowledgements

This work was carried out in the context of the DigiTherapeutX project,
which was funded by the Austrian Science Fund (FWF). The research was
conducted under the supervision of Prof. Nikolaus Fortelny, who leads
the Computational Systems Biology working group at the Paris Lodron
University of Salzburg, Austria. You can find more information about
Prof. Fortelny’s research group
[here](https://www.plus.ac.at/biowissenschaften/der-fachbereich/arbeitsgruppen/fortelny/).
