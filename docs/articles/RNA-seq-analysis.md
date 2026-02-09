# RNA-seq-analysis

## About this vignette

This tutorial intends to showcase and explain the capabilities of the
**SplineOmics** package by walking through a real and complete RNA-seq
example, from start to finish. **SplineOmics** is explained in more
detail in the **get-started** vignette, where a proteomics example is
covered. This vignette is more focused on showing how RNA-seq data can
be used, and because of this, less details about the overall package are
provided here.

### Data Overview

This dataset originates from a **time-series RNA-seq experiment**
designed to study Chinese Hamster Ovary (CHO) cells. The experiment
involved cultivating cells in **eight bioreactors**, with **four
bioreactors** subjected to a temperature shift after 146 hours
(experimental condition) and the remaining **four bioreactors**
maintained without a temperature shift (control condition).

#### Timepoints

Samples were collected at **17 distinct time points** throughout the
experiment, specifically: `"72h"`, `"76h"`, `"96h"`, `"120h"`, `"124h"`,
`"144h"`, `"148h"`, `"152h"`, `"168h"`, `"192h"`, `"216h"`, `"220h"`,
`"240h"`, `"264h"`, `"268h"`, `"288h"`, and `"312h"` after cultivation
start. Each time point was sampled from all eight bioreactors, resulting
in a total of **136 samples**.

#### Additonal effects in the experiment: Reactor and Plate

In this experiment, there are two effects to consider: **Reactor** and
**Plate**.

1.  **Reactor**:
    - This refers to the different bioreactors used for cell
      cultivation, which can exhibit substantial variability.
    - Each reactor was assigned a single condition: either constant
      temperature or temperature-shifted. As a result, **condition and
      reactor are confounded**.
    - Reactor should not be treated as a fixed effect to simply remove
      its influence. Instead, it is treated as a **random effect**,
      which allows us to model its variability appropriately.
2.  **Plate**:
    - This refers to the two separate plates used during the RNA-seq
      analysis of the samples.
    - A fully random design was employed to distribute samples across
      the two plates, ensuring no bias in their assignment.
    - Plate is considered a **batch effect** with respect to the
      condition (constant temperature vs. temperature-shifted).

Since reactor is a blocked effect, the variability due to the reactors
cannot be directly separated from the condition. Instead, **linear mixed
models (LMMs)** are used to attribute the reactor as a random effect,
allowing us to account for its variability while isolating the effects
of the condition. This approach ensures that the analysis appropriately
handles the hierarchical structure of the data and avoids incorrect
conclusions.

In this vignette, we will demonstrate how to use linear mixed models to
address these challenges and properly account for both reactor and plate
effects.

#### Further info

The data matrix comprises **genes as rows** and **samples as columns**,
providing gene expression measurements for all time points. Each sample
was initially sequenced in **three technical replicates** across **two
NovaSeq X flow cells**. These technical replicates were collapsed to
generate the final dataset used in this analysis.

The goal of this experiment is to investigate the effect of a
temperature shift during CHO cell cultivation on gene expression
dynamics over time.

Note: This is not the original dataset, as it has not yet been published
at the time of this vignette’s creation. For demonstration purposes, the
genes have been randomly shuffled, and only a subset of the data has
been included to reduce the dataset size.

### Analysis Goals

The main objectives of this analysis are:

- **Identify genes with significant temporal changes**: Among the
  thousands of genes measured, the goal is to identify those that
  exhibit significant changes in expression over time.

- **Cluster genes based on temporal patterns**: Genes showing
  significant temporal changes (hits) will be grouped into clusters
  based on their time-dependent expression patterns.

- **Perform gene set enrichment analysis**: For each cluster, a gene set
  enrichment analysis will be conducted to identify whether specific
  biological pathways or processes are up- or downregulated in response
  to feeding and how these processes are influenced by the temperature
  shift.

- **Assess the impact of temperature shifts on temporal patterns**: The
  analysis will determine whether the temporal patterns of gene
  expression are affected by the temperature shift, i.e., whether gene
  expression dynamics differ over time under temperature shift
  conditions compared to controls.

#### Note

The documentation of all the **SplineOmics** package functions can be
viewed [here](https://csbg.github.io/SplineOmics/reference)

## Load the packages

``` r
library(SplineOmics)
library(readr) # For reading the meta CSV file
library(dplyr) # For data manipulation
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
library(knitr) # For Showing the head of the data and the meta tables.
```

## Load the files

``` r
data <- readRDS(xzfile(system.file(
    "extdata",
    "rna_seq_data.rds.xz",
    package = "SplineOmics"
)))

meta <- readr::read_csv(
    system.file(
        "extdata",
        "rna_seq_meta.csv",
        package = "SplineOmics"
    ),
    show_col_types = FALSE
)
```

### Show top rows of data

``` r
kable(
    head(data),
    format = "markdown"
)
```

|  | Sample1 | Sample57 | Sample80 | Sample68 | Sample29 | Sample121 | Sample129 | Sample46 | Sample49 | Sample75 | Sample103 | Sample41 | Sample52 | Sample120 | Sample107 | Sample23 | Sample71 | Sample56 | Sample12 | Sample86 | Sample5 | Sample98 | Sample37 | Sample22 | Sample94 | Sample62 | Sample106 | Sample21 | Sample119 | Sample11 | Sample40 | Sample114 | Sample28 | Sample77 | Sample82 | Sample20 | Sample31 | Sample125 | Sample110 | Sample92 | Sample27 | Sample61 | Sample105 | Sample70 | Sample55 | Sample102 | Sample39 | Sample48 | Sample116 | Sample60 | Sample16 | Sample132 | Sample6 | Sample124 | Sample67 | Sample113 | Sample115 | Sample9 | Sample15 | Sample45 | Sample74 | Sample123 | Sample85 | Sample24 | Sample26 | Sample59 | Sample65 | Sample88 | Sample118 | Sample101 | Sample18 | Sample136 | Sample72 | Sample35 | Sample128 | Sample87 | Sample73 | Sample100 | Sample84 | Sample135 | Sample25 | Sample34 | Sample81 | Sample131 | Sample97 | Sample64 | Sample17 | Sample134 | Sample93 | Sample8 | Sample104 | Sample44 | Sample30 | Sample10 | Sample109 | Sample91 | Sample4 | Sample7 | Sample36 | Sample111 | Sample54 | Sample99 | Sample38 | Sample133 | Sample51 | Sample58 | Sample14 | Sample43 | Sample53 | Sample63 | Sample130 | Sample112 | Sample3 | Sample33 | Sample127 | Sample42 | Sample117 | Sample122 | Sample108 | Sample90 | Sample50 | Sample76 | Sample126 | Sample19 | Sample96 | Sample79 | Sample83 | Sample47 | Sample2 | Sample32 | Sample13 | Sample69 | Sample95 | Sample78 | Sample66 | Sample89 |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| LOC113838358 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 4 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| LOC113836143 | 1227 | 1288 | 1303 | 1494 | 1380 | 1679 | 1620 | 1305 | 1243 | 1526 | 1889 | 1079 | 1386 | 1612 | 1880 | 1503 | 1262 | 1590 | 1260 | 1385 | 1631 | 1976 | 1348 | 1143 | 998 | 1903 | 1579 | 1253 | 2115 | 1218 | 1309 | 1977 | 1255 | 1522 | 1362 | 1465 | 1587 | 1893 | 2166 | 1776 | 1411 | 2065 | 1678 | 1763 | 1182 | 1881 | 989 | 1512 | 1828 | 1495 | 1306 | 2169 | 1511 | 1963 | 1304 | 1947 | 1816 | 1510 | 1528 | 1483 | 1197 | 1911 | 1652 | 1248 | 1214 | 1688 | 1279 | 1802 | 3340 | 2441 | 1502 | 1208 | 1179 | 1406 | 1870 | 1903 | 1628 | 2386 | 1651 | 1650 | 1439 | 1777 | 2009 | 2521 | 1839 | 1569 | 1550 | 3951 | 2018 | 1843 | 2807 | 1675 | 1488 | 1353 | 2205 | 2261 | 1634 | 1571 | 2161 | 2230 | 1512 | 2207 | 1426 | 3082 | 1445 | 1421 | 1840 | 1908 | 1779 | 1944 | 2394 | 3384 | 1771 | 1780 | 3673 | 1601 | 2786 | 2196 | 2504 | 2384 | 1543 | 1961 | 2433 | 1403 | 3074 | 1577 | 1846 | 2136 | 1981 | 1655 | 2160 | 1542 | 2777 | 1426 | 1723 | 2858 |
| Coro7_1 | 2 | 6 | 0 | 4 | 6 | 2 | 6 | 7 | 0 | 3 | 13 | 5 | 5 | 3 | 3 | 0 | 1 | 2 | 2 | 1 | 3 | 14 | 3 | 7 | 4 | 0 | 5 | 15 | 4 | 9 | 5 | 0 | 2 | 1 | 1 | 0 | 2 | 5 | 4 | 2 | 2 | 4 | 3 | 4 | 0 | 0 | 1 | 1 | 10 | 1 | 0 | 2 | 8 | 0 | 4 | 0 | 1 | 0 | 2 | 5 | 0 | 5 | 2 | 0 | 3 | 3 | 0 | 4 | 3 | 2 | 0 | 0 | 7 | 2 | 3 | 1 | 3 | 4 | 4 | 8 | 2 | 2 | 6 | 4 | 1 | 2 | 9 | 10 | 3 | 4 | 12 | 10 | 5 | 0 | 3 | 3 | 0 | 7 | 4 | 2 | 3 | 1 | 0 | 2 | 4 | 0 | 1 | 0 | 0 | 0 | 3 | 2 | 4 | 0 | 8 | 2 | 6 | 8 | 13 | 6 | 5 | 0 | 4 | 0 | 5 | 0 | 5 | 1 | 1 | 0 | 4 | 0 | 27 | 1 | 4 | 9 |
| LOC100765316_1 | 778 | 826 | 886 | 759 | 854 | 845 | 851 | 881 | 797 | 1097 | 1119 | 880 | 564 | 704 | 1244 | 941 | 781 | 583 | 727 | 958 | 987 | 1061 | 807 | 823 | 527 | 63 | 724 | 635 | 1070 | 752 | 702 | 893 | 746 | 859 | 786 | 777 | 897 | 931 | 959 | 1014 | 760 | 128 | 628 | 55 | 628 | 736 | 304 | 645 | 902 | 55 | 623 | 663 | 655 | 674 | 837 | 803 | 784 | 505 | 594 | 665 | 440 | 777 | 473 | 618 | 329 | 196 | 313 | 566 | 717 | 494 | 281 | 414 | 556 | 354 | 499 | 768 | 205 | 106 | 264 | 161 | 285 | 182 | 341 | 468 | 162 | 171 | 224 | 224 | 441 | 342 | 447 | 454 | 171 | 169 | 240 | 186 | 279 | 263 | 171 | 486 | 624 | 273 | 156 | 178 | 798 | 85 | 51 | 167 | 409 | 285 | 151 | 60 | 201 | 85 | 80 | 197 | 102 | 227 | 101 | 51 | 1013 | 70 | 46 | 167 | 80 | 212 | 71 | 42 | 65 | 26 | 33 | 60 | 67 | 132 | 489 | 62 |
| Slc7a14_1 | 6 | 1 | 7 | 10 | 2 | 5 | 12 | 3 | 1 | 5 | 2 | 4 | 2 | 2 | 7 | 1 | 10 | 5 | 5 | 8 | 4 | 9 | 11 | 2 | 4 | 19 | 19 | 7 | 7 | 19 | 12 | 10 | 12 | 4 | 7 | 2 | 11 | 5 | 21 | 12 | 7 | 15 | 8 | 8 | 3 | 10 | 8 | 10 | 11 | 9 | 10 | 0 | 3 | 12 | 0 | 14 | 3 | 11 | 6 | 4 | 13 | 8 | 11 | 5 | 8 | 24 | 12 | 29 | 17 | 19 | 12 | 0 | 20 | 22 | 12 | 16 | 14 | 16 | 18 | 8 | 22 | 11 | 4 | 24 | 26 | 7 | 11 | 43 | 13 | 5 | 19 | 7 | 13 | 13 | 30 | 15 | 12 | 25 | 13 | 30 | 8 | 4 | 9 | 32 | 5 | 16 | 2 | 36 | 8 | 18 | 15 | 24 | 8 | 15 | 12 | 5 | 10 | 9 | 22 | 15 | 0 | 13 | 12 | 10 | 29 | 8 | 4 | 8 | 29 | 12 | 14 | 15 | 36 | 8 | 7 | 7 |
| LOC100751623_1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

**data**: A numeric matrix where each row represents a gene (features)
and each column corresponds to a sample. The row names of the matrix
contain the gene identifiers, while the columns are aligned with the
sample metadata in meta. The matrix contains expression values for 136
samples.

Note that the study was conducted in a blinded manner, with samples
randomly distributed across two plates for RNA-seq analysis. As a
result, the sample numbers (e.g., 1, 2, 3, etc.) are not in sequential
order with respect to time, condition, or plate.

For the data analysis involving splines over time, it is essential to
sort the samples based on time to establish a valid temporal sequence.
Additionally, organizing the data in this way improves clarity and
ensures consistency. Within each time point, samples were further sorted
by condition (e.g., constant before `temp_shift`) and, subsequently, by
plate (e.g., `plate_1` before `plate_2`).

### Show top rows of meta

``` r
kable(
    head(meta),
    format = "markdown"
)
```

| SampleNr | Reactor | Time | Condition  | Plate   |
|---------:|:--------|-----:|:-----------|:--------|
|        1 | E13     |   72 | constant   | plate_1 |
|       57 | E15     |   72 | constant   | plate_1 |
|       80 | E17     |   72 | constant   | plate_1 |
|       68 | E19     |   72 | constant   | plate_1 |
|       29 | E14     |   72 | temp_shift | plate_1 |
|      121 | E16     |   72 | temp_shift | plate_2 |

**meta**: A data frame containing metadata information about the samples
in data. Each row in meta corresponds to a column in data, ensuring a
1:1 alignment between metadata entries and expression data samples. The
columns in meta describe various attributes of the samples, such as
SampleNr, Reactor, Time, Condition, and Plate.

## Preprocess the data

Filter out data rows (genes) with zero counts across all samples. This
step is a standard preprocessing procedure in RNA-seq data analysis, as
genes with zero counts in all samples provide no information for
downstream analyses.

``` r
rows_before <- nrow(data)

# Filter data rows
data <- data[rowSums(data) > 0, ]

rows_after <- nrow(data)
rows_removed <- rows_before - rows_after

cat(sprintf(
    "Rows before filtering: %d\nRows after filtering: %d\nRows removed: %d\n",
    rows_before,
    rows_after,
    rows_removed
))
#> Rows before filtering: 1000
#> Rows after filtering: 944
#> Rows removed: 56
```

## Perform EDA (exploratory data analysis)

``` r
report_info <- list(
    omics_data_type = "RNA",
    data_description = "RNA-seq data of CHO cells",
    data_collection_date = "December 2024",
    analyst_name = "Thomas Rauter",
    contact_info = "thomas.rauter@plus.ac.at",
    project_name = "DGTX"
)

report_dir <- file.path(
    tempdir(),
    "results",
    "explore_data"
)
```

``` r
splineomics <- SplineOmics::create_splineomics(
    data = data,
    meta = meta,
    report_info = report_info,
    condition = "Condition", # Column of meta that contains the levels.
    meta_batch_column = "Plate" # Remove batch effect for plotting.
)
```

``` r
plots <- SplineOmics::explore_data(
    splineomics = splineomics,
    report_dir = report_dir
)
```

[Here](https://drive.google.com/uc?export=view&id=1mE7FFXFO9ZZ907in3M5Cf1w7yRxRFHiy)
you can see the HTML report of the explore_data() function with the NOT
batch-corrected data, and
[here](https://drive.google.com/uc?export=view&id=1g_Lt4kO4H79GnVh_f7KfXqIzxfTnMrmq)
the report for the batch-corrected data.

## Run limma spline analysis

In this example, we are skipping finding the best hyperparameters with
the screen_limma_hyperparams() function, because we already have a clear
idea of what do use.

Lets define our parameters and put them into the `SplineOmics` object:

``` r
spline_params <- list(
    spline_type = c("n"), # natural cubic splines
    dof = c(3L) # Degree of freedom of 3 for the splines.
)

design <- "~ 1 + Condition*Time + Plate + (1|Reactor)"

splineomics <- SplineOmics::update_splineomics(
    splineomics = splineomics,
    data = data,
    design = design,
    # dream_params = dream_params,   if we would want to adjust dream()
    # means limma "borrows" statistical power from the other levels
    mode = "integrated",
    spline_params = spline_params,
    bp_cfg = c(
        n_cores = 1,
        blas_threads = 1
    )
)
```

Preprocess the RNA-seq data with
[`limma::voom`](https://rdrr.io/pkg/limma/man/voom.html). This method
transforms raw RNA-seq counts into log-counts per million (logCPM) while
modeling the mean-variance relationship. It assigns precision weights to
each observation, ensuring accurate linear modeling of RNA-seq data,
which often exhibits heteroscedasticity (varying variance across
expression levels). To normalize the data, TMM (Trimmed Mean of
M-values) normalization is applied using a `DGEList` object from the
**edgeR** package, correcting for library size differences and
compositional biases.

When random effects are included in the design, the function
automatically uses `voomWithDreamWeights` from the **variancePartition**
package instead. This method extends `voom` by incorporating random
effects into the model, allowing for precise handling of complex
experimental designs such as repeated measures or hierarchical
structures. The calculated weights account for both fixed and random
effects, providing robust results for differential expression analysis.

``` r
splineomics <- SplineOmics::preprocess_rna_seq_data(
    splineomics = splineomics
)
#> Warning: the 'findbars' function has moved to the reformulas package. Please update your imports, or ask an upstream package maintainter to do so.
#> This warning is displayed once per session.
```

You can customize the normalization method by providing a specific
normalization function through the `normalize_func` argument in the
[`preprocess_rna_seq_data()`](https://csbg.github.io/SplineOmics/reference/preprocess_rna_seq_data.md)
function. For details on how to use this feature, please refer to the
function documentation available under ‘References’ on the website.

Additionally, the use of
[`preprocess_rna_seq_data()`](https://csbg.github.io/SplineOmics/reference/preprocess_rna_seq_data.md)
is optional for your RNA-seq data. Alternatively, you can use the
[`limma::voom`](https://rdrr.io/pkg/limma/man/voom.html) function
directly and pass the resulting `voom` object to the `rna_seq_data`
argument of
[`create_splineomics()`](https://csbg.github.io/SplineOmics/reference/create_splineomics.md)
or
[`update_splineomics()`](https://csbg.github.io/SplineOmics/reference/update_splineomics.md).
Alongside this, you must pass the `$E` data matrix as the `data`
argument.

In general, as long as the `data` argument contains the actual data
matrix and the `rna_seq_data` argument contains an object compatible
with `limma`, your data will be correctly processed.

Run the
[`run_limma_splines()`](https://csbg.github.io/SplineOmics/reference/run_limma_splines.md)
function:

``` r
splineomics <- SplineOmics::run_limma_splines(
    splineomics = splineomics
)
#> Array weights are ignored at this stage for RNA-seq data, as they were already incorporated during the preprocessing step.
#> Warning: the 'nobars' function has moved to the reformulas package. Please update your imports, or ask an upstream package maintainter to do so.
#> This warning is displayed once per session.
```

The output of the function run_limma_splines() is a named list, where
each element is a specific “category” of results. Refer to [this
document](https://csbg.github.io/SplineOmics/articles/limma_result_categories.html)
for an explanation of the different result categories. Each of those
elements is a list, containing as elements the respective limma
topTables, either for each level or each comparison between two levels.

The element “time_effect” is a list, where each element is the topTable
where the p-value for each feature for the respective level are
reported.

The element “avrg_diff_conditions” is a list that contains as elements
the topTables, that represent the comparison of the average differences
of the levels.

The element “interaction_condition_time” is a list that contains as
elements the topTables, that represent the interaction between the
levels (which includes both time and the average differences)

## Build limma report

The topTables of all three limma result categories can be used to
generate p-value histograms an volcano plots.

``` r
plots <- SplineOmics::create_limma_report(
    splineomics = splineomics,
    report_dir = withr::local_tempdir()
)
```

You can view the generated analysis report of the create_limma_report
function
[here](https://drive.google.com/uc?export=view&id=1hvKloiaHQA5CO4Y6io8qC2lOpBJ54ypJ).

## Cluster the hits (significant features)

Prepare arguments:

``` r
nr_clusters <- list(
    constant = 4,
    temp_shift = 4
)

plot_info <- list( # For the spline plots
    y_axis_label = "log2 counts",
    time_unit = "hours",
    treatment_labels = list(
        temp_shift = "temp shift"
    ),
    treatment_timepoints = list(
        temp_shift = 146
    )
)

genes <- rownames(data)
genes <- sub("_\\d+$", "", genes) # remove the _1 part from the end

plot_options <- list(
    # When meta_replicate_column is not there, all datapoints are blue.
    meta_replicate_column = "Reactor" # Colors the data points based on Reactor
)
```

Run the clustering:

``` r
clustering_results <- SplineOmics::cluster_hits(
    splineomics = splineomics,
    adj_pthresh_time_effect = 0.05,
    adj_pthresh_avrg_diff_conditions = 0.05,
    adj_pthresh_interaction_condition_time = 0.05,
    nr_clusters = nr_clusters,
    genes = genes
)
```

Generate the report:

``` r
report_dir <- file.path(tempdir(), "splineomics_report")
dir.create(report_dir, showWarnings = FALSE, recursive = TRUE)

plots <- SplineOmics::create_clustering_report(
    report_payload = clustering_results$report_payload,
    plot_info = plot_info,
    plot_options = plot_options,
    verbose = TRUE,
    max_hit_number = 5,
    report_dir = withr::local_tempdir()
)
#> Generating heatmap...
#> Generating cluster mean splines for level:  constant
#> Generating spline plots...
#> Generating cluster mean splines for level:  temp_shift
#> Generating spline plots...
#> Generating report. This takes a few seconds.
#> 
#>  Info Clustering the hits completed successfully.
#>  Your HTML reports are located in the directory:  /tmp/RtmpwxlMvB/filea6a95b4a24d1 .
#>  Please note that due to embedded files, the reports might be flagged as
#>  harmful by other software. Rest assured that they provide no harm.
```

You can view the generated analysis report of the cluster_hits function
[here](https://drive.google.com/uc?export=view&id=1EvPiEwHyGnBQFxRets00s87tochgyGK4).

As discussed before, there are three limma result categories. The
cluster_hits() report shows the results of all three, if they are
present (category 2 and 3 can only be generated when the design formula
contains an interaction effect).

## Perform gene set enrichment analysis (GSEA)

Once the clustered hits are identified, a subsequent step to gain
biological insights is to perform GSEA. For this, the respective genes
can be assigned to each clustered hit, and GSEA can be carried out. To
proceed, the Enrichr databases of choice need to be downloaded:

``` r
# Specify which databases you want to download from Enrichr
gene_set_lib <- c(
    "WikiPathways_2019_Human",
    "NCI-Nature_2016",
    "TRRUST_Transcription_Factors_2019",
    "MSigDB_Hallmark_2020",
    "GO_Cellular_Component_2018",
    "CORUM",
    "KEGG_2019_Human",
    "TRANSFAC_and_JASPAR_PWMs",
    "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
    "GO_Biological_Process_2018",
    "GO_Molecular_Function_2018",
    "Human_Gene_Atlas"
)

SplineOmics::download_enrichr_databases(
    gene_set_lib = gene_set_lib,
    output_dir = tempdir(), # output into the current working dir (default)
    filename = "databases.tsv" # just the name of the file, not the full path
)
```

Per default the file is placed in the current working directory, which
is the root dir of the R project.

To run GSEA, the downloaded database file has to be loaded as a
dataframe. Further, optionally, the clusterProfiler parameters and the
report dir can be specified. The function create_gsea_report() runs GSEA
using clusterProfiler, generates an HTML report and returns the GSEA
dotplots in R.

``` r
# Specify the filepath of the TSV file with the database info
downloaded_dbs_filepath <- file.path(
    tempdir(),
    "databases.tsv"
)

# Load the file
databases <- read.delim(
    downloaded_dbs_filepath,
    sep = "\t",
    stringsAsFactors = FALSE
)

# Specify the clusterProfiler parameters
clusterProfiler_params <- list(
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500
)

report_dir <- file.path(
    tempdir(),
    "results",
    "gsea_reports"
)
```

Run ORA:

``` r
ora_results <- SplineOmics::run_ora(
    cluster_table = clustering_results[["cluster_table"]],
    databases = databases,
    clusterProfiler_params = clusterProfiler_params,
    universe = genes,
    mapping_cfg = mapping_cfg,
    enrichGO_cfg = enrichGO_cfg
)
```

Create ORA HTML report:

``` r
plots <- SplineOmics::create_ora_report(
    report_payload = ora_results$report_payload,
    report_info = report_info,
    cluster_hits_report_name = "RNA_seq_analysis",
    verbose = TRUE,
    report_dir = results_dir
)
```

You can view the generated analysis report of the run_gsea function
[here](https://drive.google.com/uc?export=view&id=1EvPiEwHyGnBQFxRets00s87tochgyGK4).

## Session Info

    #> R version 4.5.2 (2025-10-31)
    #> Platform: x86_64-pc-linux-gnu
    #> Running under: Ubuntu 22.04.5 LTS
    #> 
    #> Matrix products: default
    #> BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0 
    #> LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0  LAPACK version 3.10.0
    #> 
    #> locale:
    #>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    #>  [3] LC_TIME=de_AT.UTF-8        LC_COLLATE=en_US.UTF-8    
    #>  [5] LC_MONETARY=de_AT.UTF-8    LC_MESSAGES=en_US.UTF-8   
    #>  [7] LC_PAPER=de_AT.UTF-8       LC_NAME=C                 
    #>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    #> [11] LC_MEASUREMENT=de_AT.UTF-8 LC_IDENTIFICATION=C       
    #> 
    #> time zone: Europe/Vienna
    #> tzcode source: system (glibc)
    #> 
    #> attached base packages:
    #> [1] stats     graphics  grDevices datasets  utils     methods   base     
    #> 
    #> other attached packages:
    #> [1] knitr_1.51        dplyr_1.2.0       readr_2.1.6       SplineOmics_0.4.3
    #> 
    #> loaded via a namespace (and not attached):
    #>   [1] RColorBrewer_1.1-3       rstudioapi_0.18.0        jsonlite_2.0.0          
    #>   [4] shape_1.4.6.1            magrittr_2.0.4           farver_2.1.2            
    #>   [7] nloptr_2.2.1             rmarkdown_2.30           GlobalOptions_0.1.3     
    #>  [10] fs_1.6.6                 ragg_1.5.0               vctrs_0.7.1             
    #>  [13] minqa_1.2.8              base64enc_0.1-6          htmltools_0.5.9         
    #>  [16] progress_1.2.3           broom_1.0.12             Formula_1.2-5           
    #>  [19] variancePartition_1.38.1 sass_0.4.10              KernSmooth_2.23-26      
    #>  [22] bslib_0.10.0             htmlwidgets_1.6.4        desc_1.4.3              
    #>  [25] pbkrtest_0.5.5           plyr_1.8.9               cachem_1.1.0            
    #>  [28] lifecycle_1.0.5          iterators_1.0.14         pkgconfig_2.0.3         
    #>  [31] Matrix_1.7-4             R6_2.6.1                 fastmap_1.2.0           
    #>  [34] rbibutils_2.4.1          clue_0.3-66              digest_0.6.39           
    #>  [37] numDeriv_2016.8-1.1      colorspace_2.1-2         S4Vectors_0.46.0        
    #>  [40] textshaping_1.0.4        labeling_0.4.3           abind_1.4-8             
    #>  [43] compiler_4.5.2           withr_3.0.2              bit64_4.6.0-1           
    #>  [46] aod_1.3.3                doParallel_1.0.17        S7_0.2.1                
    #>  [49] backports_1.5.0          BiocParallel_1.42.2      carData_3.0-6           
    #>  [52] gplots_3.3.0             MASS_7.3-65              rjson_0.2.23            
    #>  [55] corpcor_1.6.10           gtools_3.9.5             caTools_1.18.3          
    #>  [58] tools_4.5.2              otel_0.2.0               zip_2.3.3               
    #>  [61] remaCor_0.0.20           glue_1.8.0               nlme_3.1-168            
    #>  [64] grid_4.5.2               checkmate_2.3.4          cluster_2.1.8.1         
    #>  [67] reshape2_1.4.5           generics_0.1.4           gtable_0.3.6            
    #>  [70] tzdb_0.5.0               tidyr_1.3.2              hms_1.1.4               
    #>  [73] car_3.1-5                BiocGenerics_0.54.1      ggrepel_0.9.6           
    #>  [76] foreach_1.5.2            pillar_1.11.1            stringr_1.6.0           
    #>  [79] vroom_1.7.0              limma_3.64.3             circlize_0.4.17         
    #>  [82] splines_4.5.2            lattice_0.22-5           renv_1.1.7              
    #>  [85] gmp_0.7-5                bit_4.6.0                tidyselect_1.2.1        
    #>  [88] ComplexHeatmap_2.24.1    locfit_1.5-9.12          pbapply_1.7-4           
    #>  [91] reformulas_0.4.4         IRanges_2.42.0           edgeR_4.6.3             
    #>  [94] svglite_2.2.2            RhpcBLASctl_0.23-42      stats4_4.5.2            
    #>  [97] xfun_0.56                Biobase_2.68.0           statmod_1.5.1           
    #> [100] matrixStats_1.5.0        stringi_1.8.7            yaml_2.3.12             
    #> [103] boot_1.3-32              evaluate_1.0.5           codetools_0.2-19        
    #> [106] tibble_3.3.1             BiocManager_1.30.27      cli_3.6.5               
    #> [109] systemfonts_1.3.1        Rdpack_2.6.5             jquerylib_0.1.4         
    #> [112] Rcpp_1.1.1               EnvStats_3.1.0           png_0.1-8               
    #> [115] parallel_4.5.2           pkgdown_2.2.0            ggplot2_4.0.2           
    #> [118] prettyunits_1.2.0        ClusterR_1.3.6           bitops_1.0-9            
    #> [121] lme4_1.1-38              mvtnorm_1.3-3            lmerTest_3.2-0          
    #> [124] scales_1.4.0             purrr_1.2.1              crayon_1.5.3            
    #> [127] writexl_1.5.4            fANCOVA_0.6-1            GetoptLong_1.1.0        
    #> [130] rlang_1.1.7
