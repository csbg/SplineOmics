# Get started

## About this tutorial

This tutorial intends to showcase and explain the capabilities of the
**SplineOmics** package by walking through a real and complete example,
from start to finish.

#### Example Overview

The example involves a **time-series proteomics experiment**, where CHO
(chinese hamster ovary) cells were cultivated in three bioreactors
(three biological replicates). The experiment includes the following
setup:

- Samples were taken both during the **exponential** and **stationary
  growth phases**.
- Samples were collected in triplicates from each reactor at defined
  timepoints relative to cell feeding:
  - 60 minutes before feeding
  - 15, 60, 90, 120, and 240 minutes after feeding

#### Analysis Goals

The main goals of this analysis are:

- **Identify proteins with significant temporal changes**: Out of 7162
  cellular proteins, the objective is to detect which proteins show a
  significant change over time after the CHO cells were fed (i.e., the
  impact of the feeding).
- **Cluster hits based on temporal patterns**: The proteins (hits) with
  significant temporal changes will be clustered according to their
  time-based patterns.
- **Perform gene set enrichment analysis**: For each cluster, a gene set
  enrichment analysis will be performed to determine if specific
  biological processes are up- or downregulated after feeding.

#### Note

The documentation of all the **SplineOmics** package functions can be
viewed [here](https://csbg.github.io/SplineOmics/reference)

## Load the packages

``` r
library(SplineOmics)
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
library(org.Mm.eg.db) # BioConductor database
#> Loading required package: AnnotationDbi
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: 'generics'
#> The following object is masked from 'package:dplyr':
#> 
#>     explain
#> The following objects are masked from 'package:base':
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: 'BiocGenerics'
#> The following object is masked from 'package:dplyr':
#> 
#>     combine
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from 'package:base':
#> 
#>     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
#>     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
#>     get, grep, grepl, is.unsorted, lapply, Map, mapply, match, mget,
#>     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
#>     rbind, Reduce, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> Loading required package: IRanges
#> Loading required package: S4Vectors
#> 
#> Attaching package: 'S4Vectors'
#> The following objects are masked from 'package:dplyr':
#> 
#>     first, rename
#> The following object is masked from 'package:utils':
#> 
#>     findMatches
#> The following objects are masked from 'package:base':
#> 
#>     expand.grid, I, unname
#> 
#> Attaching package: 'IRanges'
#> The following objects are masked from 'package:dplyr':
#> 
#>     collapse, desc, slice
#> 
#> Attaching package: 'AnnotationDbi'
#> The following object is masked from 'package:dplyr':
#> 
#>     select
#> 
```

## Load the files

In this example, the proteomics_data.rds file contains the numeric
values (the intensities) and also the feature descriptions, such as gene
and protein name (= annotation part). Usually, you would load the data
from for example an Excel file, but the .rds file is more compressed,
which is the reason this format was chosen here to limit the size of the
SplineOmics package.

The file meta.xlsx contains the meta information, which are the
descriptions of the columns of the numeric values of data.

(These example files are part of the package and don’t have to be
present on your system).

Please note that this dataset is an actual experimental dataset, but the
annotation information, such as gene names, has been removed since it
was not yet published at the time of making the SplineOmics package
public. Instead, the dataset includes randomly generated gene symbols
and gene names corresponding to Cricetulus griseus (Chinese Hamster) for
each row. This is intended to demonstrate the functionality of the
package. Further, the dataset was subsamples to 1165 proteins to limit
the size of this package.

The left part of data contains the numeric values, and the right part
the annotation info, which can be copied in a separate dataframe, as
shown below.

``` r
data <- readRDS(xzfile(system.file(
    "extdata",
    "proteomics_data.rds.xz",
    package = "SplineOmics"
)))

meta <- read.csv(
    system.file(
        "extdata",
        "proteomics_meta.csv",
        package = "SplineOmics"
    ),
    stringsAsFactors = FALSE
)

# Extract the annotation part from the dataframe.
first_na_col <- which(is.na(data[1, ]))[1]
annotation <- data |>
    dplyr::select((first_na_col + 1):ncol(data)) |>
    dplyr::slice(-c(1:3))
```

### Show top rows of data

``` r
kable(
    head(data),
    format = "markdown"
)
```

| Sample ID | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 | 18 | 19 | 20 | 21 | 22 | 23 | 24 | 25 | 26 | 27 | 28 | 29 | 30 | 31 | 32 | 33 | 34 | 35 | 36 | …38 | Gene_symbol | Gene_name |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| Reactor | E09 | E10 | E12 | E09 | E10 | E12 | E09 | E10 | E12 | E09 | E10 | E12 | E09 | E10 | E12 | E09 | E10 | E12 | E09 | E10 | E12 | E09 | E10 | E12 | E09 | E10 | E12 | E09 | E10 | E12 | E09 | E10 | E12 | E09 | E10 | E12 | NA | NA | NA |
| Time Point | TP01 | TP01 | TP01 | TP02 | TP02 | TP02 | TP03 | TP03 | TP03 | TP04 | TP04 | TP04 | TP05 | TP05 | TP05 | TP06 | TP06 | TP06 | TP07 | TP07 | TP07 | TP08 | TP08 | TP08 | TP09 | TP09 | TP09 | TP10 | TP10 | TP10 | TP11 | TP11 | TP11 | TP12 | TP12 | TP12 | NA | NA | NA |
| Phase of Fermentation | Exponential | Exponential | Exponential | Exponential | Exponential | Exponential | Exponential | Exponential | Exponential | Exponential | Exponential | Exponential | Exponential | Exponential | Exponential | Exponential | Exponential | Exponential | Stationary | Stationary | Stationary | Stationary | Stationary | Stationary | Stationary | Stationary | Stationary | Stationary | Stationary | Stationary | Stationary | Stationary | Stationary | Stationary | Stationary | Stationary | NA | NA | NA |
| NA | 15.2908248901367 | 15.2148580551147 | 15.3249082565308 | 15.1468906402588 | 15.2205858230591 | 15.0987968444824 | 15.20960521698 | 15.2530574798584 | 15.2426490783691 | 15.1146211624146 | 15.2669916152954 | 15.3227806091309 | 15.2889413833618 | 15.2406797409058 | 15.3296279907227 | 15.1914215087891 | 15.2090301513672 | 15.2645998001099 | 16.4583435058594 | 16.5065593719482 | 16.7934131622314 | 16.60471534729 | 16.516918182373 | 16.7601490020752 | 16.7143936157227 | 16.4811534881592 | 16.9271640777588 | 16.6487922668457 | 16.5887832641602 | 16.817533493042 | 16.7007789611816 | 16.4686756134033 | 16.8468856811523 | 16.7374858856201 | 16.4771633148193 | 16.8459930419922 | NA | LOC113838844 | cone-rod homeobox protein-like |
| NA | 15.0816164016724 | 15.1488542556763 | 15.2205448150635 | 15.1780500411987 | 15.1131563186646 | 15.2813320159912 | 15.2288188934326 | 15.3091497421265 | 15.2797374725342 | 15.1128349304199 | 15.1366987228394 | 15.1883420944214 | 15.2453184127808 | 15.2480907440186 | 15.1772556304932 | 15.0709342956543 | 14.9833717346191 | 15.1247396469116 | 15.523494720459 | 15.4470376968384 | 15.3960628509521 | 15.4655313491821 | 15.4535474777222 | 15.4290637969971 | 15.2166357040405 | 15.3904008865356 | 15.2344093322754 | 15.4606895446777 | 15.4035310745239 | 15.338812828064 | 15.3438587188721 | 15.4040069580078 | 15.3590221405029 | 15.2574281692505 | 15.4244060516357 | 15.3422203063965 | NA | Wdr83os | WD repeat domain 83 opposite strand |
| NA | 14.5862588882446 | 14.7493143081665 | 14.6296472549438 | 14.5930347442627 | 14.6145257949829 | 14.6926259994507 | 14.5930233001709 | 14.6085433959961 | 14.7053565979004 | 14.5293817520142 | 14.6346235275269 | 14.7601327896118 | 14.5953159332275 | 14.6270151138306 | 14.7492513656616 | 14.6475257873535 | 14.6990842819214 | 14.6165771484375 | 14.9039011001587 | 14.8870449066162 | 15.0393772125244 | 14.8328456878662 | 14.9858465194702 | 15.0382928848267 | 14.8647422790527 | 14.7847366333008 | 15.1070947647095 | 14.8750152587891 | 15.0755500793457 | 15.021595954895 | 14.9314804077148 | 14.8654098510742 | 15.1668844223022 | 14.8610773086548 | 14.909330368042 | 14.9585390090942 | NA | Cubn | cubilin |

### Show top rows of meta

``` r
kable(
    head(meta),
    format = "markdown"
)
```

| Sample.ID            | Reactor | Time.Point | Phase       | Time |
|:---------------------|:--------|:-----------|:------------|-----:|
| E09_TP01_Exponential | E09     | TP01       | Exponential |  -60 |
| E10_TP01_Exponential | E10     | TP01       | Exponential |  -60 |
| E12_TP01_Exponential | E12     | TP01       | Exponential |  -60 |
| E09_TP02_Exponential | E09     | TP02       | Exponential |   15 |
| E10_TP02_Exponential | E10     | TP02       | Exponential |   15 |
| E12_TP02_Exponential | E12     | TP02       | Exponential |   15 |

### Show top rows of annotation

``` r
kable(
    head(annotation),
    format = "markdown"
)
```

| Gene_symbol  | Gene_name                                               |
|:-------------|:--------------------------------------------------------|
| LOC113838844 | cone-rod homeobox protein-like                          |
| Wdr83os      | WD repeat domain 83 opposite strand                     |
| Cubn         | cubilin                                                 |
| Dynlt1       | dynein light chain Tctex-type 1                         |
| Ostc         | oligosaccharyltransferase complex non-catalytic subunit |
| Unc5cl       | unc-5 family C-terminal like                            |

Three comments about the characteristics the input data should have:

- The data must not contain any NA values or other special values, and
  must consist only of numbers. For example, the original proteomics
  data contained some NA values, which were resolved in this case by
  imputation (replacing NA values with numbers).

- All features of the data should ideally be normally distributed when
  analyzed with `limma`, which fits a linear model to each feature.
  These models rely on statistical tests that assume normality. Although
  `limma` can still function if the data is not normally distributed,
  the resulting p-values may become less reliable. For this reason, it
  is strongly recommended to transform your data using techniques such
  as log2 transformation when the features deviate from normality.
  Proper transformation helps ensure that the assumptions underlying the
  statistical tests are met, leading to more accurate and trustworthy
  results.

- The samples in the data should be independent of each other. Linear
  models,  
  such as those used in limma, assume that the observations (samples)
  are  
  independent. If there is a dependency between samples (e.g.,
  repeated  
  measurements of the same subject), this assumption is violated, which
  can lead to incorrect statistical inferences.

### Bring the Inputs into the Standardized Format

Since `data` is not in the format required by the **SplineOmics**
package, it needs some processing. The SplineOmics package requires data
to be a numeric matrix, so no element is allowed to be anything else
than a number. This can be done with a few commands in R, but if your
file has a specific structure, the function
[`extract_data()`](https://csbg.github.io/SplineOmics/reference/extract_data.md)
can handle this more easily

#### Usage of the extract_data() function

[`extract_data()`](https://csbg.github.io/SplineOmics/reference/extract_data.md)
can:

- **Extract the data matrix field** by specifying the location of the
  corners of the matrix.
- **Create column headers** from the information written in the cells
  above the respective columns of the data matrix field.
- **Assign rowheaders**:
  - If no annotation columns are specified, rowheaders will be
    increasing numbers.
  - If annotation columns are specified (like
    `"First.Protein.Description"` and `"ID"` in this example), these
    will be combined to form the rowheaders (feature names).

#### Usage in Plotting

The generated rowheaders will be used to label any plots where a feature
is shown individually, such as:

- **Spline plots** with the datapoints from an individual feature.

``` r
data <- SplineOmics::extract_data(
    # The dataframe with the numbers on the left and info on the right.
    data = data,
    # Use this annotation column for the feature names.
    feature_name_columns = c("Gene_name"),
    use_row_index = TRUE,
    top_row = 4,
    bottom_row = 1165,
    right_col = 37,
    left_col = 2
)
```

## Perform EDA (exploratory data analysis)

Now that we have the data in the required format (numeric matrix) we can
go on.

The first step in analyzing data is typically **Exploratory Data
Analysis (EDA)**. EDA involves summarizing the main characteristics of
the data, often through visualizations.

#### Common EDA Plots

Some common types of EDA plots include:

- **Density distributions**
- **Boxplots**
- **PCA (Principal Component Analysis)**
- **Correlation heatmaps**

Again, you can generate those plots yourself with a few lines of R code.
However, if you prefer, for convenience, the
[`explore_data()`](https://csbg.github.io/SplineOmics/reference/explore_data.md)
function can handle this for you.

#### Using `explore_data()` for EDA

The **SplineOmics** package provides the function
[`explore_data()`](https://csbg.github.io/SplineOmics/reference/explore_data.md)
to perform EDA. This function requires the following arguments:

- **data**: The numeric data matrix.
- **meta**: The metadata table.
- **condition**: The name of the column in the metadata that contains
  the levels of the experiment (e.g., “Exponential” and “Stationary”).
- **report_info**: A list that contains general information about the
  analysis, such as the name of the analyst and the datatype (e.g.
  proteomics)

#### Optional Arguments

In addition to the required arguments,
[`explore_data()`](https://csbg.github.io/SplineOmics/reference/explore_data.md)
offers several optional arguments:

- **meta_batch_column**: The name of the column that contains the first
  batch effect.

- **meta_batch2_column**: The name of the column that contains the
  second batch effect.

  If at least one batch column is provided, the function will:

  - Use the `removeBatchEffect()` function from **limma** to remove the
    batch effect from the data before plotting.
  - Generate two EDA HTML reports: one for the **uncorrected data** and
    one for the **batch-corrected data**.

#### Output and Report Options

- By default, the reports are saved in the **current working
  directory**, but this location can be changed using the `report_dir`
  argument.
- The function also **returns all plots** generated during the analysis,
  so that you can modify them according to your own needs.
- If you do not want a report to be generated, you can set the `report`
  argument to `FALSE` (when you for example just want the figures in the
  R environment)

``` r
# Those fields are mandatory, because we believe that when such a report is
# opened after half a year, those infos can be very helpful.
report_info <- list(
    omics_data_type = "PTX",
    data_description = "Proteomics data of CHO cells",
    data_collection_date = "February 2024",
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

### SplineOmics Object

In the SplineOmics package, multiple functions take the same arguments
as input. To make this easier and to avoid errors, we decided that those
arguments are not provided individually to the functions, but are all
stored in an R6 object (which is of type ‘SplineOmics’) and then this
object is passed to the functions. Additionally, some functions generate
intermediate output, which is just necessary for the next function in
the workflow, which is then also just passed along by updating the
SplineOmics object. But you don’t have to worry about this.

#### Functionality

The SplineOmics object can be seen as a container where all necessary
arguments are stored. Each function retrieves the required arguments
from the object and potentially adds new data or results back into it.

#### Documentation

The documentation of the function that creates the SplineOmics object
can be found
[here](https://csbg.github.io/SplineOmics/reference/create_splineomics.html)
and the documentation of the function that updates it
\[[here](https://csbg.github.io/SplineOmics/reference/update_splineomics.html)

The documentation for each function that takes the SplineOmics object as
input specifies which arguments must be present in the SplineOmics
object when it is passed to the respective function.

### Required Arguments `create_splineomics()`

- **data**: A matrix with the data
- **meta**: Metadata associated with the data.
- **condition**: Meta column name of the levels (e.g., Exponential and
  Stationary). In this analysis we illustrate a case with two
  conditions, but SplineOmics supports any number of conditions. All
  conditions are modeled jointly, and SplineOmics automatically performs
  all pairwise comparisons. For each pair of conditions, it computes
  both the average difference (limma Category 2) and the condition–time
  interaction (limma Category 3), providing a complete set of contrast
  results regardless of how many conditions are included.

### Optional Arguments `create_splineomics()`

- **rna_seq_data**: An object containing the preprocessed RNA-seq data,
  such as the output from
  [`limma::voom`](https://rdrr.io/pkg/limma/man/voom.html) function.
- **annotation**: A dataframe with the feature descriptions of data.
- **report_info**: A list containing general information about the
  analysis.
- **meta_batch_column**: Column for meta batch information.
- **meta_batch2_column**: Column for secondary meta batch information.
- **design**: A limma design formula
- **spline_params**: Parameters for the spline functions.

``` r
# splineomics now contains the SplineOmics object.
splineomics <- SplineOmics::create_splineomics(
    data = data,
    meta = meta,
    annotation = annotation,
    report_info = report_info,
    condition = "Phase", # Column of meta that contains the levels.
    meta_batch_column = "Reactor" # For batch effect removal
)

# Special print.SplineOmics function leads to selective printing
print(splineomics)
#> data:SplineOmics Object
#> -------------------
#> Number of features (rows): 1162 
#> Number of samples (columns): 36 
#> Meta data columns: 5 
#> First few meta columns:
#>              Sample.ID Reactor Time.Point       Phase Time
#> 1 E09_TP01_Exponential     E09       TP01 Exponential  -60
#> 2 E10_TP01_Exponential     E10       TP01 Exponential  -60
#> 3 E12_TP01_Exponential     E12       TP01 Exponential  -60
#> Condition: Phase 
#> No RNA-seq data provided.
#> Annotation provided with 1162 entries.
#> No spline parameters set.
#> P-value adjustment method: BH
```

Now that we have the SplineOmics object defined, we can perform our
exploratory data analysis.

``` r
plots <- SplineOmics::explore_data(
    splineomics = splineomics, # SplineOmics object
    report_dir = withr::local_tempdir()
)
#> Making density plots...
#> Making violin plots...
#> Making PCA plots...
#> Making MDS plots...
#> Making correlation heatmaps...
#> Subsampled to top 1000 most variable features (after filtering rows with > 20% missing) for correlation heatmap.
#> Making mean correlation with time plots...
#> Making lag1 differences plots...
#> Making first lag auto-correlation with time plots...
#> Making cv plots...
#> Making density plots...
#> Making violin plots...
#> Making PCA plots...
#> Making MDS plots...
#> Making correlation heatmaps...
#> Subsampled to top 1000 most variable features (after filtering rows with > 20% missing) for correlation heatmap.
#> Making mean correlation with time plots...
#> Making lag1 differences plots...
#> Making first lag auto-correlation with time plots...
#> Making cv plots...
#> 
#>  Info Exploratory data analysis completed successfully.
#>  Your HTML reports are located in the directory:  /tmp/RtmpPca44v/filea4f8217c711a .
#>  Please note that due to embedded files, the reports might be flagged as
#>  harmful by other software. Rest assured that they provide no harm.
```

[Here](https://csbg.github.io/SplineOmics_html_reports/explore_data_PTX.html)
you can see the HTML report of the explore_data() function with the NOT
batch-corrected data, and
[here](https://csbg.github.io/SplineOmics_html_reports/explore_batch_corrected_data_PTX.html)
the report for the batch-corrected data.

Note that any report linked in this vignette might have been created
with an older version of SplineOmics, and not updated yet.

The EDA plots can tell you a range of things. The plots in the HTML
report are grouped into three categories: Distribution and Variability
Analysis, Time Series Analysis, and Dimensionality Reduction and
Clustering.

If you look at the correlation heatmaps in the HTML report, you can see
that the samples E12_TP05_Exponential and E10_TP10_Stationary stick out.
Seeing this, you might want to remove them from the data.

## Run limma spline analysis

Now we decide on the parameters and run the limma spline analysis.

For the design formula, you must specify either ‘isolated’ or
‘integrated’. Isolated means limma determines the results for each level
using only the data from that level. Integrated means limma determines
the results for all levels using the full dataset (from all levels). For
the integrated mode, the condition column (here Phase) must be included
in the design. Isolated means that limma only uses the part of the
dataset that belongs to a level to obtain the results for that level.

To generate the limma result categories 2 and 3 ()

``` r
splineomics <- SplineOmics::update_splineomics(
    splineomics = splineomics,
    design = "~ 1 + Phase*Time + Reactor", # best design formula
    mode = "integrated", # means limma uses the full data for each condition.
    # States explicitly that there is no problem of heteroscedasticity and
    # therefore, this does not need to be adressed. Setting it to TRUE would
    # mean
    # the opposite, and when setting it to NULL, it means it should be handled
    # implicitly. For details, see Reference
    # documentation of the create_splineomics() function.
    use_array_weights = FALSE,
    spline_params = list(
        spline_type = c("n"), # natural cubic splines (take these if unsure)
        # If you are unsure about which dof, start with 2 and increase
        dof = c(2L) 
    )
)
```

Run the
[`run_limma_splines()`](https://csbg.github.io/SplineOmics/reference/run_limma_splines.md)
function with the updated SplineOmics object:

``` r
splineomics <- SplineOmics::run_limma_splines(
    splineomics = splineomics
)
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
[here](https://csbg.github.io/SplineOmics_html_reports/create_limma_report_PTX.html).

This report contains p-value histograms for all three limma result
categories and a volcano plot for category 2. Embedded in this file are
the downloadable limma topTables with the results for category 1 if the
mode was ‘isolated’ and also with the results of category 2 and 3 if the
mode was ‘integrated’.

Note that in the upcoming cluster_hits() function report, the embedded
file will only contain the clustered significant features from result
category 1.

## Cluster the hits (significant features)

After we obtained the limma spline results, we can cluster the hits
based on their temporal pattern (their spline shape). We define what a
hit is by setting an adj. p-value threshold for every level. Hits are
features (e.g. proteins) that have an adj. p-value below the threshold.
K-means clustering is used to place every hit in one of as many clusters
as we have specified for that specific level.

``` r
# The amount of clusters can be a fixed number (e.g. 6) or a range. When you
# specify a range (e.g. 2:3, which corresponds to 2 3 in the vector) then the
# cluster_hits() function tries all those cluster numbers and picks the one with
# the highest silhouette score (automatic cluster number identification). When
# you don't want to have a clustering for a level, write 1 for the cluster
# number for that level.
nr_clusters <- list(
    Exponential = 6, # specifically 6 clusters for the exponential phase level
    Stationary = 2:3 # range of cluster numbers for the stationary phase level
)

plot_info <- list( # For the spline plots
    y_axis_label = "log2 intensity",
    time_unit = "min", # our measurements were in minutes
    treatment_labels = list(
        Exponential = "feeding",
        Stationary = "feeding"
    ),
    treatment_timepoints = list(
        Exponential = 0,
        Stationary = 0
    )
)


# Get all the gene names. They are used for generating files
# which contents can be directly used as the input for the Enrichr webtool,
# if you prefer to manually perform the enrichment. Those files are
# embedded in the output HTML report and can be downloaded from there.
gene_column_name <- "Gene_symbol"
genes <- annotation[[gene_column_name]]

plot_options <- list(
    # When meta_replicate_column is not there, all datapoints are blue.
    meta_replicate_column = "Reactor", # Colors the data points based on Reactor
    cluster_heatmap_columns = FALSE # Per default FALSE, just for demonstration
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
plots <- SplineOmics::create_clustering_report(
    report_payload = clustering_results$report_payload,
    plot_info = plot_info,
    plot_options = plot_options,
    verbose = TRUE,
    max_hit_number = 5,
    report_dir = withr::local_tempdir()
)
#> Generating heatmap...
#> Generating cluster mean splines for level:  Exponential
#> Generating spline plots...
#> Generating cluster mean splines for level:  Stationary
#> Generating spline plots...
#> Generating report. This takes a few seconds.
#> 
#>  Info Clustering the hits completed successfully.
#>  Your HTML reports are located in the directory:  /tmp/RtmpPca44v/filea4f81b2d6c9c .
#>  Please note that due to embedded files, the reports might be flagged as
#>  harmful by other software. Rest assured that they provide no harm.
```

You can view the generated analysis report of the clustering
[here](https://csbg.github.io/SplineOmics_html_reports/report_clustered_hits_PTX.html).

As discussed before, there are three limma result categories. The
cluster_hits() report shows the results of all three, if they are
present (category 2 and 3 can only be generated when the design formula
contains an interaction effect).

## Perform overrepresentation analysis (ORA)

Once the clustered hits are identified, a subsequent step to gain
biological insights is to perform ORA For this, the respective genes can
be assigned to each clustered hit, and ORA can be carried out. To
proceed, the Enrichr databases of choice need to be downloaded:

``` r
# Create a temporary directory for R CMD check
results_dir <- file.path(tempdir(), "ora")
dir.create(report_dir, showWarnings = FALSE, recursive = TRUE)

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
    output_dir = results_dir,
    filename = "databases.tsv"
)
```

Per default the file is placed in the current working directory, which
is the root dir of the R project.

To run ORA, the downloaded database file has to be loaded as a
dataframe. Further, optionally, the clusterProfiler parameters and the
report dir can be specified. The function run_ora() runs ORA using
clusterProfiler, generates an HTML report and returns the ORA dotplots
in R.

``` r
# Specify the filepath of the TSV file with the database info
downloaded_dbs_filepath <- file.path(results_dir, "databases.tsv")

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
    maxGSSize = 500,
    qvalueCutoff = 0.2
)
```

Run ORA:

``` r
ora_results <- SplineOmics::run_ora(
    cluster_table = clustering_results[["cluster_table"]],
    databases = databases,
    clusterProfiler_params = clusterProfiler_params,
    universe = genes
)
```

Create ORA HTML report:

``` r
plots <- SplineOmics::create_ora_report(
    report_payload = ora_results$report_payload,
    report_info = report_info,
    cluster_hits_report_name = "get_started",
    verbose = TRUE,
    report_dir = results_dir
)
```

You can view the generated analysis report of the run_ora function
[here](https://csbg.github.io/SplineOmics_html_reports/run_ora_report_PTX.html).

This report first shows all enrichment results, where more than 2 genes
supported a term, in a tabular format. The table with all the terms with
\< 2 genes supporting it can be downloaded by clicking on a button below
that table.

For the dotplots below that, every row is a term from a specific
database, and the columns are the respective clusters. The color scale
contains the info about the odds ratio and the size the -log10 adj.
p-value. Only terms that have \> 2 genes as support are included in the
plot. Further, for each cluster, just maximally 5 terms are shown (the
terms with the highest odds ratios). Note that when for example cluster
1 already has 5 terms, and cluster 2 does not, and gets a term which was
also found for cluster 1, than this term would be included as the sixth
term for cluster 1, so this is a way the maximum of 5 can be exceeded.

If a phase, like stationary here, does not lead to any enrichment
results, that is stated with a red message.

## Perform ORA with the Bioconductor database

The ORA can also be performed with any other database that is in the
format DB, Geneset, Gene (see documentation of the run_ora function).
The BioConductor databases are an example for that, and SplineOmics also
contains a function to conveniently download these.

``` r
SplineOmics::extract_gene_sets(
    organism_db = "org.Mm.eg.db",
    output_dir = results_dir,
    filename = "bioconductor_database.tsv"
)
```

``` r
# Specify the filepath of the TSV file with the database info
downloaded_dbs_filepath <- file.path(
    results_dir,
    "bioconductor_database.tsv"
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

mapping_cfg <- list(
    method = "gprofiler",
    from_species = "mmusculus",
    to_species = "cgchok1gshd"  
)

genes_clean <- genes[!grepl("-", genes)]

enrichGO_cfg <- list(
    GO_BP = list(
        OrgDb = org.Mm.eg.db,
        keyType = "SYMBOL",
        ontology = "BP"
    ),
    GO_MF = list(
        OrgDb = org.Mm.eg.db,
        keyType = "SYMBOL",
        ontology = "MF"
    ),
    GO_CC = list(
        OrgDb = org.Mm.eg.db,
        keyType = "SYMBOL",
        ontology = "CC"
    )
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
    cluster_hits_report_name = "get_started",
    verbose = TRUE,
    report_dir = results_dir
)
```

## Conclusion

This example showed most functionalities of the SplineOmics package. You
can also run other datatypes with it, including timeseries RNA-seq and
glycan data (for those, refer to the documentation in the README file on
the GitHub page under Usage/RNA-seq and Glycan Data).

We hope that the SplineOmics package makes your scientific data analysis
easier. If you face any problems (bugs in the code) or are not satisfied
with the documentation, open an issue on GitHub or check out the other
options under the Feedback section of the README on GitHub. Thank you!

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
    #> [1] stats4    stats     graphics  grDevices datasets  utils     methods  
    #> [8] base     
    #> 
    #> other attached packages:
    #>  [1] org.Mm.eg.db_3.21.0  AnnotationDbi_1.70.0 IRanges_2.42.0      
    #>  [4] S4Vectors_0.46.0     Biobase_2.68.0       BiocGenerics_0.54.1 
    #>  [7] generics_0.1.4       knitr_1.51           dplyr_1.2.0         
    #> [10] SplineOmics_0.4.3   
    #> 
    #> loaded via a namespace (and not attached):
    #>   [1] RColorBrewer_1.1-3       rstudioapi_0.18.0        jsonlite_2.0.0          
    #>   [4] shape_1.4.6.1            magrittr_2.0.4           farver_2.1.2            
    #>   [7] nloptr_2.2.1             rmarkdown_2.30           GlobalOptions_0.1.3     
    #>  [10] fs_1.6.6                 ragg_1.5.0               vctrs_0.7.1             
    #>  [13] memoise_2.0.1            minqa_1.2.8              base64enc_0.1-6         
    #>  [16] htmltools_0.5.9          progress_1.2.3           broom_1.0.12            
    #>  [19] Formula_1.2-5            variancePartition_1.38.1 sass_0.4.10             
    #>  [22] KernSmooth_2.23-26       bslib_0.10.0             htmlwidgets_1.6.4       
    #>  [25] desc_1.4.3               pbkrtest_0.5.5           plyr_1.8.9              
    #>  [28] cachem_1.1.0             lifecycle_1.0.5          iterators_1.0.14        
    #>  [31] pkgconfig_2.0.3          Matrix_1.7-4             R6_2.6.1                
    #>  [34] fastmap_1.2.0            GenomeInfoDbData_1.2.14  rbibutils_2.4.1         
    #>  [37] clue_0.3-66              digest_0.6.39            numDeriv_2016.8-1.1     
    #>  [40] colorspace_2.1-2         textshaping_1.0.4        RSQLite_2.4.5           
    #>  [43] labeling_0.4.3           httr_1.4.7               abind_1.4-8             
    #>  [46] compiler_4.5.2           withr_3.0.2              bit64_4.6.0-1           
    #>  [49] aod_1.3.3                doParallel_1.0.17        S7_0.2.1                
    #>  [52] backports_1.5.0          BiocParallel_1.42.2      carData_3.0-6           
    #>  [55] DBI_1.2.3                gplots_3.3.0             MASS_7.3-65             
    #>  [58] rjson_0.2.23             corpcor_1.6.10           gtools_3.9.5            
    #>  [61] caTools_1.18.3           tools_4.5.2              otel_0.2.0              
    #>  [64] zip_2.3.3                remaCor_0.0.20           glue_1.8.0              
    #>  [67] nlme_3.1-168             grid_4.5.2               checkmate_2.3.4         
    #>  [70] cluster_2.1.8.1          reshape2_1.4.5           gtable_0.3.6            
    #>  [73] tidyr_1.3.2              hms_1.1.4                XVector_0.48.0          
    #>  [76] car_3.1-5                ggrepel_0.9.6            foreach_1.5.2           
    #>  [79] pillar_1.11.1            stringr_1.6.0            limma_3.64.3            
    #>  [82] circlize_0.4.17          splines_4.5.2            lattice_0.22-5          
    #>  [85] renv_1.1.7               gmp_0.7-5                bit_4.6.0               
    #>  [88] tidyselect_1.2.1         ComplexHeatmap_2.24.1    Biostrings_2.76.0       
    #>  [91] pbapply_1.7-4            reformulas_0.4.4         svglite_2.2.2           
    #>  [94] RhpcBLASctl_0.23-42      xfun_0.56                statmod_1.5.1           
    #>  [97] matrixStats_1.5.0        UCSC.utils_1.4.0         stringi_1.8.7           
    #> [100] yaml_2.3.12              boot_1.3-32              evaluate_1.0.5          
    #> [103] codetools_0.2-19         tibble_3.3.1             BiocManager_1.30.27     
    #> [106] cli_3.6.5                systemfonts_1.3.1        Rdpack_2.6.5            
    #> [109] jquerylib_0.1.4          GenomeInfoDb_1.44.3      Rcpp_1.1.1              
    #> [112] EnvStats_3.1.0           png_0.1-8                parallel_4.5.2          
    #> [115] pkgdown_2.2.0            ggplot2_4.0.2            blob_1.3.0              
    #> [118] prettyunits_1.2.0        ClusterR_1.3.6           bitops_1.0-9            
    #> [121] lme4_1.1-38              mvtnorm_1.3-3            lmerTest_3.2-0          
    #> [124] scales_1.4.0             purrr_1.2.1              crayon_1.5.3            
    #> [127] writexl_1.5.4            fANCOVA_0.6-1            GetoptLong_1.1.0        
    #> [130] rlang_1.1.7              KEGGREST_1.48.1
