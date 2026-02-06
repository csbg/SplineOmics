# methylation-data-analysis

## About this vignette

This tutorial intends to showcase and explain the capabilities of the
**SplineOmics** package by walking through a real and complete
methylation data analysis example, from start to finish. **SplineOmics**
is explained in more detail in the **get-started** vignette, where a
proteomics example is covered. This vignette is more focused on showing
how methylation data can be used, and because of this, less details
about the overall package are provided here. Further, this vignette is
more minimalistic than the other ones, because the same functionalities
than in the other vignettes are used and the explanations are not
repeated here.

### Data Overview

This dataset originates from a **time-series methylation data
experiment** designed to study Chinese Hamster Ovary (CHO) cells. The
experiment involved cultivating cells in **eight bioreactors**, with
**four bioreactors** subjected to a temperature shift after 146 hours
(experimental condition) and the remaining **four bioreactors**
maintained without a temperature shift (control condition).

#### Timepoints

Samples were collected at **19 distinct time points** throughout the
experiment, specifically: `"48h"`, `"72h"`, `"96h"`, `"120h"`, `"144h"`,
`"168h"`, `"192h"`, `"216h"`, `"240h"`, `"264h"`, `"268h"`, `"288h"`,
`"312h"`, `"336h"`, `"360h"`, `"384h"`, `"408h"`, `"432h"`, `"456h"`,
after cultivation start. Most time point were sampled from all eight
bioreactors, resulting in a total of **130 samples**.

#### Additonal effects in the experiment: Reactor and Group

In this experiment, there are two effects to consider: **Reactor** and
**Group**.

1.  **Reactor**:
    - This refers to the different bioreactors used for cell
      cultivation, which can exhibit substantial variability.
    - Each reactor was assigned a single condition: either constant
      temperature or temperature-shifted. As a result, reactor is a
      blocked effect.
    - Reactor should not be treated as a fixed effect to simply remove
      its influence. Instead, it is treated as a **random effect**,
      which allows us to model its variability appropriately.
2.  **Group**:
    - This refers to the two separate groups in which the bioreactors
      where run. Because the apparatus containes 4 bioreactors, only
      that amount can be run simultaneosly. The other 4 were run after
      the initial 4.
    - Group is considered a **batch effect** with respect to the
      condition (constant temperature vs. temperature-shifted).

Since reactor is a blocked effect, the variability due to the reactors
cannot be directly separated from the condition. Instead, **linear mixed
models (LMMs)** are used to attribute the reactor as a random effect,
allowing us to account for its variability while isolating the effects
of the condition. This approach ensures that the analysis appropriately
handles the hierarchical structure of the data and avoids incorrect
conclusions.

#### Note

The documentation of all the **SplineOmics** package functions can be
viewed [here](https://csbg.github.io/SplineOmics/reference)

## Load the packages

``` r
library(SplineOmics)
library(readr) # For reading the meta CSV file
library(here) # For managing filepaths
#> here() starts at /home/thomas/Documents/PhD/projects/DGTX/SplineOmics_hub/SplineOmics
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
data_df <- readRDS(xzfile(system.file(
    "extdata",
    "methylation_data.rds.xz",
    package = "SplineOmics"
)))

meta <- readr::read_csv(
    system.file(
        "extdata",
        "methylation_meta.csv",
        package = "SplineOmics"
    ),
    show_col_types = FALSE
)
meta$Condition[meta$Condition == "Temp. Shifted"] <- "Temp_shift"
```

### Show top rows of data

``` r
kable(
    head(data),
    format = "markdown"
)
```

|                                                                       |
|:----------------------------------------------------------------------|
| function (…, list = character(), package = NULL, lib.loc = NULL,      |
| verbose = getOption(“verbose”), envir = .GlobalEnv, overwrite = TRUE) |
| {                                                                     |
| fileExt \<- function(x) {                                             |
| db \<- grepl(“\\\[^.\]+\\(gz\|bz2\|xz)\$”, x)                         |
| ans \<- sub(“.\*\\”, ““, x)                                           |

### Show top rows of meta

``` r
kable(
    head(meta),
    format = "markdown"
)
```

| No. | TP | Time | Reactor | Condition | Group | Cell Line | Type | Date Isolated | Tube Label | A260/280 | Tapestation Concentration | DIN |
|---:|---:|---:|:---|:---|:---|:---|:---|:---|:---|---:|---:|---:|
| 1 | 5 | 48 | E13 | Constant | CharRun1 | NISTCHO | 1x10^6 cells | 2024-11-04 | NIST CHO BR24 48H E13 | 1.888 | 45.4 | 9.9 |
| 3 | 5 | 48 | E15 | Constant | CharRun1 | NISTCHO | 1x10^6 cells | 2024-11-04 | NIST CHO BR24 48H E15 | 1.869 | 50.7 | 9.9 |
| 65 | 5 | 48 | E17 | Constant | CharRun2 | NISTCHO | 1x10^6 cells | 2024-11-13 | NIST CHO CHARRUN2 48H E17 | 1.881 | 68.0 | 9.8 |
| 67 | 5 | 48 | E19 | Constant | CharRun2 | NISTCHO | 1x10^6 cells | 2024-11-13 | NIST CHO CHARRUN2 48H E19 | 1.885 | 39.7 | 9.9 |
| 2 | 5 | 48 | E14 | Temp_shift | CharRun1 | NISTCHO | 1x10^6 cells | 2024-11-04 | NIST CHO BR24 48H E14 | 1.874 | 60.2 | 9.9 |
| 4 | 5 | 48 | E16 | Temp_shift | CharRun1 | NISTCHO | 1x10^6 cells | 2024-11-04 | NIST CHO BR24 48H E16 | 1.892 | 30.9 | 9.9 |

``` r
data <- SplineOmics::extract_data(
    # The dataframe with the numbers on the left and info on the right.
    data = data_df,
    # Use this annotation column for the feature names.
    feature_name_columns = c("feature_names"),
    top_row = 1,
    bottom_row = 300,
    right_col = 131,
    left_col = 2
)
```

## Perform EDA (exploratory data analysis)

``` r
report_info <- list(
    omics_data_type = "EPI",
    data_description = "Methylation data of CHO cells",
    data_collection_date = "December 2024",
    analyst_name = "Thomas Rauter",
    contact_info = "thomas.rauter@plus.ac.at",
    project_name = "DGTX"
)

report_dir <- here::here(
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
    meta_batch_column = "Group" # Remove batch effect for plotting.
)
```

``` r
plots <- SplineOmics::explore_data(
    splineomics = splineomics,
    report_dir = withr::local_tempdir()
)
#> Making density plots...
#> Making violin plots...
#> Making PCA plots...
#> 25 rows with missing values removed before PCA.
#> Making MDS plots...
#> 25 rows with missing values removed before MDS
#> Making correlation heatmaps...
#> Making mean correlation with time plots...
#> Making lag1 differences plots...
#> Making first lag auto-correlation with time plots...
#> 25 rows with missing values removed before firstlag autocorrelation plot
#> Making cv plots...
#> Warning: Removed 367 rows containing non-finite outside the scale range
#> (`stat_density()`).
#> Warning: Removed 177 rows containing non-finite outside the scale range
#> (`stat_density()`).
#> Warning: Removed 190 rows containing non-finite outside the scale range
#> (`stat_density()`).
#> Warning: Removed 177 rows containing non-finite outside the scale range
#> (`stat_ydensity()`).
#> Warning: Removed 177 rows containing non-finite outside the scale range
#> (`stat_boxplot()`).
#> Warning: Removed 190 rows containing non-finite outside the scale range
#> (`stat_ydensity()`).
#> Warning: Removed 190 rows containing non-finite outside the scale range
#> (`stat_boxplot()`).
#> Making density plots...
#> Making violin plots...
#> Making PCA plots...
#> 25 rows with missing values removed before PCA.
#> Making MDS plots...
#> 25 rows with missing values removed before MDS
#> Making correlation heatmaps...
#> Making mean correlation with time plots...
#> Making lag1 differences plots...
#> Making first lag auto-correlation with time plots...
#> 25 rows with missing values removed before firstlag autocorrelation plot
#> Making cv plots...
#> Warning: Removed 367 rows containing non-finite outside the scale range
#> (`stat_density()`).
#> Warning: Removed 177 rows containing non-finite outside the scale range
#> (`stat_density()`).
#> Warning: Removed 190 rows containing non-finite outside the scale range
#> (`stat_density()`).
#> Warning: Removed 177 rows containing non-finite outside the scale range
#> (`stat_ydensity()`).
#> Warning: Removed 177 rows containing non-finite outside the scale range
#> (`stat_boxplot()`).
#> Warning: Removed 190 rows containing non-finite outside the scale range
#> (`stat_ydensity()`).
#> Warning: Removed 190 rows containing non-finite outside the scale range
#> (`stat_boxplot()`).
#> 
#>  Info Exploratory data analysis completed successfully.
#>  Your HTML reports are located in the directory:  /tmp/RtmpNQV3EO/file45aa6c6dd68a .
#>  Please note that due to embedded files, the reports might be flagged as
#>  harmful by other software. Rest assured that they provide no harm.
```

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

design <- "~ 1 + Condition*Time + Group + (1|Reactor)"

splineomics <- SplineOmics::update_splineomics(
    splineomics = splineomics,
    data = data,
    design = design,
    mode = "integrated",
    spline_params = spline_params,
    bp_cfg = c( # For parallel computing for variancePartition::dream (lmm)
        n_cores = 1, # For the vignette, just use one core.
        blas_threads = 1
    )
)
```

Run the
[`run_limma_splines()`](https://csbg.github.io/SplineOmics/reference/run_limma_splines.md)
function:

``` r
splineomics <- SplineOmics::run_limma_splines(
    splineomics = splineomics
)
#> Hint: The data contains missing values (NA). Ensure that imputation or handling of missing values aligns with your analysis requirements. Note that limma can handle missing values (it just ignores them), and therefore SplineOmics can also handle them.
#> 
#> Ensure the design has no Condition*Time interaction for mode == 'isolated', and includes it for mode == 'integrated'.
#> 
#> Fitting global model...
#> 
#> NOTE: If you manually stop run_limma_splines() in RStudio and  used parallelization for variancePartition::dream(), then those parallelized processes may continue running. Use your system's process manager to terminate them manually!
#> Warning in variancePartition::dream(exprObj = data, formula = stats::as.formula(design), : Model failed for 8 responses.
#>   See errors with attr(., 'errors')
#> Info Finished limma spline analysis in 0.2 min
```

## Build limma report

The topTables of all three limma result categories can be used to
generate p-value histograms an volcano plots.

``` r
plots <- SplineOmics::create_limma_report(
    splineomics = splineomics,
    report_dir = withr::local_tempdir()
)
#> 
#>  Info Limma report generation completed successfully.
#>  Your HTML reports are located in the directory:  /tmp/RtmpNQV3EO/file45aa4b767150 .
#>  Please note that due to embedded files, the reports might be flagged as
#>  harmful by other software. Rest assured that they provide no harm.
```

## Cluster the hits (significant features)

Prepare arguments:

``` r
nr_clusters <- list(
    Constant   = 2:6,
    Temp_shift = 2:6
)

report_dir <- here::here(
    "results",
    "clustering_reports"
)

plot_info <- list( # For the spline plots
    y_axis_label = "beta value",
    time_unit = "hours",
    treatment_labels = list(
        Temp_shift = "temp shift"
    ),
    treatment_timepoints = list(
        Temp_shift = 146
    )
)

genes <- rownames(data)

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
#> Hint: The data contains missing values (NA). Ensure that imputation or handling of missing values aligns with your analysis requirements. Note that limma can handle missing values (it just ignores them), and therefore SplineOmics can also handle them.
#> 
#> Ensure the design has no Condition*Time interaction for mode == 'isolated', and includes it for mode == 'integrated'.
#> 
#>  Performing the clustering...
#> For the level:  Constant
#> For the level:  Temp_shift
#> Running this function took 0.0 min
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
#> Generating cluster mean splines for level:  Constant
#> Generating spline plots...
#> Generating cluster mean splines for level:  Temp_shift
#> Generating spline plots...
#> Generating report. This takes a few seconds.
#> 
#>  Info Clustering the hits completed successfully.
#>  Your HTML reports are located in the directory:  /tmp/RtmpNQV3EO/file45aa1da11751 .
#>  Please note that due to embedded files, the reports might be flagged as
#>  harmful by other software. Rest assured that they provide no harm.
```

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
    #> [1] knitr_1.50        dplyr_1.1.4       here_1.0.1        readr_2.1.5      
    #> [5] SplineOmics_0.4.2
    #> 
    #> loaded via a namespace (and not attached):
    #>   [1] RColorBrewer_1.1-3       rstudioapi_0.17.1        jsonlite_2.0.0          
    #>   [4] shape_1.4.6.1            magrittr_2.0.3           farver_2.1.2            
    #>   [7] nloptr_2.2.1             rmarkdown_2.29           GlobalOptions_0.1.2     
    #>  [10] fs_1.6.6                 ragg_1.4.0               vctrs_0.6.5             
    #>  [13] minqa_1.2.8              base64enc_0.1-3          htmltools_0.5.8.1       
    #>  [16] progress_1.2.3           broom_1.0.8              Formula_1.2-5           
    #>  [19] variancePartition_1.38.0 sass_0.4.10              KernSmooth_2.23-26      
    #>  [22] bslib_0.9.0              htmlwidgets_1.6.4        desc_1.4.3              
    #>  [25] pbkrtest_0.5.4           plyr_1.8.9               cachem_1.1.0            
    #>  [28] lifecycle_1.0.4          iterators_1.0.14         pkgconfig_2.0.3         
    #>  [31] Matrix_1.7-3             R6_2.6.1                 fastmap_1.2.0           
    #>  [34] rbibutils_2.3            clue_0.3-66              digest_0.6.37           
    #>  [37] numDeriv_2016.8-1.1      colorspace_2.1-1         S4Vectors_0.46.0        
    #>  [40] rprojroot_2.1.1          textshaping_1.0.1        labeling_0.4.3          
    #>  [43] abind_1.4-8              compiler_4.5.2           bit64_4.6.0-1           
    #>  [46] aod_1.3.3                withr_3.0.2              doParallel_1.0.17       
    #>  [49] backports_1.5.0          BiocParallel_1.42.1      carData_3.0-5           
    #>  [52] gplots_3.2.0             MASS_7.3-65              rjson_0.2.23            
    #>  [55] corpcor_1.6.10           gtools_3.9.5             caTools_1.18.3          
    #>  [58] tools_4.5.2              zip_2.3.3                remaCor_0.0.18          
    #>  [61] glue_1.8.0               nlme_3.1-168             grid_4.5.2              
    #>  [64] checkmate_2.3.3          cluster_2.1.8.1          reshape2_1.4.4          
    #>  [67] generics_0.1.4           gtable_0.3.6             tzdb_0.5.0              
    #>  [70] tidyr_1.3.1              hms_1.1.3                car_3.1-3               
    #>  [73] BiocGenerics_0.54.0      ggrepel_0.9.6            foreach_1.5.2           
    #>  [76] pillar_1.11.0            stringr_1.5.1            vroom_1.6.5             
    #>  [79] limma_3.64.1             circlize_0.4.16          splines_4.5.2           
    #>  [82] lattice_0.22-5           renv_1.1.5               gmp_0.7-5               
    #>  [85] bit_4.6.0                tidyselect_1.2.1         ComplexHeatmap_2.24.1   
    #>  [88] pbapply_1.7-2            reformulas_0.4.1         IRanges_2.42.0          
    #>  [91] svglite_2.2.1            RhpcBLASctl_0.23-42      stats4_4.5.2            
    #>  [94] xfun_0.52                Biobase_2.68.0           statmod_1.5.0           
    #>  [97] matrixStats_1.5.0        stringi_1.8.7            yaml_2.3.10             
    #> [100] boot_1.3-31              evaluate_1.0.4           codetools_0.2-19        
    #> [103] tibble_3.3.0             BiocManager_1.30.26      cli_3.6.5               
    #> [106] systemfonts_1.2.3        Rdpack_2.6.4             jquerylib_0.1.4         
    #> [109] Rcpp_1.1.0               EnvStats_3.1.0           png_0.1-8               
    #> [112] parallel_4.5.2           pkgdown_2.1.3            ggplot2_3.5.2           
    #> [115] prettyunits_1.2.0        ClusterR_1.3.3           bitops_1.0-9            
    #> [118] lme4_1.1-37              mvtnorm_1.3-3            lmerTest_3.1-3          
    #> [121] scales_1.4.0             purrr_1.1.0              crayon_1.5.3            
    #> [124] writexl_1.5.4            fANCOVA_0.6-1            GetoptLong_1.0.5        
    #> [127] rlang_1.1.6
