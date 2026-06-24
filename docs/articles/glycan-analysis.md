# Glycan-analysis

## About this vignette

This tutorial intends to showcase and explain the capabilities of the
**SplineOmics** package by walking through a real and complete glycan
analysis example, from start to finish. **SplineOmics** is explained in
more detail in the **get-started** vignette, where a proteomics example
is covered. This vignette is more focused on showing how glycan data can
be used, and because of this, less details about the overall package are
provided here.

### Data Overview

This dataset originates from a **time-series glycan experiment**
designed to study Chinese Hamster Ovary (CHO) cells. The data is
compositional, meaning that the observed glycoform abundances represent
parts of a whole and are constrained by a constant-sum constraint,
making standard statistical analyses inappropriate without appropriate
log-ratio transformations. The experiment involved cultivating cells in
**eight bioreactors**, with **four bioreactors** subjected to a
temperature shift after 146 hours (experimental condition) and the
remaining **four bioreactors** maintained without a temperature shift
(control condition).

#### Timepoints

Samples were collected at **7 distinct time points** throughout the
experiment, specifically: `"120h"`, `"144h"`, `"168h"`, `"192h"`,
`"216h"`, `"288h"`, and `"336h"` after cultivation start. Each time
point was sampled from all eight bioreactors, but E17_336 is missing,
therefore resulting in a total of **55 samples**.

#### Effects in the Experiment: Reactor and Batch

In this experiment, there are two effects to consider: **Reactor** and
**Batch**.

1.  **Reactor**:
    - This refers to the different bioreactors used for cell
      cultivation, which can exhibit substantial variability.
    - Each reactor was assigned a single condition: either constant
      temperature or temperature-shifted. As a result, **condition and
      reactor are confounded**.
    - Reactor should not be treated as a fixed effect to simply remove
      its influence. Instead, it is treated as a **random effect**,
      which allows us to model its variability appropriately.
2.  **Batch**:
    - This refers to the two separate glycan analysis batches.
    - Batch is considered a **batch effect** with respect to the
      condition (constant temperature vs. temperature-shifted).

Since Condition and Reactor are confounded, the variability due to the
reactors cannot be directly separated from the condition. Instead,
**linear mixed models (LMMs)** are used to attribute the reactor as a
random effect, allowing us to account for its variability while
isolating the effects of the condition. This approach ensures that the
analysis appropriately handles the hierarchical structure of the data
and avoids incorrect conclusions.

In this vignette, we will demonstrate how to use linear mixed models to
address these challenges and properly account for both reactor and batch
effects.

#### Further info

The data matrix comprises **glycoforms as rows** (such as G0/G0 and
none/G0F) and **samples as columns**, providing glycoform measurements
for all time points. Each glycoform as a row stands for a combination of
sugars (glycans) that can be attached to the left and right side of the
product antibody, that we produced in our CHO cell cultivation. For
example G0/G0 means that the glycan G0 was attached to both sides, and
none/G0F means at the left side, there was no glycan, and on the right
side, there was the G0F glycan.

The goal of this experiment is to investigate the effect of a
temperature shift during CHO cell cultivation on the antibody glycan
dynamics over time.

Note: This is not the original dataset, as it has not yet been published
at the time of this vignette’s creation. For demonstration purposes, the
glycans have been randomly shuffled, and 2% has been randomly added or
substracted off each value.

### Analysis Goals

The main objectives of this analysis are:

- **Identify glycans with significant temporal changes**: Among the
  glycoforms measured, the goal is to identify those that exhibit
  significant changes in abundance over time.

- **Cluster glycans based on temporal patterns**: Glycoforms showing
  significant temporal changes (hits) will be grouped into clusters
  based on their time-dependent expression patterns.

- **Assess the impact of temperature shifts on temporal patterns**: The
  analysis will determine whether the temporal patterns of glycoform
  abundance are affected by the temperature shift, i.e., whether
  glycoform abundance dynamics differ over time under temperature shift
  conditions compared to controls.

#### Note

The documentation of all the **SplineOmics** package functions can be
viewed [here](https://csbg.github.io/SplineOmics/reference)

## Load the packages

``` r

library(SplineOmics)
#> Registered S3 method overwritten by 'lme4':
#>   method           from
#>   na.action.merMod car
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
library(compositions) # For clr transforming the glycan data
#> Welcome to compositions, a package for compositional data analysis.
#> Find an intro with "? compositions"
#> 
#> Attaching package: 'compositions'
#> The following objects are masked from 'package:stats':
#> 
#>     anova, cor, cov, dist, var
#> The following object is masked from 'package:graphics':
#> 
#>     segments
#> The following objects are masked from 'package:base':
#> 
#>     %*%, norm, scale, scale.default
```

## Load the files

``` r

data <- read.csv(
    system.file(
        "extdata",
        "glycan_data.csv",
        package = "SplineOmics"
    ),
    stringsAsFactors = FALSE
)

meta <- read.csv(
    system.file(
        "extdata",
        "glycan_meta.csv",
        package = "SplineOmics"
    ),
    stringsAsFactors = FALSE
)

# Set the first column as row names and remove it from the dataframe
rownames(data) <- data[[1]] # Assign first column as row names
data <- data[, -1] # Remove the first column


# Make data to a numeric matrix (required by SplineOmics)
data <- data.matrix(data)
```

### Show top rows of data

``` r

knitr::kable(
    head(data),
    format = "markdown"
)
```

|  | E13_120 | E13_144 | E13_168 | E13_192 | E13_216 | E13_288 | E13_336 | E14_120 | E14_144 | E14_168 | E14_192 | E14_216 | E14_288 | E14_336 | E15_120 | E15_144 | E15_168 | E15_192 | E15_216 | E15_288 | E15_336 | E16_120 | E16_144 | E16_168 | E16_192 | E16_216 | E16_288 | E16_336 | E17_120 | E17_144 | E17_168 | E17_192 | E17_216 | E17_288 | E18_120 | E18_144 | E18_168 | E18_192 | E18_216 | E18_288 | E18_336 | E19_120 | E19_144 | E19_168 | E19_192 | E19_216 | E19_288 | E19_336 | E20_120 | E20_144 | E20_168 | E20_192 | E20_216 | E20_288 | E20_336 |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| none/G0F | 1.0558960 | 1.6535404 | 0.6493199 | 0.4020102 | 0.2991287 | 0.2952674 | 0.4744168 | 1.0534692 | 1.093678 | 0.5601851 | 0.5350867 | 0.3788699 | 0.4035115 | 0.3757274 | 0.8211157 | 0.8911038 | 0.3801903 | 0.4426778 | 0.3843382 | 0.3318367 | 0.3398791 | 0.7996867 | 1.571295 | 1.4926070 | 0.4777732 | 0.6603077 | 0.3305801 | 0.2346606 | 0.7880709 | 1.0201295 | 0.7577885 | 0.1905728 | 0.3328567 | 0.2802764 | 0.7150925 | 0.7753927 | 0.9080968 | 0.3811781 | 0.2896830 | 0.3272081 | 0.2351659 | 1.3833084 | 0.9376732 | 1.1173998 | 0.6469200 | 0.3043375 | 0.2228878 | 0.3791463 | 1.0784873 | 0.8032912 | 0.3756631 | 0.3559814 | 0.2675088 | 0.3162241 | 0.2643130 |
| none/G1F | 0.5323374 | 0.1912681 | 0.3820570 | 0.1920724 | 0.3369441 | 0.6384174 | 0.8548695 | 0.5733259 | 0.567492 | 0.7163753 | 0.4688813 | 0.4682832 | 0.5383108 | 0.6551322 | 0.4758381 | 0.3954215 | 0.3811519 | 0.3056876 | 0.4411230 | 0.5827720 | 0.6590168 | 0.4662977 | 0.621623 | 0.4692532 | 0.4053633 | 0.4627089 | 0.5262178 | 0.4367194 | 0.2626777 | 0.4807747 | 0.4591137 | 0.0997751 | 0.2675994 | 0.7528157 | 0.3142561 | 0.2728702 | 0.4900672 | 0.3237094 | 0.3904656 | 0.5668250 | 0.4722005 | 0.5748934 | 0.6408565 | 0.5941058 | 0.3463255 | 0.2640809 | 0.4779125 | 0.8081467 | 0.3307047 | 0.3004908 | 0.2121393 | 0.3430701 | 0.3647285 | 0.5845660 | 0.6514764 |
| none/G2F | 1.5324508 | 1.8521343 | 0.8379884 | 0.7293874 | 0.4902723 | 0.5382350 | 0.6436627 | 1.5629840 | 1.289147 | 0.9681731 | 0.8381506 | 0.7008038 | 0.6532061 | 0.6719115 | 1.2138739 | 1.1464261 | 0.8146055 | 0.7811525 | 0.6912260 | 0.5626390 | 0.4371626 | 1.2651715 | 1.989103 | 1.5335954 | 0.7930346 | 0.9061705 | 0.4856366 | 0.4323194 | 1.4655321 | 1.3453025 | 1.1834043 | 0.5523239 | 0.5450371 | 0.5750546 | 1.0783679 | 1.1487218 | 1.1215321 | 0.6694769 | 0.5519309 | 0.5116042 | 0.4199132 | 1.7598936 | 1.1849219 | 1.3638929 | 0.7197493 | 0.5763007 | 0.2844577 | 0.5380468 | 1.5920806 | 1.3495185 | 0.8517591 | 0.6509213 | 0.6691747 | 0.6211062 | 0.6178701 |
| G0/G0 | 12.4576400 | 10.2160346 | 10.2302200 | 6.5323299 | 4.5557592 | 4.6899314 | 4.7203765 | 12.3929924 | 10.510602 | 11.4393911 | 9.4008800 | 9.4694548 | 8.1328422 | 7.6767012 | 12.1137799 | 11.9553981 | 11.2695201 | 9.3657469 | 8.0426315 | 6.5983448 | 5.2429438 | 12.0570776 | 10.660712 | 11.3814227 | 10.8993364 | 9.5825590 | 5.9291533 | 4.4084289 | 10.3195660 | 8.3433472 | 9.5614397 | 7.4586037 | 6.6427808 | 8.4710449 | 11.0041856 | 9.8884089 | 10.6157702 | 8.8851596 | 7.7355069 | 6.4930713 | 6.3931685 | 10.3100238 | 8.9498072 | 9.7648715 | 7.6610699 | 6.2008262 | 3.8324832 | 3.4344767 | 10.5438400 | 8.6848225 | 9.8761378 | 8.7988318 | 8.1356305 | 7.7008655 | 8.5356119 |
| G0/G0F | 34.1841189 | 38.3372702 | 41.0179628 | 49.0643714 | 57.0485903 | 55.2748321 | 51.4400433 | 35.7867219 | 39.951440 | 41.3072610 | 45.0540699 | 45.5449446 | 46.7403271 | 46.3506974 | 37.0730784 | 39.6951101 | 40.6190535 | 46.7785275 | 48.3297400 | 55.0964653 | 55.9160088 | 34.3909667 | 35.142868 | 38.1036796 | 41.7931083 | 44.5759546 | 56.0226674 | 63.5428148 | 36.5324214 | 40.7206353 | 41.4229164 | 48.7461735 | 52.2302885 | 49.1391718 | 40.4162712 | 43.1801464 | 42.3152053 | 47.0698711 | 48.8438299 | 50.1411749 | 51.0158271 | 35.9886203 | 40.0629781 | 40.8438051 | 46.0248557 | 52.3444655 | 58.5929353 | 57.2609213 | 36.2618340 | 41.6920256 | 42.9619517 | 45.8755552 | 45.5816057 | 45.9113688 | 46.1812359 |
| G0F/G0F | 1.5218235 | 0.7862833 | 1.1317116 | 0.3237543 | -0.1691198 | 0.9544022 | 1.2641091 | 1.8948004 | 1.654238 | 1.2341346 | 1.2849757 | 1.1123587 | 1.1001650 | 0.6284302 | 1.5736525 | 1.3480386 | 1.6721168 | 1.0828775 | 0.8432043 | 0.8164667 | 1.0990374 | 2.2216793 | 1.773861 | 1.5704807 | 1.3096889 | 1.2557005 | 0.8205026 | 0.5531342 | 1.1133060 | 0.8061578 | 1.1689967 | 0.3151786 | 0.2796022 | -0.1211154 | 1.1292473 | 0.9209109 | 1.3967621 | 1.0388824 | 0.7251202 | 0.7248511 | 0.5504109 | 1.3079025 | 0.7601992 | 1.1123447 | 0.5752250 | 0.4825608 | 0.0937523 | 0.2893026 | 1.2427219 | 0.7371158 | 0.8688332 | 0.8953497 | 0.6032681 | 0.8568143 | 1.0963033 |

### Show top rows of meta

``` r

knitr::kable(
    head(meta),
    format = "markdown"
)
```

| sample_name | Reactor | Time | Condition | Batch |
|:------------|:--------|-----:|:----------|------:|
| E13_120     | E13     |  120 | constant  |     1 |
| E13_144     | E13     |  144 | constant  |     1 |
| E13_168     | E13     |  168 | constant  |     1 |
| E13_192     | E13     |  192 | constant  |     1 |
| E13_216     | E13     |  216 | constant  |     1 |
| E13_288     | E13     |  288 | constant  |     1 |

## Perform EDA (exploratory data analysis)

``` r

# Those fields are mandatory, because we believe that when such a report is
# opened after half a year, those infos can be very helpful.
report_info <- list(
    omics_data_type = "Glycan",
    data_description = "clr transformed timeseries fractional abundance
    glycoform data of CHO cells, batch corrected",
    data_collection_date = "September 2024",
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

Compositional data, such as glycan profiles expressed as relative
abundances, are inherently constrained by a constant sum and therefore
reside in a simplex rather than Euclidean space. This violates the
assumptions of standard statistical methods. To address this, we applied
a centered log-ratio (CLR) transformation, which appropriately maps the
data into real space by accounting for the compositional structure.
While some studies apply a simple log-transformation to stabilize
variance, this approach ignores the relative nature of the data and may
lead to biased or misleading results. CLR transformation is thus the
more appropriate and statistically sound choice for downstream modeling.

``` r

# TRANSPOSE: make samples rows, features columns
data_t <- t(data)

# Handle zeros (simple pseudocount, crude but common)
data_t[data_t == 0] <- 1e-6

# Apply CLR
clr_data_t <- compositions::clr(data_t)

# Transpose back to original shape: features as rows, samples as columns
clr_data <- t(clr_data_t)
clr_data <- unclass(clr_data)
```

``` r

# splineomics now contains the SplineOmics object.
splineomics <- SplineOmics::create_splineomics(
    data = clr_data,
    meta = meta,
    report_info = report_info,
    condition = "Condition", # Column of meta that contains the levels.
    meta_batch_column = "Batch" # For batch effect removal in the plots
)

# Special print.SplineOmics function leads to selective printing
print(splineomics)
#> data:SplineOmics Object
#> -------------------
#> Number of features (rows): 10 
#> Number of samples (columns): 55 
#> Meta data columns: 5 
#> First few meta columns:
#>   sample_name Reactor Time Condition Batch
#> 1     E13_120     E13  120  constant     1
#> 2     E13_144     E13  144  constant     1
#> 3     E13_168     E13  168  constant     1
#> Condition: Condition 
#> No RNA-seq data provided.
#> No annotation provided.
#> No spline parameters set.
#> P-value adjustment method: BH
```

``` r

plots <- SplineOmics::explore_data(
    splineomics = splineomics, # SplineOmics object
    report_dir = report_dir
)
```

[Here](https://drive.google.com/uc?id=1h4aOSwxUHAp5KAoAYVOlL_IRotV4EtY0)
you can see the HTML report of the explore_data() function with the NOT
batch-corrected data, and
[here](https://drive.google.com/uc?id=1op3oCt7eNhrUaofIkAyNJrSqpschU78b)
the report for the batch-corrected data.

## Run limma spline analysis

We set spline parameters manually here (`dof = 2`). Alternatively, set
`spline_params$dof = 0` in
[`update_splineomics()`](https://csbg.github.io/SplineOmics/reference/update_splineomics.md)
to let
[`run_limma_splines()`](https://csbg.github.io/SplineOmics/reference/run_limma_splines.md)
choose the degrees of freedom by leave-one-out cross-validation (see the
**get_started** vignette).

Define the parameters and store them in the `SplineOmics` object:

``` r

splineomics <- SplineOmics::update_splineomics(
    splineomics = splineomics,
    use_array_weights = FALSE,
    # Reactor as random effect.
    design = "~ 1 + Condition*Time + Batch + (1|Reactor)",
    mode = "integrated", # means limma uses the full data for each condition.
    spline_params = list(
        spline_type = c("n"),
        dof = c(2L) # use 0L for automatic selection via LOOCV
    )
)
```

The term ‘(1\|Reactor)’ in the design formula are the random effects.
Therefore, the linear mixed models that we use will model Reactor as a
random effect.

Run the
[`run_limma_splines()`](https://csbg.github.io/SplineOmics/reference/run_limma_splines.md)
function with the updated SplineOmics object:

``` r

splineomics <- SplineOmics::run_limma_splines(
    splineomics = splineomics
)
```

Heteroscedasticity (unequal variance across samples) can affect CLR-
transformed glycan data as well. In this example we set
`use_array_weights = FALSE` for a fixed robustness choice. For your own
data, consider `use_array_weights = NULL` so SplineOmics applies array
weights automatically when the Levene test indicates heteroscedasticity,
or set it to `TRUE` to always use `limma` array weights with robust
`eBayes`. See the
[FAQ](https://csbg.github.io/SplineOmics/articles/faq.md) and
[`create_splineomics()`](https://csbg.github.io/SplineOmics/reference/create_splineomics.md)
for details.

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
[here](https://drive.google.com/uc?id=1jgE8d5Gs12xnMATaitXxp28j2nAKCX1Y).

## Cluster the hits (significant features)

Prepare arguments:

``` r

nr_clusters <- list(
    constant = 2,
    tshifted = 2
)

plot_info <- list(
    y_axis_label = "Glycan fractional abundance",
    time_unit = "hours",
    treatment_labels = list(
        tshifted = "temp shift"
    ),
    treatment_timepoints = list(
        tshifted = 146
    )
)

# Those are not genes for this glycan analysis here, but this argument expects
# the feature names (which usually are gene names, which is why it is called
# like this.
genes <- rownames(data)

plot_options <- list(
    # When meta_replicate_column is not there, all datapoints are blue.
    meta_replicate_column = "Reactor", # Colors the data points based on Reactor
    cluster_heatmap_columns = FALSE
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
```

You can view the generated analysis report of the cluster_hits function
[here](https://drive.google.com/uc?id=1HXKoyOniO6-tDC-uCiGgpaQcjUQ9-zyl).

As discussed before, there are three limma result categories. The
cluster_hits() report shows the results of all three, if they are
present (category 2 and 3 can only be generated when the design formula
contains an interaction effect).

## Session Info

    #> R version 4.6.0 (2026-04-24)
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
    #> [1] compositions_2.0-9 knitr_1.51         dplyr_1.2.1        SplineOmics_0.4.4 
    #> 
    #> loaded via a namespace (and not attached):
    #>   [1] Rdpack_2.6.6             bitops_1.0-9             pbapply_1.7-4           
    #>   [4] writexl_1.5.4            rlang_1.2.0              magrittr_2.0.5          
    #>   [7] clue_0.3-68              GetoptLong_1.1.1         otel_0.2.0              
    #>  [10] matrixStats_1.5.0        compiler_4.6.0           reshape2_1.4.5          
    #>  [13] png_0.1-9                systemfonts_1.3.2        vctrs_0.7.3             
    #>  [16] stringr_1.6.0            pkgconfig_2.0.3          shape_1.4.6.1           
    #>  [19] crayon_1.5.3             fastmap_1.2.0            backports_1.5.1         
    #>  [22] caTools_1.18.3           rmarkdown_2.31           nloptr_2.2.1            
    #>  [25] ragg_1.5.2               purrr_1.2.2              xfun_0.59               
    #>  [28] cachem_1.1.0             jsonlite_2.0.0           progress_1.2.3          
    #>  [31] EnvStats_3.1.0           remaCor_0.0.20           gmp_0.7-5.1             
    #>  [34] BiocParallel_1.46.0      broom_1.0.13             parallel_4.6.0          
    #>  [37] prettyunits_1.2.0        cluster_2.1.8.2          R6_2.6.1                
    #>  [40] stringi_1.8.7            bslib_0.11.0             RColorBrewer_1.1-3      
    #>  [43] limma_3.68.4             boot_1.3-32              car_3.1-5               
    #>  [46] ClusterR_1.3.6           numDeriv_2016.8-1.1      jquerylib_0.1.4         
    #>  [49] Rcpp_1.1.1-1.1           iterators_1.0.14         base64enc_0.1-6         
    #>  [52] IRanges_2.46.0           Matrix_1.7-5             splines_4.6.0           
    #>  [55] tidyselect_1.2.1         rstudioapi_0.19.0        abind_1.4-8             
    #>  [58] yaml_2.3.12              doParallel_1.0.17        gplots_3.3.0            
    #>  [61] codetools_0.2-19         plyr_1.8.9               lmerTest_3.2-1          
    #>  [64] lattice_0.22-9           tibble_3.3.1             Biobase_2.72.0          
    #>  [67] S7_0.2.2                 evaluate_1.0.5           desc_1.4.3              
    #>  [70] bayesm_3.1-7             zip_3.0.0                circlize_0.4.18         
    #>  [73] pillar_1.11.1            BiocManager_1.30.27      tensorA_0.36.2.1        
    #>  [76] carData_3.0-6            KernSmooth_2.23-26       checkmate_2.3.4         
    #>  [79] renv_1.2.3               foreach_1.5.2            stats4_4.6.0            
    #>  [82] reformulas_0.4.4         generics_0.1.4           S4Vectors_0.50.1        
    #>  [85] hms_1.1.4                ggplot2_4.0.3            scales_1.4.0            
    #>  [88] aod_1.3.3                minqa_1.2.8              gtools_3.9.5            
    #>  [91] RhpcBLASctl_0.23-42      glue_1.8.1               tools_4.6.0             
    #>  [94] fANCOVA_0.6-1            robustbase_0.99-7        variancePartition_1.42.0
    #>  [97] lme4_2.0-1               mvtnorm_1.4-1            fs_2.1.0                
    #> [100] grid_4.6.0               tidyr_1.3.2              rbibutils_2.4.1         
    #> [103] colorspace_2.1-2         nlme_3.1-169             Formula_1.2-5           
    #> [106] cli_3.6.6                textshaping_1.0.5        svglite_2.2.2           
    #> [109] ComplexHeatmap_2.28.0    corpcor_1.6.10           DEoptimR_1.2-0          
    #> [112] gtable_0.3.6             sass_0.4.10              digest_0.6.39           
    #> [115] BiocGenerics_0.58.1      pbkrtest_0.5.5           ggrepel_0.9.8           
    #> [118] rjson_0.2.23             htmlwidgets_1.6.4        farver_2.1.2            
    #> [121] htmltools_0.5.9          pkgdown_2.2.0            lifecycle_1.0.5         
    #> [124] GlobalOptions_0.1.4      statmod_1.5.2            MASS_7.3-65
