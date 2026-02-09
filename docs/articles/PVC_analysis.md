# PVC_analysis

## About this tutorial

This tutorial demonstrates the **excursion-analysis** capabilities of
the `SplineOmics` package by walking through a complete, real-world
phosphoproteomics example from start to finish.

In the context of time series data, an **excursion** refers to a
**transient deviation from a baseline** — a concept borrowed from signal
processing. These short-term, local shifts can reveal biologically
meaningful events that are not captured by global trends.

Details of how excursions are detected are described in the section on
the **PVC test**, a method specifically designed to identify such local
patterns.

While the core strength of `SplineOmics` lies in modeling **global
temporal patterns** using smoothing splines, many biologically important
phenomena present as **localized events**. These may not be well
captured by standard spline fits due to their smooth and global nature.

For example, in phosphoproteomics data, you might observe a **sharp,
temporary phosphorylation spike** at a specific timepoint following a
stimulus — a pattern that returns to baseline shortly afterward. Such a
transient peak can indicate a key signaling activation event that is
**functionally critical**, but would be smoothed out or diluted in
traditional spline-based analyses.

To handle these scenarios, we developed the **PVC test** (described
below), which is specifically designed to detect such **sharp, local
excursions** in time-resolved omics data.

### Note 1

The documentation of all the **SplineOmics** package functions can be
viewed [here](https://csbg.github.io/SplineOmics/reference)

### Note 2

This vignette focuses on the **excursion-analysis**. The general
functionalities of `SplineOmics` are explained in the **get-started**
vignette.

## The PVC Test

**PVC** stands for **P**eak, **V**alley, and **C**liff. The PVC test is
a method to identify these distinct local patterns (called excursions)
in time series data using a compound contrast approach within the
**limma** framework.

In this context:

- A **peak** occurs when a timepoint is significantly **higher** than
  both its immediate neighbors.
- A **valley** occurs when a timepoint is significantly **lower** than
  both its neighbors.
- A **cliff** describes a sharp directional shift: when one neighbor is
  similar in value, but the other is **significantly different** (higher
  or lower).

### Method

The PVC analysis is implemented as a **rolling, compound contrast**
applied across the time series. Specifically, a window of three
consecutive timepoints is evaluated at a time:  
`(Tᵢ₋₁, Tᵢ, Tᵢ₊₁)`  
where `Tᵢ` is the center of the window.

The first and last timepoints are excluded from this analysis, since
they lack both left and right neighbors.

For each triplet, a contrast is defined as:

2\*Tᵢ - Tᵢ₋₁ - Tᵢ₊₁

This contrast measures whether the center point `Tᵢ` deviates from both
its neighbors in the same direction:

- If `Tᵢ` is **greater than both neighbors** → large positive contrast →
  indicates a **peak**
- If `Tᵢ` is **lower than both neighbors** → large negative contrast →
  indicates a **valley**
- If `Tᵢ` lies **between** its neighbors → effects cancel → contrast is
  near zero → **no signal**

This contrast is tested using `limma`, and the result is a p-value that
indicates whether the observed pattern is statistically significant for
each feature (e.g., gene, phosphosite).

### Multiple Testing Correction

Because the test is performed independently at each center timepoint
across all features, **multiple testing correction** is applied at the
**timepoint level**. This ensures that significance is not inflated due
to repeated testing across features.

Importantly, PVC does **not** assess whether an entire feature contains
a full peak/valley/cliff pattern over time — it evaluates **each
timepoint independently** for such local patterns.

## Details about the dataset

The example dataset used here involves a **time-series phosphoproteomics
experiment**, where CHO (chinese hamster ovary) cells were cultivated in
three bioreactors (three biological replicates). The experiment includes
the following setup:

- Samples were taken both during the **exponential** and **stationary
  growth phases**.
- Samples were collected in triplicates from each reactor at defined
  timepoints relative to cell feeding:
  - 60 minutes before feeding
  - 15, 60, 90, 120, and 240 minutes after feeding

Note that the dataset was truncated to 1000 rows for file size reasons
and the annotation info, such as gene name, was randomly shuffled.

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
package.

The left part of data contains the numeric values, and the right part
the annotation info, which can be copied in a separate dataframe, as
shown below.

``` r
data <- readRDS(xzfile(system.file(
    "extdata",
    "phosphoproteomics_data.rds.xz",
    package = "SplineOmics"
)))

meta <- read.csv(
    system.file(
        "extdata",
        "phosphoproteomics_meta.csv",
        package = "SplineOmics"
    ),
    stringsAsFactors = FALSE
)

# Extract the annotation part from the dataframe.
first_na_col <- which(is.na(data[1, ]))[1]
annotation <- data |> dplyr::select((first_na_col + 1):ncol(data))
```

### Show top rows of data

``` r
kable(
    head(data),
    format = "markdown"
)
```

| TP01_10 MaxLFQ Intensity | TP01_12 MaxLFQ Intensity | TP01_9 MaxLFQ Intensity | TP02_10 MaxLFQ Intensity | TP02_12 MaxLFQ Intensity | TP02_9 MaxLFQ Intensity | TP03_10 MaxLFQ Intensity | TP03_12 MaxLFQ Intensity | TP03_9 MaxLFQ Intensity | TP04_10 MaxLFQ Intensity | TP04_12 MaxLFQ Intensity | TP04_9 MaxLFQ Intensity | TP05_10 MaxLFQ Intensity | TP05_12 MaxLFQ Intensity | TP05_9 MaxLFQ Intensity | TP06_10 MaxLFQ Intensity | TP06_12 MaxLFQ Intensity | TP06_9 MaxLFQ Intensity | TP07_10 MaxLFQ Intensity | TP07_12 MaxLFQ Intensity | TP07_9 MaxLFQ Intensity | TP08_10 MaxLFQ Intensity | TP08_12 MaxLFQ Intensity | TP08_9 MaxLFQ Intensity | TP09_10 MaxLFQ Intensity | TP09_12 MaxLFQ Intensity | TP09_9 MaxLFQ Intensity | TP10_10 MaxLFQ Intensity | TP10_12 MaxLFQ Intensity | TP10_9 MaxLFQ Intensity | TP11_10 MaxLFQ Intensity | TP11_12 MaxLFQ Intensity | TP11_9 MaxLFQ Intensity | TP12_10 MaxLFQ Intensity | TP12_12 MaxLFQ Intensity | TP12_9 MaxLFQ Intensity | …37 | N..Best.Localization.Probability | T..Protein | T..Index | T..Gene | T..Protein.ID | T..Peptide |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|---:|:---|:---|:---|:---|:---|
| 18.8809 | 18.7333 | 18.736 | 18.6562 | 18.767 | 18.7793 | 18.6554 | 18.7401 | 18.5699 | 18.5877 | 18.7569 | 18.6083 | 18.7291 | 18.7104 | 18.5587 | 18.7192 | 18.7356 | 18.5746 | 18.7935 | 18.6219 | 18.2898 | 18.324 | 18.8198 | 18.4836 | 16.39 | 16.5426 | 18.3209 | 18.395 | 19.0248 | 18.4475 | 18.4724 | 18.7529 | 18.3047 | 18.5518 | 18.6099 | 18.3153 | NA | 0.9186 | tr\|A0A3L7HQ14\|A0A3L7HQ14_CRIGR | A0A3L7HUN2_S117 | Top2b | A0A3L7IN31 | NsPMAQTPPCHLPNIPGVTSPGTLIEDSK;NsPMAQtPPCHLPNIPGVTSPGTLIEDSK;RNsPMAQTPPCHLPNIPGVTSPGTLIEDSK;RNsPMAQTPPCHLPNIPGVTSPGtLIEDSK;RNsPMAQtPPCHLPNIPGVTSPGTLIEDSK |
| 15.3527 | NaN | NaN | NaN | NaN | NaN | 15.7263 | NaN | NaN | 15.4816 | 15.5787 | 15.598 | NaN | 14.9219 | NaN | 15.5459 | NaN | NaN | NaN | NaN | NaN | NaN | NaN | NaN | 15.6507 | 15.3509 | NaN | NaN | NaN | NaN | NaN | 15.3698 | NaN | NaN | NaN | NaN | NA | 1.0000 | tr\|A0A3L7HYQ7\|A0A3L7HYQ7_CRIGR | A0A061I228_S641 | NA | A0A098KXG8 | NFQDsPK |
| 16.3484 | NaN | NaN | NaN | NaN | 16.5751 | 16.6528 | NaN | 16.0512 | NaN | NaN | NaN | NaN | NaN | NaN | NaN | NaN | 16.4891 | 16.9102 | 17.0032 | 16.7853 | NaN | NaN | NaN | 17.1439 | 17.163 | 16.2732 | 16.5507 | 17.2955 | 16.9251 | 16.5331 | 17.6113 | NaN | 16.3595 | NaN | 17.1642 | NA | 1.0000 | tr\|A0A061IGN9\|A0A061IGN9_CRIGR | A0A3L7GUZ0_S522 | Snip1 | A0A3L7HFU7 | QHLENDPGsNEDTDIPK |
| 16.2122 | 15.7589 | 15.8026 | 15.9001 | 15.4823 | 15.7732 | 16.452 | 15.884 | 15.8153 | 16.1033 | 16.2624 | 16.5246 | 15.773 | 15.9233 | 16.1753 | 15.6983 | 15.5666 | 16.251 | 15.7702 | 16.2542 | 16.307 | 15.5049 | 15.7808 | 16.5767 | 15.6377 | 16.0753 | 15.7783 | 15.6542 | 15.9858 | 16.0209 | 15.6547 | 16.1474 | 16.4154 | 16.1201 | 16.0911 | 16.1712 | NA | 0.9601 | tr\|A0A061HXV6\|A0A061HXV6_CRIGR | A0A3L7HQF8_S50 | Tomm34 | A0A3L7HPC6 | ASFITDEEQHAsPR;ASFITDEEQHAsPRPAPQLSR |
| 17.0699 | NaN | 17.3155 | NaN | 16.9908 | 16.9393 | 17.6292 | NaN | 16.9251 | 17.3094 | 17.4428 | 17.6153 | 17.301 | NaN | NaN | 17.2152 | 16.9523 | NaN | 16.1536 | NaN | 16.2584 | 15.6612 | 16.263 | 16.5353 | NaN | NaN | 16.0288 | NaN | NaN | 16.3069 | NaN | NaN | NaN | 16.3798 | 16.2436 | NaN | NA | 0.9561 | tr\|A0A8C2L7B4\|A0A8C2L7B4_CRIGR | A0A061IGN9_S52 | Pus3 | A0A098KXG8 | sRSPLHPTNEVPR;sRsPLHPTNEVPR |
| 15.0928 | NaN | 15.3313 | NaN | NaN | 15.3163 | NaN | 15.3458 | 15.8878 | NaN | NaN | NaN | NaN | NaN | 15.9505 | 15.6368 | NaN | 15.683 | 15.467 | 15.0704 | 15.6614 | 15.2963 | 15.7213 | 15.7003 | 15.0605 | 15.4058 | 16.3211 | 15.7134 | NaN | 15.9043 | 16.0438 | 15.7761 | 15.8954 | 15.5926 | 15.5007 | 16.1923 | NA | 0.9562 | tr\|A0A061HY55\|A0A061HY55_CRIGR | A0A061HV22_S13 | Ints1 | A0A8C2L7B4 | GGYLQGNVsGR |

### Show top rows of meta

``` r
kable(
    head(meta),
    format = "markdown"
)
```

| Sample.ID | Corrected.Sample.ID | Reactor | Time.Point | Phase | Time |
|:---|---:|:---|:---|:---|---:|
| TP01_9 MaxLFQ Intensity | 1 | E09 | TP01 | Exponential | -60 |
| TP01_10 MaxLFQ Intensity | 2 | E10 | TP01 | Exponential | -60 |
| TP01_12 MaxLFQ Intensity | 3 | E12 | TP01 | Exponential | -60 |
| TP02_9 MaxLFQ Intensity | 4 | E09 | TP02 | Exponential | 15 |
| TP02_10 MaxLFQ Intensity | 5 | E10 | TP02 | Exponential | 15 |
| TP02_12 MaxLFQ Intensity | 6 | E12 | TP02 | Exponential | 15 |

### Show top rows of annotation

``` r
kable(
    head(annotation),
    format = "markdown"
)
```

| N..Best.Localization.Probability | T..Protein | T..Index | T..Gene | T..Protein.ID | T..Peptide |
|---:|:---|:---|:---|:---|:---|
| 0.9186 | tr\|A0A3L7HQ14\|A0A3L7HQ14_CRIGR | A0A3L7HUN2_S117 | Top2b | A0A3L7IN31 | NsPMAQTPPCHLPNIPGVTSPGTLIEDSK;NsPMAQtPPCHLPNIPGVTSPGTLIEDSK;RNsPMAQTPPCHLPNIPGVTSPGTLIEDSK;RNsPMAQTPPCHLPNIPGVTSPGtLIEDSK;RNsPMAQtPPCHLPNIPGVTSPGTLIEDSK |
| 1.0000 | tr\|A0A3L7HYQ7\|A0A3L7HYQ7_CRIGR | A0A061I228_S641 | NA | A0A098KXG8 | NFQDsPK |
| 1.0000 | tr\|A0A061IGN9\|A0A061IGN9_CRIGR | A0A3L7GUZ0_S522 | Snip1 | A0A3L7HFU7 | QHLENDPGsNEDTDIPK |
| 0.9601 | tr\|A0A061HXV6\|A0A061HXV6_CRIGR | A0A3L7HQF8_S50 | Tomm34 | A0A3L7HPC6 | ASFITDEEQHAsPR;ASFITDEEQHAsPRPAPQLSR |
| 0.9561 | tr\|A0A8C2L7B4\|A0A8C2L7B4_CRIGR | A0A061IGN9_S52 | Pus3 | A0A098KXG8 | sRSPLHPTNEVPR;sRsPLHPTNEVPR |
| 0.9562 | tr\|A0A061HY55\|A0A061HY55_CRIGR | A0A061HV22_S13 | Ints1 | A0A8C2L7B4 | GGYLQGNVsGR |

### Bring the Inputs into the Standardized Format

``` r
data <- SplineOmics::extract_data(
    # The dataframe with the numbers on the left and info on the right.
    data = data,
    # Use this annotation column for the feature names.
    feature_name_columns = c("T..Gene"),
    use_row_index = TRUE, # makes the feature names unique with row index
    top_row = 1,
    bottom_row = 1000,
    right_col = 36,
    left_col = 1
)
```

### SplineOmics Object

``` r
# Those fields are mandatory, because we believe that when such a report is
# opened after half a year, those infos can be very helpful.
report_info <- list(
    omics_data_type = "PPTX",
    data_description = "Phosphoproteomics data of CHO cells",
    data_collection_date = "February 2024",
    analyst_name = "Thomas Rauter",
    contact_info = "thomas.rauter@plus.ac.at",
    project_name = "DGTX"
)
```

``` r
# splineomics now contains the SplineOmics object.
splineomics <- SplineOmics::create_splineomics(
    data = data,
    meta = meta,
    annotation = annotation,
    report_info = report_info,
    condition = "Phase", # Column of meta that contains the levels.
    meta_batch_column = "Reactor", # For batch effect removal
)
```

## Run the PVC-analysis

``` r
# Check out the documentation of the function under the Reference tab.
pvc_results <- SplineOmics::find_pvc(
    splineomics = splineomics,
    alphas = 0.025,
    padjust_method = "BH",
    support = 1
)
#> design matrix of interest not specified. Assuming a one-group experiment.
#> Warning: Partial NA coefficients for 33 probe(s)
#> Warning: Partial NA coefficients for 195 probe(s)
#> design matrix of interest not specified. Assuming a one-group experiment.
#> Warning: Partial NA coefficients for 55 probe(s)
#> Warning: Partial NA coefficients for 199 probe(s)

pvc_results
#> $Exponential
#> $Exponential$pvc_adj_pvals
#>                  15_vs_neighbors 60_vs_neighbors 90_vs_neighbors
#> 1_Top2b                0.9996591      0.89228318    0.9529197807
#> 2_NA                   1.0000000      1.00000000    0.6654892405
#> 3_Snip1                0.7976077      1.00000000    1.0000000000
#> 4_Tomm34               0.5975031      0.95301168    0.3458807846
#> 5_Pus3                 0.8241813      0.99844582    0.5103974918
#> 6_Ints1                0.7391776      1.00000000    1.0000000000
#> 7_Mlh1                 0.9996591      0.57340016    0.4231947052
#> 8_LOC100750437         0.9996591      0.91476163    0.4987772572
#> 9_Pabpc1               0.5731500      0.46668044    0.3669755656
#> 10_Top2b               0.9136090      0.89453090    0.4358371316
#> 11_Gorasp1             0.7710625      0.76519496    0.3357324904
#> 12_Ints1               0.7239466      0.57955504    0.5486871521
#> 13_Syvn1               0.9996591      0.63814186    0.8053523228
#> 14_Znf280b             0.7528198      0.98300959    0.9701805996
#> 15_Mrnip               0.8241813      0.57340016    0.8085069881
#> 16_Rragc               0.8915063      0.77220390    0.4872552464
#> 17_Gorasp1             0.9957487      0.77660746    0.3154653432
#> 18_Tomm34              0.9173402      0.92922311    0.6987798252
#> 19_LOC100757430        0.6485597      0.89453090    0.6803315820
#> 20_Ubxn1               0.9173402      0.78620962    0.9177365750
#> 21_H671_1g1131         0.6904936      0.40322706    1.0000000000
#> 22_Luzp1               0.9428566      0.94155079    0.8085069881
#> 23_Efs                 0.9996591      0.96609855    0.9271815649
#> 24_Mta2                0.9918588      0.96609855    0.5089919729
#> 25_Nedd1               0.8241813      0.61641970    0.1343299522
#> 26_Gigyf1              0.8539081      0.84244807    0.9562143923
#> 27_Myh9                0.8711111      0.64588748    0.6368831503
#> 28_Caskin2             0.5975031      0.89453090    0.3635388719
#> 29_Papolg              0.5975031      0.89453090    0.6353934056
#> 30_Tfg                 0.9996591      0.91274220    0.6992386870
#> 31_Rpl34               0.4335248      0.39726695    0.3277644801
#> 32_Mideas              0.9996591      0.74871175    0.1583039647
#> 33_Gys1                0.6481524      0.57340016    0.6460689701
#> 34_Arhgef6             0.9085372      0.89453090    0.2796179498
#> 35_Ctdspl2             0.9406578      0.96679510    0.9174607580
#> 36_Ptpn14              0.9206541      0.91266136    0.6803315820
#> 37_Raly                1.0000000      1.00000000    1.0000000000
#> 38_Znhit3              0.8775637      0.99844582    0.8179717078
#> 39_LOC113833392        1.0000000      1.00000000    1.0000000000
#> 40_Luc7l3              1.0000000      1.00000000    0.8930375367
#> 41_Rplp0               1.0000000      1.00000000    1.0000000000
#> 42_Gys1                0.9862244      0.96609855    0.2738360002
#> 43_Rpl22l1             0.9206541      0.65495110    0.2090589517
#> 44_Eif3b               0.8478867      0.89228318    0.5837544074
#> 45_Med26               0.8478867      0.97469270    0.2491244798
#> 46_Mepce               0.9996591      0.96707635    0.5299158386
#> 47_Pdcd11              0.8818738      0.91266136    0.7448404619
#> 48_Twf1                0.7976077      1.00000000    1.0000000000
#> 49_LOC100759640        1.0000000      1.00000000    1.0000000000
#> 50_Wrnip1              0.9937178      0.84503288    0.4987772572
#> 51_Poldip3             0.8340945      0.54140969    0.5103974918
#> 52_Ampd2               0.9409124      0.90654204    0.4184815759
#> 53_Mea1                0.9937178      0.60337996    0.3770190551
#> 54_Dbn1                0.3045423      0.32541914    0.3380338135
#> 55_Snip1               0.7140663      0.62588442    0.4391933166
#> 56_Srsf6               0.9307018      0.89453090    0.8591498257
#> 57_LOC113834282        0.7528198      0.97469270    0.1374908882
#> 58_Map9                0.8702312      0.89453090    0.9464121243
#> 59_Cdc42ep1            0.7528198      0.84244807    0.4084783215
#> 60_Poldip3             0.7619059      0.99844582    0.9187143328
#> 61_LOC100764225        1.0000000      1.00000000    1.0000000000
#> 62_Epb41l2             1.0000000      1.00000000    1.0000000000
#> 63_H671_4g11480        0.4884807      0.41581591    0.3669729209
#> 64_Nbn                 1.0000000      1.00000000    1.0000000000
#> 65_U2surp              1.0000000      1.00000000    1.0000000000
#> 66_Gigyf1              0.5975031      0.87597326    0.4240325835
#> 67_NA                  0.4335248      0.58893964    0.4007519987
#> 68_Luc7l3              1.0000000      0.76571006    0.4034772289
#> 69_LOC100752363        1.0000000      0.15911620    0.3352199052
#> 70_Ampd2               0.9996591      0.39726695    0.1163090667
#> 71_LOC100759640        0.9996591      0.39726695    0.1163090667
#> 72_Stam                0.9996591      0.84244807    0.7566755605
#> 73_Nsfl1c              0.9865255      0.98666147    0.9164066845
#> 74_Pfkfb3              0.9996591      0.89979241    0.9164066845
#> 75_Rad23a              0.9409124      0.89453090    0.4034772289
#> 76_Elf2                0.9996591      0.89453090    0.8085069881
#> 77_Crem                0.8540623      0.99844582    0.6437184324
#> 78_Rragc               0.9996591      0.84518166    0.6516309777
#> 79_Lrrfip2             0.9996591      0.96609855    0.9797810656
#> 80_Zyx                 0.7976077      0.39726695    0.1163090667
#> 81_Lrrfip2             1.0000000      1.00000000    1.0000000000
#> 82_Gatad2b             0.6473547      0.97469270    0.3587160453
#> 83_Bcar1               0.9426037      0.89148887    0.9993185527
#> 84_Ehd1                0.8608189      0.67222253    0.6468176079
#> 85_LOC113834282        0.8241813      0.84503288    0.4945708006
#> 86_Tmem230             0.4884807      0.75726680    0.4263010375
#> 87_Ncbp1               0.8711111      0.91266136    0.7814977065
#> 88_Mllt1               1.0000000      1.00000000    1.0000000000
#> 89_Stk17b              0.9085372      0.51929656    0.1163090667
#> 90_Dlgap4              0.4884807      0.75726680    0.4263010375
#> 91_Papolg              0.6933908      0.60337996    0.3186191377
#> 92_Cyld                0.5971522      0.58893964    0.3321813438
#> 93_Gigyf1              0.9937178      0.89453090    0.7448404619
#> 94_Lrrfip2             0.7027612      0.76571006    0.9798872443
#> 95_Lrrfip2             0.9996591      0.95179760    0.6629729761
#> 96_Rlim                0.7911358      0.92096884    0.9914761578
#> 97_Eif3b               0.8751329      0.94478306    0.9456073392
#> 98_Mphosph10           1.0000000      0.76511779    0.6044440309
#> 99_Gatad2b             0.8539081      0.57340016    0.4945708006
#> 100_Srsf6              0.5731500      0.40573021    0.1542353843
#> 101_Zyx                0.9865255      0.97469270    0.9271815649
#> 102_Mphosph10          0.7765806      0.39726695    0.1175948545
#> 103_Psip1              0.9996591      0.41581591    0.3669755656
#> 104_Fbl                1.0000000      1.00000000    1.0000000000
#> 105_H671_1g2680        0.9996591      0.76519496    0.3442993309
#> 106_Sgtb               0.8350348      0.97580669    0.6199397072
#> 107_Gnl3               0.7619059      0.99478987    0.8364333175
#> 108_Eif3b              0.7619059      0.65495110    0.8053523228
#> 109_Serpinb1           0.7528198      0.60337996    0.3444986722
#> 110_N4bp1              0.9409124      0.98261381    0.5103974918
#> 111_Snip1              0.8540623      0.98444277    0.6647202777
#> 112_Psip1              1.0000000      1.00000000    1.0000000000
#> 113_Mlh1               1.0000000      1.00000000    1.0000000000
#> 114_Bsg                0.9996591      0.39726695    0.2700018411
#> 115_Tnpo1              1.0000000      1.00000000    1.0000000000
#> 116_H671_1g2680        0.9530303      0.44163335    0.1163090667
#> 117_Cbx8               0.9362984      0.57340016    0.2646284031
#> 118_Mideas             0.9530303      0.44163335    0.1163090667
#> 119_Mideas             0.5975031      0.89453090    0.9628343032
#> 120_Dcun1d3            0.5975031      0.89453090    0.9628343032
#> 121_Dlg1               0.8478867      0.57340016    0.2221736132
#> 122_Rad23a             0.7708242      0.80301683    0.5641930786
#> 123_Srsf6              0.3045423      0.59567667    0.8287783448
#> 124_Stx7               0.9206541      0.89827722    0.8267719777
#> 125_Pdcd11             0.5601500      0.76519496    0.9973033256
#> 126_Kiaa1958           0.6657483      0.57955504    0.2543712380
#> 127_Pwp1               0.9377457      0.97864177    0.5103974918
#> 128_Txlng              0.5975031      0.68059529    0.6987798252
#> 129_Junb               0.8241813      0.92070000    0.7874312574
#> 130_LOC100759640       0.7224935      0.99844582    0.8318410441
#> 131_Dbn1               0.7391776      0.84244807    0.8474656873
#> 132_Top2b              0.9307018      0.72336026    0.1975058450
#> 133_Rusc2              0.8443625      0.89453090    0.8140711759
#> 134_NA                 1.0000000      1.00000000    1.0000000000
#> 135_LOC113837251       0.9996591      0.95561081    0.7126478180
#> 136_Fam76b             0.9173402      0.95301168    0.7566755605
#> 137_Ptpn14             0.9206541      0.99844582    0.4987772572
#> 138_Chmp4b             0.7391776      0.93310639    0.6987798252
#> 139_Prpf4b             0.5975031      0.94478306    0.6629729761
#> 140_Eif3b              0.4884807      1.00000000    1.0000000000
#> 141_Nsfl1c             0.8617286      0.84244807    0.8053523228
#> 142_Pdlim7             0.4884807      0.94478306    0.5103974918
#> 143_Rnf113a            0.7619059      0.58902875    0.4014644565
#> 144_Epb41l2            0.8447753      0.63814186    0.3587160453
#> 145_Hnrnpc             0.9409124      0.79747244    0.2480080014
#> 146_LOC113834282       0.5601500      0.97469270    0.8026677103
#> 147_Plekho2            0.9937178      0.76519496    0.2738360002
#> 148_Med26              1.0000000      0.84244807    0.8140711759
#> 149_Arhgef40           0.8327078      0.88113882    0.4945708006
#> 150_NA                 0.9937178      0.78620962    0.6368831503
#> 151_Phf8               0.8241813      0.86446114    0.9464121243
#> 152_Minar1             1.0000000      1.00000000    1.0000000000
#> 153_H671_21690         0.9361973      0.99844582    0.8473517850
#> 154_Arhgef40           0.9865255      0.83239038    0.3458807846
#> 155_Chaf1b             0.9805267      0.91266136    0.7100903489
#> 156_Prpf4b             0.9937178      0.89453090    0.9405569286
#> 157_Znf367             0.8447753      0.41581591    0.2738360002
#> 158_Luzp1              0.9996224      0.99844582    0.9164066845
#> 159_LOC113833882       0.7765806      0.84244807    0.8474656873
#> 160_Hnrnpc             0.4884807      0.59208712    0.6987798252
#> 161_Mepce              0.7976077      0.88610890    0.4888053927
#> 162_Ubxn1              0.8443625      0.87450350    0.7448404619
#> 163_Mllt1              0.9996591      0.81902350    0.3075246846
#> 164_Chaf1b             1.0000000      1.00000000    1.0000000000
#> 165_Raly               1.0000000      1.00000000    1.0000000000
#> 166_Gas2l1             0.7239466      0.41853713    0.2038598949
#> 167_Dlg1               0.9996591      0.58893964    0.1163090667
#> 168_Hoxc10             1.0000000      1.00000000    1.0000000000
#> 169_Gigyf1             1.0000000      1.00000000    1.0000000000
#> 170_Luzp1              0.9937178      0.89453090    0.8474656873
#> 171_Srp72              0.7976077      0.84244807    0.3989353266
#> 172_LOC100771461       0.7619059      0.57340016    0.4052192524
#> 173_Chaf1b             0.8818738      0.78620962    0.9271815649
#> 174_C3H11orf58         0.9996591      0.77083169    0.8053523228
#> 175_Pdcd11             0.8078319      0.65561346    0.7977575454
#> 176_Psip1              0.8540623      0.96609855    0.9848227291
#> 177_Prpf4b             0.9996591      0.89771461    0.3669755656
#> 178_Rnf113a            0.6473547      0.98840974    0.2729399625
#> 179_Irf3               1.0000000      0.64588748    0.5103974918
#> 180_Smim13             0.5682360      0.81938063    0.7074776148
#> 181_Gnl3               0.4884807      0.84244807    0.5981207601
#> 182_Psma5              0.8499250      0.67222253    0.2203807883
#> 183_Ptpn14             0.8818738      0.75726680    0.7809050456
#> 184_Prpf4b             0.4335248      0.39726695    0.5356628179
#> 185_Top2b              0.7528198      0.89228318    0.9271815649
#> 186_Prpf38b            0.4763454      0.91266136    0.7448404619
#> 187_Epb41l2            0.9426037      0.61874884    1.0000000000
#> 188_Eif3b              0.9996591      0.32541914    0.2562757161
#> 189_Hnrnpc             0.8327078      0.73152943    0.2562114509
#> 190_LOC100758278       0.4884807      0.89453090    0.5016274831
#> 191_Prpf4b             0.8241813      0.57955504    0.3183586939
#> 192_Caskin2            0.5731500      0.01144135    0.0001024892
#> 193_LOC100752363       0.9996591      0.84503288    0.6352867542
#> 194_Septin6            0.8540623      0.71646090    0.3154653432
#> 195_Max                0.8789016      0.57340016    0.1975058450
#> 196_Mid1ip1            0.7391776      0.94044689    0.6352867542
#> 197_NA                 0.8540623      0.94208942    0.9201386737
#> 198_Hsph1              1.0000000      1.00000000    0.7068411866
#> 199_Nol7               0.5682360      0.57340016    0.6353934056
#> 200_Raly               0.6473547      0.74503582    0.8094385315
#> 201_Smim13             0.7528198      0.90047417    0.5457472072
#> 202_LOC100757535       0.4884807      0.85102646    0.2611499299
#> 203_Net1               0.4884807      0.93743843    0.3937197943
#> 204_LOC100754077       1.0000000      1.00000000    1.0000000000
#> 205_Snip1              0.5601500      0.96609855    0.2116134333
#> 206_Hnrnpc             1.0000000      1.00000000    0.0237123566
#> 207_Ldlrap1            0.9085372      0.99844582    0.9906710175
#> 208_Luzp1              0.8540623      0.81938063    0.3154653432
#> 209_Rpl26              0.7987719      0.57340016    0.5243063397
#> 210_Epb41l2            0.8179315      0.89771461    0.2738360002
#> 211_Znf367             0.7239466      0.92434829    0.2992400844
#> 212_Dlgap4             0.6213394      0.97469270    0.3390104519
#> 213_Plekho2            1.0000000      1.00000000    1.0000000000
#> 214_Zpr1               0.8443625      0.86446114    0.3169935031
#> 215_Dlgap4             0.9996591      1.00000000    1.0000000000
#> 216_Def6               0.9865255      0.87450350    0.3243591041
#> 217_Eif4ebp2           0.9996224      0.63814186    0.3498488866
#> 218_Eef1b2             0.6485597      0.72336026    0.2728600300
#> 219_Rad23a             0.9362126      0.63814186    0.3974824450
#> 220_Morf4l2            0.6317540      0.63814186    0.8873517778
#> 221_Arhgef40           0.7911358      0.85102646    0.9973033256
#> 222_NA                 0.8443625      0.78296882    0.8179717078
#> 223_LOC100773565       0.5731500      0.62629746    0.4945708006
#> 224_Dus2               0.7585314      0.57340016    0.2682605900
#> 225_Pip4p2             0.9996591      0.89228318    0.3183586939
#> 226_Top2b              0.4937630      0.32541914    0.1492654263
#> 227_Znf280b            0.7365385      0.89453090    0.9690539503
#> 228_Pdcd11             0.5975031      0.62629746    0.7229251756
#> 229_Bckdk              0.7528198      0.81778149    0.2895665722
#> 230_Arhgef40           1.0000000      1.00000000    1.0000000000
#> 231_Mepce              0.5682360      0.32541914    0.1163090667
#> 232_Ccnd3              0.6473547      0.95648901    0.6803315820
#> 233_Phf8               0.7138456      0.75201157    0.9271815649
#> 234_H671_1g2680        0.6971286      1.00000000    1.0000000000
#> 235_Ell                0.9996224      0.70094781    0.4034772289
#> 236_U2surp             0.8540623      0.99844582    0.2941869072
#> 237_Rps10              0.8327078      0.39726695    0.3974824450
#> 238_Ctdspl2            0.7391776      0.50159122    0.1492654263
#> 239_Top2b              1.0000000      1.00000000    0.4028573087
#> 240_Msantd3            0.9551112      0.66173205    0.5329495286
#> 241_Fam76b             1.0000000      1.00000000    0.6498365687
#> 242_Ppp4r3a            0.4884807      0.77660746    0.4014644565
#> 243_Gpatch4            1.0000000      1.00000000    1.0000000000
#> 244_Nudc               0.9996591      0.66227289    0.7292065106
#> 245_Nol7               1.0000000      1.00000000    0.4987772572
#> 246_Plekho2            1.0000000      1.00000000    1.0000000000
#> 247_Prpf4b             0.9862244      0.52311323    0.1975058450
#> 248_Mta2               1.0000000      1.00000000    1.0000000000
#> 249_U2surp             0.7976077      0.89453090    0.8262137970
#> 250_Ubxn1              0.7595219      0.76934373    0.5103974918
#> 251_Rlim               0.8540623      0.89309712    0.8474656873
#> 252_Atat1              0.4335248      0.63814186    0.4987772572
#> 253_Ubxn1              0.7797809      0.97469270    0.8140711759
#> 254_H671_1g2680        0.8241813      0.97469270    0.7032969196
#> 255_eIF2aK2            0.4937630      0.57340016    0.2758043564
#> 256_Skiv2l             0.5975031      0.57340016    0.3480771364
#> 257_Rpl28              0.8447753      0.32541914    0.1896719181
#> 258_LOC100759640       0.8818738      0.39726695    0.1163090667
#> 259_Gatad2b            0.5731500      0.70094781    0.4028573087
#> 260_NA                 0.9996591      0.97469270    0.5302084872
#> 261_Gprasp1            0.8711111      0.88113882    0.3937197943
#> 262_Luzp1              0.4884807      0.90654204    0.3498488866
#> 263_Slc1a5             0.9996591      0.89453090    0.6416169647
#> 264_LOC113834282       0.8540623      0.39726695    0.3757441469
#> 265_Srsf6              0.8179315      0.77660746    0.9562143923
#> 266_Cdc42ep1           0.7976077      0.59650136    0.2992400844
#> 267_Net1               1.0000000      1.00000000    1.0000000000
#> 268_Caskin2            0.9996591      0.84518166    0.5103974918
#> 269_LOC100759640       0.5731500      0.39726695    0.4590470303
#> 270_Mideas             0.9173402      0.99844582    0.6803315820
#> 271_Luzp1              0.5731500      1.00000000    1.0000000000
#> 272_Emd                0.9996591      0.97469270    0.9291830695
#> 273_Plpp6              0.4937630      0.57340016    0.6180316715
#> 274_LOC100759640       0.7528198      0.89453090    0.2720197870
#> 275_Rps7               0.8478867      0.96609855    0.9178825492
#> 276_Fkbp1a             0.5842771      0.99844582    0.4638862180
#> 277_Gatad2b            0.9409124      0.57340016    0.1706411676
#> 278_Znf385a            0.9937178      0.73042151    0.2728600300
#> 279_Arhgef6            0.9173402      0.58893964    0.3458807846
#> 280_Slirp              0.9996591      0.57340016    0.3183586939
#> 281_Skiv2l             0.7528008      0.62588442    0.3390104519
#> 282_H671_21690         0.9937178      0.91266136    0.9273415365
#> 283_Kat8               1.0000000      1.00000000    0.9174607580
#> 284_Nkap               0.9173402      0.89148887    0.4363021752
#> 285_Gsk3b              0.7987719      0.94004182    0.8085069881
#> 286_Ints1              0.9307018      0.74691171    0.4007519987
#> 287_Gas2l1             0.6968896      0.39726695    1.0000000000
#> 288_LOC100759640       0.9996591      0.78620962    0.3095333237
#> 289_Top2b              0.9996591      0.90047417    0.6353934056
#> 290_Kif20b             0.9996591      0.77660746    0.8232359030
#> 291_Phf8               0.9996591      1.00000000    1.0000000000
#> 292_Snip1              0.5731500      0.94208942    0.2738360002
#> 293_Gsk3b              0.9996591      0.88700774    0.7448404619
#> 294_Caskin2            0.9996591      0.88700774    0.7448404619
#> 295_C3H11orf58         0.4884807      1.00000000    1.0000000000
#> 296_Lrch4              1.0000000      1.00000000    1.0000000000
#> 297_LOC113834282       0.9937178      0.85102646    0.8085069881
#> 298_LOC100750407       0.9967931      0.87647129    0.5329495286
#> 299_LOC113833392       0.9996591      0.57340016    0.2738360002
#> 300_LOC113833882       0.8540623      0.89453090    0.2728600300
#> 301_Ldlrap1            0.9173402      0.91476163    0.9798872443
#> 302_Wee1               0.9996224      0.70094781    0.3437425283
#> 303_Caap1              1.0000000      1.00000000    1.0000000000
#> 304_Eif4ebp2           0.9546909      0.98666147    0.6946987247
#> 305_Ripk2              0.9206541      0.99844582    0.6897028754
#> 306_Srp72              1.0000000      1.00000000    0.1691476201
#> 307_Taok2              0.9634942      1.00000000    1.0000000000
#> 308_Nr2f6              0.9206541      0.65495110    0.4987772572
#> 309_Arhgef40           0.9428566      0.66022793    0.3770190551
#> 310_Gys1               0.7987719      0.84244807    0.9733407167
#> 311_Dlg1               0.8241813      0.85826672    0.9909934070
#> 312_Vapb               0.6904936      0.63814186    0.4262173942
#> 313_LOC100757535       1.0000000      1.00000000    1.0000000000
#> 314_Mkrn2              0.9926304      0.85430905    0.8053523228
#> 315_Eif3b              0.5601500      0.39726695    0.4084783215
#> 316_Isyna1             0.9996591      0.54140969    0.1896719181
#> 317_Prpf4b             0.7461617      0.57340016    0.3624569089
#> 318_LOC113833882       0.4884807      0.96609855    0.8351075753
#> 319_Lrch4              0.4884807      0.39726695    0.8053523228
#> 320_Dbn1               0.5975031      0.89228318    0.9470511571
#> 321_Abcf1              0.9362893      0.57340016    0.2116134333
#> 322_Ints1              0.4884807      0.39726695    0.4518029413
#> 323_C3H11orf58         0.8987058      0.58893964    0.3390788295
#> 324_Psma5              0.5337685      0.50159122    1.0000000000
#> 325_Fundc1             0.5975031      0.42185492    0.6724947438
#> 326_Papolg             0.8192612      0.99844582    0.9562143923
#> 327_Mideas             0.9996591      0.98873241    0.6990690966
#> 328_Ubxn1              0.6748108      0.74404228    0.2988269514
#> 329_Synm               0.9996591      0.85102646    0.6417500446
#> 330_Arhgef6            1.0000000      1.00000000    1.0000000000
#> 331_Ptpn14             1.0000000      1.00000000    1.0000000000
#> 332_Pgrmc1             1.0000000      1.00000000    1.0000000000
#> 333_Myh9               0.9996224      0.96826768    0.8085069881
#> 334_Etv3               1.0000000      0.74796187    0.4358371316
#> 335_Ip6k1              1.0000000      1.00000000    0.4919700587
#> 336_Luzp1              0.9996591      0.93838860    0.9771317229
#> 337_Ptpn14             0.8241813      0.76519496    0.8474656873
#> 338_Caskin2            0.8353686      0.77660746    0.9628343032
#> 339_Chaf1b             1.0000000      1.00000000    1.0000000000
#> 340_Ubxn1              1.0000000      1.00000000    1.0000000000
#> 341_Ube2c              0.6485597      0.89453090    0.5329495286
#> 342_Gins2              0.4884807      0.41581591    0.4945708006
#> 343_Nlgn2              0.5601500      0.74796187    0.6182392042
#> 344_Nf2                0.9996591      0.84244807    0.7100903489
#> 345_Pip4p2             0.6971286      0.65688382    0.7566755605
#> 346_Emd                0.4884807      0.57340016    0.9562143923
#> 347_Top2b              0.9361973      0.84244807    0.7774012636
#> 348_Trim35             0.9173402      0.99844582    0.9798872443
#> 349_NA                 0.8818738      0.97469270    0.8140711759
#> 350_NA                 0.9996591      1.00000000    1.0000000000
#> 351_Mideas             0.9996591      0.61840397    0.3093690624
#> 352_Gas2l1             0.5731500      0.90047417    0.4240325835
#> 353_Ampd2              0.8540623      0.77660746    0.9177365750
#> 354_Calu               1.0000000      0.98168118    0.8053523228
#> 355_Fam76b             1.0000000      0.99844582    0.8364333175
#> 356_Dlg1               0.9996591      0.91476163    0.9562143923
#> 357_Srsf6              0.9551112      0.70094781    0.3458807846
#> 358_Chaf1b             0.9551112      0.70094781    0.3458807846
#> 359_Dbn1               0.9865255      0.97469270    0.8474656873
#> 360_Tcf25              0.8617286      0.70125132    0.3154653432
#> 361_Psip1              0.9996591      0.89453090    0.4459274268
#> 362_Cnpy3              0.7987719      0.73046398    0.1896719181
#> 363_LOC100759640       0.9409124      0.57340016    0.3587160453
#> 364_Zyx                0.8818738      0.90581890    1.0000000000
#> 365_Lrch4              0.5975031      0.63315255    0.7781328172
#> 366_Bola1              1.0000000      1.00000000    1.0000000000
#> 367_Znf385a            0.9996591      0.64588748    0.4945708006
#> 368_Kif20b             0.8443625      0.57340016    0.2895665722
#> 369_Ell                0.8617286      1.00000000    1.0000000000
#> 370_Ell                0.9307018      0.93310639    1.0000000000
#> 371_Srsf6              0.9996591      0.95236696    0.4116522570
#> 372_Pwp1               0.9996591      0.61840397    0.4034772289
#> 373_Def6               0.9556619      0.70094781    0.5512163402
#> 374_Cbx8               0.9996591      1.00000000    1.0000000000
#> 375_Ddx51              0.9362984      0.94208942    0.5274679042
#> 376_Psip1              0.9937178      0.49257001    0.5689471290
#> 377_Arhgef40           0.9937178      0.39726695    0.1175948545
#> 378_Raly               0.5556805      0.63814186    0.3587160453
#> 379_NA                 0.1475220      0.15274857    0.2090589517
#> 380_Lrrfip2            0.4884807      0.32541914    0.1682782731
#> 381_Gnl3               0.9996591      0.85430905    0.9562143923
#> 382_Caskin2            0.5731500      0.39726695    0.6897028754
#> 383_Rragc              0.7619059      0.63814186    0.4945708006
#> 384_Caskin2            0.1475220      0.15911620    0.1650552773
#> 385_Bcar1              0.9996591      0.77220390    0.7589135521
#> 386_Homer3             0.9996591      0.84244807    0.8664983034
#> 387_Luzp1              0.9996591      0.95561081    0.9464121243
#> 388_N4bp1              1.0000000      1.00000000    1.0000000000
#> 389_Ppp4r3a            0.9996591      0.97692950    0.7613779935
#> 390_H671_1g2680        1.0000000      1.00000000    1.0000000000
#> 391_Gnl3               0.7710625      0.62588442    0.1956873478
#> 392_Top2b              0.9996591      0.99844582    0.9562143923
#> 393_Oser1              1.0000000      0.85102646    0.3907241794
#> 394_Snrk               0.9996591      0.99844582    0.4987772572
#> 395_Kat8               1.0000000      1.00000000    1.0000000000
#> 396_Raver1             1.0000000      1.00000000    1.0000000000
#> 397_Pdcd11             0.8447753      0.90047417    0.8268537217
#> 398_Rps20              0.9996591      0.58893964    0.2738360002
#> 399_Bsg                1.0000000      1.00000000    1.0000000000
#> 400_Raly               0.9996591      0.89453090    0.2796179498
#> 401_Pdcd2              0.8447753      0.91560726    0.1706411676
#> 402_Caskin2            0.8540623      0.91476163    0.5302084872
#> 403_LOC100773571       0.9362893      0.76934373    0.2441304415
#> 404_Papolg             0.7528198      0.89309712    0.7046018595
#> 405_LOC100757535       0.8540623      0.83239038    0.9909934070
#> 406_Caap1              0.9996591      0.57340016    0.1163090667
#> 407_Psip1              0.8540623      0.76701587    0.9271815649
#> 408_Dbn1               0.9996591      0.42243195    0.1650552773
#> 409_Mta2               0.7671830      0.84244807    0.5902417364
#> 410_Abcf1              0.9937178      0.89148887    1.0000000000
#> 411_LOC100754108       0.9206541      0.98168118    0.9185882056
#> 412_Slirp              0.9996591      0.89228318    0.4111144282
#> 413_Nelfa              0.9173402      0.98168118    0.9426420690
#> 414_Aggf1              0.9996591      0.62588442    0.4945708006
#> 415_Bap1               0.9173402      0.98261381    0.7448404619
#> 416_Luc7l3             0.6213394      0.57340016    0.3669729209
#> 417_Rrp1               0.5731500      0.57955504    0.8873517778
#> 418_Wrnip1             0.9926304      1.00000000    1.0000000000
#> 419_NA                 1.0000000      1.00000000    0.9973033256
#> 420_Abcf1              0.9996591      0.95871390    0.6987798252
#> 421_Cluap1             0.3045423      0.89309712    0.8474656873
#> 422_Hnrnpc             0.9937178      0.84503288    0.2738360002
#> 423_Ptpn1              0.9504850      0.81115498    0.1706411676
#> 424_Myh9               0.4335248      0.50159122    0.5103974918
#> 425_U2surp             0.7359970      0.63814186    0.5103974918
#> 426_NA                 0.9996591      0.96609855    0.9291830695
#> 427_Arhgef40           1.0000000      1.00000000    0.7448404619
#> 428_Chaf1b             0.9206541      0.80342508    0.3989353266
#> 429_Prpf4b             0.9996591      0.73710357    0.6838225235
#> 430_Epb41l2            1.0000000      1.00000000    1.0000000000
#> 431_Eif3b              0.9862244      0.76511779    0.2441304415
#> 432_Isyna1             0.9996591      0.70125132    0.1542353843
#> 433_U2surp             0.7528198      0.57340016    0.4187510969
#> 434_LOC100765020       0.9996591      0.84244807    0.3676390944
#> 435_Arhgef6            0.9996224      0.84503288    0.3444986722
#> 436_Ptpn1              0.9173402      0.58038856    0.9733407167
#> 437_Prpf4b             0.7619059      0.75201157    0.7347336557
#> 438_Rpl35a             0.5731500      0.65688382    0.8085069881
#> 439_Prpf4b             1.0000000      1.00000000    1.0000000000
#> 440_Zyx                1.0000000      1.00000000    0.7046018595
#> 441_Dbn1               0.4884807      0.76672626    0.4945708006
#> 442_Chaf1b             1.0000000      0.47674732    0.5902417364
#> 443_LOC113834282       0.7976077      0.60337996    0.4927937173
#> 444_Gpsm2              0.9428566      0.90374179    0.7384669751
#> 445_LOC100757535       0.5731500      0.84503288    0.6146985018
#> 446_Cfap410            0.7619059      0.93310639    0.7977575454
#> 447_Epb41l2            0.5975031      0.57340016    0.2125739710
#> 448_Ncbp1              0.9996591      0.72336026    0.3444986722
#> 449_Pacsin1            0.9307018      0.65440610    0.9623469453
#> 450_Cstf2              0.9173402      0.90964926    0.8094385315
#> 451_LOC100769437       0.9996591      0.89453090    0.4391933166
#> 452_eIF2aK2            0.7911358      0.75726680    0.7589135521
#> 453_Kiaa1191           0.9996591      0.76519496    0.7977575454
#> 454_Mepce              0.9361973      0.93838860    0.8094385315
#> 455_Cbx8               0.4884807      0.57340016    0.8474656873
#> 456_Eed                0.7619059      0.61924693    0.1583039647
#> 457_Cdc42ep1           0.3045423      0.76519496    0.3474399848
#> 458_Lrrfip2            0.7619059      0.83239038    0.5103974918
#> 459_Pacsin1            0.8327078      0.83239038    0.6368831503
#> 460_Gpatch4            0.8818738      0.61840397    0.4111144282
#> 461_Plin4              0.8241813      0.74871175    0.4889794954
#> 462_NA                 0.7710625      0.44694320    0.9798872443
#> 463_Snip1              0.9307018      0.84436118    0.2738360002
#> 464_Cyld               0.8540623      0.73042151    0.6654892405
#> 465_Plin4              0.9996591      0.74871175    0.3557894948
#> 466_Twf1               1.0000000      1.00000000    1.0000000000
#> 467_LOC113834282       1.0000000      1.00000000    1.0000000000
#> 468_Snip1              0.9996591      0.88113882    0.9253735310
#> 469_Ppp4r3a            1.0000000      1.00000000    1.0000000000
#> 470_Psip1              0.8702312      0.96609855    0.5423506258
#> 471_Dnajc5             0.8478867      0.78620962    0.4945708006
#> 472_Phf8               0.9996591      0.41581591    0.1163090667
#> 473_Bola1              0.9085372      0.93310639    0.9973033256
#> 474_Cdc42ep1           0.4884807      0.32541914    0.3390104519
#> 475_Eif4ebp2           0.6160334      0.91266136    0.7083822981
#> 476_Prpf38b            0.7976077      0.63814186    0.6187384142
#> 477_Klhl26             0.9307018      0.89453090    0.3908588040
#> 478_Hsph1              1.0000000      1.00000000    1.0000000000
#> 479_Snip1              0.5654106      0.95648901    0.3474399848
#> 480_Caskin2            0.9996591      0.69844769    0.3458807846
#> 481_Plpp6              0.8915063      0.39726695    0.2738360002
#> 482_NA                 0.8478867      0.39726695    1.0000000000
#> 483_Mlh1               0.4884807      0.62588442    0.8053523228
#> 484_Gys1               0.4884807      0.76519496    0.9973033256
#> 485_Tfg                0.6473547      0.99844582    0.4034772289
#> 486_Arhgef6            1.0000000      1.00000000    0.3277644801
#> 487_Mphosph10          1.0000000      1.00000000    0.3277644801
#> 488_Hoxc10             0.6267469      0.32541914    0.1272404753
#> 489_LOC100759640       0.9771509      0.96826768    0.9914761578
#> 490_Arhgef40           1.0000000      1.00000000    1.0000000000
#> 491_Dnajc5             0.8711111      0.47176584    0.4919700587
#> 492_Tbc1d23            0.9996591      1.00000000    1.0000000000
#> 493_Ubxn1              0.4884807      0.81301115    0.2738360002
#> 494_Rab1a              1.0000000      1.00000000    1.0000000000
#> 495_Eif3b              0.8192612      0.89148887    0.7046018595
#> 496_Tceal8             0.9996591      0.99844582    0.7942831879
#> 497_Dlgap4             0.6968896      0.99844582    0.3770190551
#> 498_Smim13             0.9996591      0.89148887    0.2649667117
#> 499_NA                 0.6968896      0.57340016    0.4240325835
#> 500_Lrch4              0.9937178      0.86446114    0.1492654263
#> 501_Bola1              0.7619059      0.95179760    0.2895665722
#> 502_NA                 0.9926304      0.44694320    0.1343299522
#> 503_Ptpn14             0.9996591      0.96707922    0.7448404619
#> 504_LOC100759640       0.4884807      0.84138511    0.3779801995
#> 505_Rps10              0.5601500      0.61840397    0.8611757851
#> 506_Top2b              0.8156605      0.87450350    0.4159943815
#> 507_Ssr3               1.0000000      0.99844582    0.5329495286
#> 508_Homer3             0.8327078      0.98261381    0.6803315820
#> 509_Phf8               0.7976077      0.74796187    0.3937197943
#> 510_LOC100767716       0.7976077      0.78620962    0.2225177202
#> 511_Xpa                0.5731500      0.88113882    0.5696238511
#> 512_H671_21690         0.9996591      0.94208942    0.4945708006
#> 513_LOC100769471       0.8128341      0.77660746    0.3154653432
#> 514_Gas2l1             0.4884807      0.57340016    0.6353934056
#> 515_Luzp1              0.3045423      0.50709661    0.2090589517
#> 516_Gpbp1              0.5731500      0.76519496    0.9623341869
#> 517_Gatad2b            1.0000000      1.00000000    1.0000000000
#> 518_Gys1               0.8632347      0.95016268    0.5025865483
#> 519_Top2b              0.9996591      0.97469270    0.8262137970
#> 520_LOC100757535       0.5971522      1.00000000    1.0000000000
#> 521_Lpcat4             0.5975031      0.75726680    0.8967754084
#> 522_Arhgef6            0.4251762      0.87647129    0.8611757851
#> 523_Cavin3             0.9996591      0.97411817    0.9642694022
#> 524_Gpatch4            0.9996591      0.71801079    0.7514243941
#> 525_Prpf38b            0.4884807      0.63814186    0.8591498257
#> 526_Timm8a             0.4116137      0.49183911    0.2738360002
#> 527_Cavin3             0.8241813      0.88113882    0.5103974918
#> 528_Mkrn2              0.8039534      0.75726680    0.3458807846
#> 529_Oser1              0.5601500      0.57340016    0.3587160453
#> 530_Gsk3b              0.4116137      0.49183911    0.2738360002
#> 531_Eef1b2             0.6933908      1.00000000    1.0000000000
#> 532_Ampd2              0.8241813      0.99652195    0.7361793136
#> 533_Lrrfip2            0.9085372      0.74796187    0.6324435208
#> 534_Ring1              0.4884807      0.89453090    0.6654892405
#> 535_Rlim               1.0000000      1.00000000    1.0000000000
#> 536_LOC100759640       0.9409124      0.99844582    0.5512163402
#> 537_LOC100759640       0.7619059      0.84764375    0.9630822526
#> 538_Atp5pf             0.9146177      0.96609855    0.6314147294
#> 539_Max                0.6499875      0.76519496    0.8094385315
#> 540_Bap1               0.9337346      0.93190716    0.9234152431
#> 541_Nsfl1c             0.4884807      0.57955504    0.1163090667
#> 542_Prpf4b             1.0000000      1.00000000    0.9354251433
#> 543_LOC100757535       0.7911358      0.84244807    0.9273415365
#> 544_Mtmr10             0.9937178      0.84764375    0.4028573087
#> 545_Hoxc10             1.0000000      1.00000000    1.0000000000
#> 546_Trim35             1.0000000      1.00000000    1.0000000000
#> 547_Eif4ebp2           0.7987719      0.44808277    0.2941869072
#> 548_Dlgap4             0.4884807      0.57340016    1.0000000000
#> 549_Gys1               0.5975031      0.54140969    0.4084783215
#> 550_Sgtb               0.9996591      0.89309712    0.5302084872
#> 551_Eri2               0.7391776      0.89309712    0.9460896228
#> 552_Ccnd3              0.9996591      1.00000000    1.0000000000
#> 553_Smim13             0.6485597      0.89453090    0.3390104519
#> 554_Snrk               0.6738408      0.90442509    0.2738360002
#> 555_Caskin2            0.4884807      0.39726695    0.6266024511
#> 556_Pdcd11             1.0000000      1.00000000    1.0000000000
#> 557_Pgam5              0.9996591      0.96609855    0.8473517850
#> 558_Mphosph10          0.9085372      0.88246087    0.8474656873
#> 559_Mideas             0.9996591      0.99844582    0.7865438845
#> 560_Top2b              0.9996591      0.84503288    0.9291830695
#> 561_LOC100763014       0.7057658      0.47674732    0.3169935031
#> 562_Snip1              0.9173402      0.84244807    0.5778503376
#> 563_Ubxn1              1.0000000      0.78620962    0.5457472072
#> 564_LOC100750407       1.0000000      1.00000000    1.0000000000
#> 565_Morf4l2            1.0000000      1.00000000    1.0000000000
#> 566_Ctdspl2            0.9996591      0.89453090    0.9291830695
#> 567_Cwf19l1            0.5193881      0.50159122    0.7764469293
#> 568_Eef1b2             0.8751329      1.00000000    1.0000000000
#> 569_C1H12orf45         0.8751329      1.00000000    1.0000000000
#> 570_Znf367             1.0000000      1.00000000    1.0000000000
#> 571_Ankrd34a           0.7886993      0.57340016    0.5319516213
#> 572_Mllt11             0.7391776      0.97692950    0.9701805996
#> 573_LOC100774792       0.9996591      0.98300959    0.8179717078
#> 574_NA                 0.5975031      0.99652195    0.3610206615
#> 575_Cbx8               0.7461617      0.97469270    0.6182392042
#> 576_Bckdk              0.7391776      0.74703751    0.3365248466
#> 577_Snip1              1.0000000      1.00000000    1.0000000000
#> 578_Nsfl1c             1.0000000      1.00000000    1.0000000000
#> 579_Gas2l1             1.0000000      0.61641970    0.6928203220
#> 580_Nudc               0.7528198      0.63315255    0.7604572285
#> 581_Epb41l2            0.9307018      1.00000000    1.0000000000
#> 582_Mtmr6              0.6726106      0.69493031    0.8094385315
#> 583_Znf668             0.7619059      0.32541914    0.2116134333
#> 584_Hsph1              1.0000000      1.00000000    0.5329495286
#> 585_LOC113834282       0.8447753      0.97469270    0.9271815649
#> 586_Ctdspl2            1.0000000      1.00000000    1.0000000000
#> 587_Foxf1              0.8711111      0.69097390    0.3183586939
#> 588_Luzp1              0.8539081      0.98226136    0.7083822981
#> 589_Xpa                0.9937178      0.93770530    0.9701805996
#> 590_Psip1              0.8241813      0.76934373    0.3974824450
#> 591_Rbm7               0.4884807      0.50159122    0.8473517850
#> 592_Mtrex              0.8705424      0.44694320    0.9174607580
#> 593_Arhgef40           1.0000000      1.00000000    1.0000000000
#> 594_Plekho2            0.8818738      0.88113882    0.9623341869
#> 595_Bckdk              0.5975031      0.78620962    0.6629729761
#> 596_Dut                0.4884807      0.57340016    0.9104099835
#> 597_Abcf1              0.4884807      0.57340016    0.9104099835
#> 598_Txnl1              0.9996591      0.94478306    0.8085069881
#> 599_Nudc               0.8818738      0.94517210    0.8311625987
#> 600_Sh3gl1             0.7994190      0.39726695    0.1709825636
#> 601_Gatad2b            0.9996591      0.89309712    0.9973033256
#> 602_Homer3             1.0000000      1.00000000    0.6182392042
#> 603_Septin6            0.9996591      0.93838860    0.8271678824
#> 604_Smim13             0.6971286      0.34703515    0.3175586501
#> 605_Arhgef40           0.9996591      0.99480848    0.5457472072
#> 606_Rpl32              0.9996591      1.00000000    1.0000000000
#> 607_Tomm34             0.6968896      0.46668044    0.5329495286
#> 608_Mlh1               0.7619059      0.68895932    0.2555563503
#> 609_Tbcc               1.0000000      1.00000000    1.0000000000
#> 610_Eif3d              1.0000000      1.00000000    0.7589135521
#> 611_Snrk               1.0000000      0.96609855    0.6803315820
#> 612_Bckdk              1.0000000      1.00000000    0.9456073392
#> 613_Wdr3               1.0000000      1.00000000    1.0000000000
#> 614_LOC100757535       0.9634942      0.81115498    0.8794786048
#> 615_Dlg1               0.4116137      0.64588748    0.4987772572
#> 616_LOC100767716       0.9996591      0.93998435    0.4945708006
#> 617_Hnrnpc             0.4335248      0.53918279    0.9456073392
#> 618_Mphosph10          0.5731500      0.58902875    0.9273415365
#> 619_Eif3b              0.4884807      0.57340016    0.6724947438
#> 620_Emd                0.5975031      0.84244807    0.4987772572
#> 621_Txlng              1.0000000      1.00000000    1.0000000000
#> 622_Prpf4b             0.7391776      0.65495110    0.5103974918
#> 623_Rlim               0.9862244      0.81115498    0.6241079841
#> 624_Eef1b2             0.9996591      0.84503288    0.8094385315
#> 625_Def6               0.4335248      0.97469270    0.3344749025
#> 626_LOC100765020       0.9996591      0.95016268    0.3989353266
#> 627_U2surp             0.9937178      0.94004182    0.4919700587
#> 628_Elf2               0.9926304      0.93144523    0.4987772572
#> 629_Slc1a5             0.7528198      0.97469270    0.5696238511
#> 630_NA                 0.8540623      0.65495110    0.9252842391
#> 631_Tfg                0.7391776      0.97469270    0.9529197807
#> 632_Top2b              0.5842771      0.89228318    0.5319516213
#> 633_Pip4p2             0.9362984      0.90047417    0.9672759026
#> 634_Cdc42ep1           1.0000000      1.00000000    1.0000000000
#> 635_Hsph1              0.9937178      0.89453090    1.0000000000
#> 636_Twf1               0.9085372      1.00000000    1.0000000000
#> 637_Nbn                0.5975031      0.84244807    0.8377870644
#> 638_Psmd4              0.9206541      0.39726695    0.1492654263
#> 639_Bap1               0.9996591      0.61840397    0.3147010633
#> 640_Mepce              0.5842771      0.99844582    0.9863023699
#> 641_Mideas             1.0000000      1.00000000    1.0000000000
#> 642_LOC100759640       0.5731500      0.39726695    0.6352867542
#> 643_Epb41l2            0.6738408      0.93309878    0.7589135521
#> 644_Sav1               0.5731500      0.44694320    0.4014644565
#> 645_Prpf4b             0.5975031      0.57340016    0.4834221588
#> 646_Gnas               0.7619059      0.78620962    0.4823171066
#> 647_Mllt1              0.7585314      0.97469270    0.4987772572
#> 648_Poldip3            0.4335248      0.90047417    0.5655994197
#> 649_Aldoa              0.7391776      0.57340016    0.4830442104
#> 650_Rbbp8              0.9206541      0.89453090    0.9914761578
#> 651_LOC113834282       0.5731500      0.44694320    0.9071163669
#> 652_Gys1               0.9996591      0.93450401    0.9456073392
#> 653_Hnrnpc             0.8577026      0.51929656    0.3458807846
#> 654_Vps35              0.6445819      0.39726695    0.1870155117
#> 655_Miga2              0.9738333      0.81115498    0.8474656873
#> 656_Epb41l2            0.9996591      0.96609855    0.9973033256
#> 657_Tob2               1.0000000      1.00000000    1.0000000000
#> 658_Lamtor1            0.9996591      0.89188709    0.9984999026
#> 659_LOC100759640       0.6473547      0.65495110    0.6606817456
#> 660_Epb41l2            1.0000000      1.00000000    1.0000000000
#> 661_Rlim               1.0000000      1.00000000    0.9174607580
#> 662_Gys1               0.4884807      0.69361567    0.2758043564
#> 663_LOC100750437       0.9085372      0.84244807    0.3175586501
#> 664_NA                 0.7528198      0.50159122    0.1163090667
#> 665_Nbn                0.9996591      0.63814186    0.6956285583
#> 666_Tyw3               0.9085372      0.96685495    0.6796661112
#> 667_Gas2l1             0.5975031      0.92644080    0.3458807846
#> 668_Fus                0.9173402      1.00000000    1.0000000000
#> 669_Prpf38b            0.9937178      0.89453090    0.4945708006
#> 670_Calu               1.0000000      1.00000000    1.0000000000
#> 671_Rras2              0.9361973      0.91299684    0.3390104519
#> 672_Prpf4b             0.4884807      0.41581591    0.8967754084
#> 673_Nelfa              0.5731500      0.88113882    0.7036922531
#> 674_LOC100754077       0.5536475      0.89228318    0.6992386870
#> 675_Rbm28              0.8711111      0.99844582    0.8004961183
#> 676_Nsfl1c             0.6317540      0.86446114    0.8726332264
#> 677_Rnf126             0.3686822      0.81576416    0.9906710175
#> 678_Eme1               0.6471239      0.75726680    0.9174607580
#> 679_Nbn                0.4335248      0.94208942    0.1163090667
#> 680_Eif4ebp2           0.8540623      0.81647594    0.2738360002
#> 681_Wee1               0.5536475      0.42243195    0.1492654263
#> 682_Prpf38b            0.3045423      0.34703515    0.8094385315
#> 683_Luzp1              0.5975031      0.94478306    0.6164629145
#> 684_Gas2l1             0.7976077      0.39726695    0.1170839101
#> 685_Pdcd11             0.5975031      0.91266136    0.8053523228
#> 686_Chaf1b             0.3045423      0.22149083    0.6803315820
#> 687_Pycr1              0.4116137      0.39726695    0.6647202777
#> 688_Phf8               0.8789937      0.49995158    0.5105655186
#> 689_Raver1             0.9206541      0.97469270    0.8754200386
#> 690_Dbn1               0.7057658      0.72336026    0.8094385315
#> 691_Dut                0.8327078      0.41581591    0.1163090667
#> 692_Prpf4b             0.9996591      0.97469270    0.6790956907
#> 693_Prpf4b             0.9996591      0.89148887    0.6187384142
#> 694_Efs                0.4884807      0.88113882    0.6199397072
#> 695_NA                 0.5731500      0.97469270    0.3390104519
#> 696_Ppp2r5b            0.6485597      0.99844582    0.7448404619
#> 697_Caskin2            0.4335248      0.58893964    0.2007870458
#> 698_Arhgef40           1.0000000      1.00000000    1.0000000000
#> 699_Zyx                0.3533344      0.57340016    0.1163090667
#> 700_Mphosph10          0.9996591      0.58893964    0.1395399686
#> 701_LOC113833392       1.0000000      0.84503288    0.8474656873
#> 702_Cdc42ep1           0.9206541      0.74345593    0.4014644565
#> 703_Snrpa1             1.0000000      1.00000000    1.0000000000
#> 704_Ncbp1              0.7976077      0.86641680    0.9973033256
#> 705_Gas2l1             0.6904936      1.00000000    1.0000000000
#> 706_Gas2l1             0.9379294      0.83239038    0.8677435804
#> 707_Bap1               0.7239466      0.83239038    0.6724077964
#> 708_LOC100759640       0.9937178      0.69571122    0.7229251756
#> 709_Cherp              1.0000000      0.97469270    0.9909934070
#> 710_Nbn                0.0104809      0.16977096    0.5183821396
#> 711_LOC100759640       0.9996591      1.00000000    1.0000000000
#> 712_NA                 0.8540623      0.83239038    0.4014644565
#> 713_Eif3b              1.0000000      1.00000000    1.0000000000
#> 714_Miga2              0.4884807      0.84764375    0.6458262979
#> 715_Prpf4b             0.7765806      0.76519496    0.1542353843
#> 716_Dbn1               0.8443625      0.58893964    0.1272404753
#> 717_Ppp2r5b            0.7765806      0.89453090    0.3138611909
#> 718_Exosc9             1.0000000      1.00000000    0.1492654263
#> 719_Eif3b              0.9173402      0.96609855    0.4184815759
#> 720_Ripk2              0.9996591      0.97469270    0.6801679133
#> 721_Dlg1               0.6473547      0.95912280    0.8394579780
#> 722_N4bp1              0.7528198      0.57340016    0.5509669344
#> 723_Nudc               0.9996591      0.84503288    0.8744753671
#> 724_Znf367             1.0000000      0.83022221    1.0000000000
#> 725_Ring1              0.9307018      0.83239038    0.8187916031
#> 726_Snrpa1             0.7976077      0.92922311    0.9174607580
#> 727_U2surp             0.5975031      0.75726680    0.8620480770
#> 728_LOC100764225       0.7976077      0.78620962    0.9460896228
#> 729_Cdc42ep1           0.9996591      0.76519496    0.9446298214
#> 730_Znf385a            0.9150685      0.58893964    0.2441304415
#> 731_Ints1              0.9150685      0.58893964    0.2441304415
#> 732_LOC113833392       0.9150685      0.58893964    0.2441304415
#> 733_Lrch4              1.0000000      0.70094781    0.8515168795
#> 734_Ctdspl2            0.9996591      0.83239038    0.2384727084
#> 735_Prpf4b             0.8540623      0.83239038    0.5457472072
#> 736_Luzp1              0.7129418      0.84244807    0.9797810656
#> 737_Eif3b              0.1475220      0.77220390    0.3458807846
#> 738_Ptpn14             0.7528198      0.85102646    0.9701805996
#> 739_Rrp1               0.9937178      0.97808534    0.5959542038
#> 740_Lrrfip2            0.3045423      0.39726695    0.2738360002
#> 741_Nsfl1c             0.8913195      0.73710357    0.4463656788
#> 742_Ddx51              0.9996591      0.96609855    0.8474656873
#> 743_Prpf38b            0.5629666      0.52179637    0.2738360002
#> 744_Eef1b2             0.9307018      0.97469270    0.7962693514
#> 745_Znf385a            1.0000000      0.70094781    0.4314319337
#> 746_Map9               0.7239466      0.84674366    0.3669729209
#> 747_Rflnb              0.7239466      0.84674366    0.3669729209
#> 748_NA                 0.4884807      0.98226136    0.5103974918
#> 749_C1H12orf45         0.4884807      0.98226136    0.5103974918
#> 750_U2surp             0.7528198      0.91266136    0.3807330765
#> 751_Caskin2            0.8443625      0.97469270    0.5103974918
#> 752_Eri1               0.9173402      0.90047417    0.4358371316
#> 753_Gsk3b              0.9206541      0.98261381    0.6647202777
#> 754_LOC100766946       0.8789016      0.91266136    0.9240238739
#> 755_Cnpy3              1.0000000      1.00000000    0.2729399625
#> 756_Hnrnpc             0.7391776      0.83239038    1.0000000000
#> 757_Ptpn14             0.7391776      0.83239038    1.0000000000
#> 758_Slc7a11            0.6968896      0.90756350    0.5619464051
#> 759_Hnrnpc             0.7461617      0.95912280    0.5693079184
#> 760_Cdc37l1            0.4116137      0.32541914    0.3587160453
#> 761_LOC100768405       0.4884807      0.63814186    0.8311625987
#> 762_Rragc              0.4884807      0.63814186    0.8311625987
#> 763_LOC113834282       0.9996591      0.69493031    0.4028573087
#> 764_Fus                1.0000000      1.00000000    0.6629729761
#> 765_Ubxn1              0.4884807      0.84244807    0.9906710175
#> 766_Mmut               0.4884807      0.95236696    0.3911676182
#> 767_Pdcd11             0.4470790      0.70094781    0.6897028754
#> 768_LOC100757535       0.8353686      0.10784466    0.1406775311
#> 769_Eif3b              0.8039534      0.83239038    0.2895665722
#> 770_Rnf113a            0.9996591      0.89309712    0.3770190551
#> 771_Sytl4              0.9996591      0.58893964    0.3183586939
#> 772_Tlnrd1             0.8608189      0.96609855    0.9456073392
#> 773_H671_1g1131        0.9173402      0.74796187    0.5457472072
#> 774_Neurl1             0.7619059      0.61840397    0.3449452914
#> 775_Zyx                1.0000000      1.00000000    1.0000000000
#> 776_Ctdspl2            0.7619059      0.76519496    0.6272688737
#> 777_Chaf1b             0.9996591      0.57340016    0.1163090667
#> 778_Rragc              0.9937178      0.96609855    0.6724947438
#> 779_Srfbp1             0.4335248      0.51929656    0.1896719181
#> 780_Gys1               0.7585314      0.81647594    0.9292673283
#> 781_Usp15              0.6473547      0.99844582    0.5521486360
#> 782_Arhgef40           0.6485597      0.58893964    0.8120905602
#> 783_Gigyf1             1.0000000      1.00000000    1.0000000000
#> 784_Minar1             0.8540623      0.80301683    0.3458807846
#> 785_Dus2               0.9996591      0.57340016    0.4927937173
#> 786_Gatad2b            0.9173402      0.41853713    0.1662906723
#> 787_Eif5               1.0000000      0.27486659    0.1163090667
#> 788_Epb41l2            1.0000000      1.00000000    1.0000000000
#> 789_Arl6ip4            0.9173402      0.94517210    0.7759084747
#> 790_Plin4              0.9361973      0.91767997    0.9591634617
#> 791_Elf2               0.9996591      0.90654204    0.9090160656
#> 792_Plin4              1.0000000      0.57340016    0.3183586939
#> 793_Snip1              0.6933908      0.83239038    0.5423506258
#> 794_Txlng              0.5731500      0.84244807    0.7530007652
#> 795_LOC100769437       0.9307018      0.34703515    0.1163090667
#> 796_Caskin2            0.9996591      0.67796309    0.3183586939
#> 797_NA                 0.7987719      0.97462156    0.3557894948
#> 798_Synm               0.9428566      0.97469270    0.6803315820
#> 799_Synm               0.9206541      0.98168118    0.9690539503
#> 800_Ube2c              0.8447753      1.00000000    1.0000000000
#> 801_Sgtb               0.7528198      0.57340016    0.6275469607
#> 802_Prpf4b             0.5601500      0.84764375    0.5457472072
#> 803_Epb41l2            0.8241813      0.49257001    0.8453895336
#> 804_Mllt1              0.5975031      0.68895932    0.5778503376
#> 805_LOC100759640       1.0000000      1.00000000    1.0000000000
#> 806_Epb41l2            0.8478867      0.78620962    0.3814379204
#> 807_Znf280b            0.9996591      0.90654204    0.8473517850
#> 808_Kiaa1143           1.0000000      1.00000000    1.0000000000
#> 809_Gas2l1             0.5975031      0.32541914    0.8474656873
#> 810_Srp72              0.4884807      0.54140969    0.6460689701
#> 811_Tomm22             0.5975031      0.92922311    0.3498488866
#> 812_Psip1              0.4884807      0.97469270    0.3548546154
#> 813_Arhgef37           0.8478867      0.89453090    0.2796179498
#> 814_Bckdk              0.9615458      0.77660746    0.4927937173
#> 815_Strip1             0.8340945      0.74703751    0.1975058450
#> 816_Usp15              0.9173402      0.57955504    0.2738360002
#> 817_Ssr3               0.9996591      0.63814186    0.3243591041
#> 818_Strip1             0.7619059      0.70094781    0.3154653432
#> 819_Eif3b              0.4884807      0.98261381    0.1761106619
#> 820_U2surp             0.5682360      0.84518166    0.2102464538
#> 821_Bend3              0.9485913      0.57340016    0.2441304415
#> 822_Rps10              1.0000000      1.00000000    1.0000000000
#> 823_Rpl23a             0.8340945      0.34703515    0.1395399686
#> 824_Nbn                0.7992518      0.57340016    0.3879479941
#> 825_Rpap3              1.0000000      0.88113882    0.9460896228
#> 826_LOC100759640       0.9937178      0.89979241    0.2928150105
#> 827_Ric8a              0.9240383      0.80301683    0.2738360002
#> 828_Hsph1              0.9937178      0.66216815    0.4188878049
#> 829_LOC100759640       0.9409124      0.84764375    0.9273415365
#> 830_LOC100757535       1.0000000      1.00000000    1.0000000000
#> 831_Gigyf1             0.6507661      0.57955504    0.1272404753
#> 832_Dbn1               1.0000000      1.00000000    1.0000000000
#> 833_Snrk               1.0000000      1.00000000    1.0000000000
#> 834_Prpf38b            0.7224935      0.59567667    0.3458807846
#> 835_LOC100766868       1.0000000      1.00000000    1.0000000000
#> 836_LOC100766868       0.9206541      0.96609855    0.1732450585
#> 837_Wbp11              0.7671830      0.95236696    0.5619464051
#> 838_Rusc2              0.7976077      0.91629651    0.5103974918
#> 839_Eif3b              0.8711111      0.80301683    0.7764469293
#> 840_Ptpn14             0.9996591      0.89148887    0.6928203220
#> 841_Rlim               0.9996591      0.65495110    0.6654892405
#> 842_Ints1              0.6473547      0.84244807    0.6359380798
#> 843_Chaf1b             0.9206541      0.97469270    0.2496049280
#> 844_Dlg1               0.7239466      1.00000000    1.0000000000
#> 845_Lamtor1            0.7619059      0.76519496    0.2738360002
#> 846_Tab1               0.7671830      0.72336026    0.3558907560
#> 847_Dbn1               0.9362893      0.44694320    0.1896719181
#> 848_Psip1              0.9362893      0.44694320    0.1896719181
#> 849_Dbn1               0.9362893      0.84244807    0.4945708006
#> 850_Pabpc1             0.9524993      0.90047417    0.6987798252
#> 851_Hnrnpc             0.9361973      0.61641970    0.6595569957
#> 852_Emd                0.8049314      0.76519496    0.5319516213
#> 853_LOC100764225       0.8540623      0.99844582    0.8085069881
#> 854_Nup50              0.9442833      0.97469270    0.6619954261
#> 855_Ctcf               0.9173402      0.89883798    0.8444428354
#> 856_Raly               0.9307018      0.99844582    0.5103974918
#> 857_Bard1              0.9996591      0.84138511    0.6595569957
#> 858_Ptpn14             0.8540623      0.91548160    0.6262551451
#> 859_LOC100757535       0.5536475      0.62588442    0.1691476201
#> 860_Psmd2              1.0000000      1.00000000    1.0000000000
#> 861_Junb               0.9937178      0.81938063    0.3444986722
#> 862_C1qbp              0.8241813      0.73046398    0.5696238511
#> 863_Lrch4              0.8327078      0.89453090    0.6353934056
#> 864_CUNH14orf93        1.0000000      1.00000000    1.0000000000
#> 865_U2surp             0.9996224      0.63315255    0.3879479941
#> 866_Raly               1.0000000      1.00000000    1.0000000000
#> 867_LOC100774417       0.5842771      0.57340016    0.6803315820
#> 868_Srp72              0.9085372      1.00000000    1.0000000000
#> 869_LOC100764225       0.9409124      0.91266136    1.0000000000
#> 870_Morf4l2            0.7976077      0.58893964    0.3669755656
#> 871_CUNH9orf40         0.5975031      0.41581591    0.2441304415
#> 872_Gas2l1             0.9996591      0.70194592    0.8179717078
#> 873_Atp5pf             0.9996591      0.73152943    0.7126478180
#> 874_Lrrfip2            0.7976077      0.59567667    0.7962693514
#> 875_Prpf4b             1.0000000      0.84244807    0.2738360002
#> 876_Top2b              0.9996591      0.75726680    0.1542353843
#> 877_Mepce              0.4335248      0.75726680    0.8085069881
#> 878_Ptpn14             0.5731500      0.89453090    0.8474656873
#> 879_Dnajc25            0.6473547      0.90654204    0.1492654263
#> 880_Cbx8               0.9361973      0.57340016    0.2738360002
#> 881_Synm               0.6982773      0.96609855    0.5103974918
#> 882_Def6               1.0000000      1.00000000    1.0000000000
#> 883_Gys1               0.9996591      0.92123506    0.9773291621
#> 884_Luzp1              0.9996591      0.84764375    0.3437425283
#> 885_Synm               0.8241813      0.65688382    0.4655348841
#> 886_Snip1              0.7239466      0.39726695    0.2738360002
#> 887_Top2b              0.5601500      0.75726680    0.9273415365
#> 888_NA                 0.8566323      0.84518166    1.0000000000
#> 889_Trim35             1.0000000      1.00000000    1.0000000000
#> 890_Znf385a            0.9996591      0.73046398    0.4209451907
#> 891_Chaf1b             0.6971286      0.84244807    0.6801679133
#> 892_Abcf1              1.0000000      1.00000000    1.0000000000
#> 893_Pdcd11             0.5601500      0.68895932    0.8744753671
#> 894_Dlg1               0.9361973      0.68895932    0.4021681980
#> 895_Dbn1               0.9996224      0.90047417    0.8041289506
#> 896_LOC100752363       0.9173402      0.89453090    0.9562143923
#> 897_Ppp4r3a            0.9996591      0.89453090    0.3955081305
#> 898_Gas2l1             0.7239466      0.57340016    0.4945708006
#> 899_Mtmr10             1.0000000      1.00000000    1.0000000000
#> 900_Cyld               0.7528198      0.97469270    0.9642694022
#> 901_NA                 0.9937178      0.63814186    0.8138575422
#> 902_Rnf113a            0.9996591      0.99844582    0.4945708006
#> 903_Nelfa              1.0000000      0.99844582    0.7068411866
#> 904_Zkscan1            0.9996591      0.99844582    0.7448404619
#> 905_Chaf1b             0.9206541      0.77660746    0.2441304415
#> 906_Eif3b              0.7585314      0.96609855    0.2796179498
#> 907_Top2b              0.8312584      0.77660746    0.4823171066
#> 908_Chaf1b             1.0000000      1.00000000    0.3458807846
#> 909_Epb41l2            1.0000000      1.00000000    1.0000000000
#> 910_C3H11orf58         1.0000000      1.00000000    0.3154653432
#> 911_Top2b              0.8540623      0.89228318    0.4987772572
#> 912_Wee1               0.5975031      0.74590123    0.9413836018
#> 913_Raly               0.9173402      0.57340016    0.4184815759
#> 914_H671_1g2680        0.9996591      0.89148887    0.6926345811
#> 915_Eef1b2             0.9085372      0.72336026    0.5103974918
#> 916_Gas2l1             0.7619059      0.47674732    0.1163090667
#> 917_Epb41l2            0.8540623      0.40322706    0.1395399686
#> 918_Rpl23a             0.8540623      0.84503288    0.3342221489
#> 919_Chmp2b             1.0000000      0.34703515    0.3669729209
#> 920_Lrrfip2            0.7585314      0.99844582    0.8611757851
#> 921_Aldoa              0.9634942      0.70868902    0.2090589517
#> 922_Cby1               0.4884807      0.49995158    0.5805572519
#> 923_LOC100759640       0.9426037      0.57340016    0.4187969375
#> 924_Rbm28              0.3045423      0.32541914    0.8053523228
#> 925_Skiv2l             1.0000000      1.00000000    0.1761106619
#> 926_Ints1              0.9996591      0.89797917    0.8474656873
#> 927_Ehd1               0.8617286      0.39726695    0.4889794954
#> 928_Nr2f6              0.8093895      0.74503582    0.8053523228
#> 929_Top2b              0.9206541      0.83239038    0.6803315820
#> 930_Lrrfip2            0.9937178      0.86446114    0.6320453800
#> 931_Pip4p2             0.9996591      0.76519496    0.3458807846
#> 932_Srp72              0.7789744      0.95648901    0.7448404619
#> 933_Mtmr9              0.9996591      0.57340016    0.1552486432
#> 934_Gigyf1             0.9996591      0.80301683    0.7347336557
#> 935_Rbm7               0.9173402      0.62588442    0.2738360002
#> 936_LOC100773565       0.8156605      0.77660746    0.4919700587
#> 937_Trim35             0.8327078      0.59876431    0.3676390944
#> 938_Cbx8               0.9996591      0.92123506    0.4830442104
#> 939_Rplp0              0.9918258      1.00000000    1.0000000000
#> 940_Aldoa              0.9996591      0.97469270    0.9460896228
#> 941_NA                 0.6213394      0.50159122    0.3277644801
#> 942_Zyx                1.0000000      0.94208942    0.7566755605
#> 943_Psip1              0.5731500      0.57340016    0.5532033310
#> 944_Slc7a11            0.5682360      0.76519496    0.8179717078
#> 945_Miga2              0.5805890      0.39726695    0.1552486432
#> 946_Arhgef6            0.7886993      0.96609855    0.2988269514
#> 947_Dlgap4             0.9937178      0.73152943    0.8873517778
#> 948_Ampd2              1.0000000      1.00000000    1.0000000000
#> 949_Luzp1              0.5445966      0.80301683    0.9174607580
#> 950_Camlg              1.0000000      1.00000000    1.0000000000
#> 951_Pfkfb3             0.8705424      0.85102646    0.3458807846
#> 952_NA                 1.0000000      1.00000000    1.0000000000
#> 953_Raly               1.0000000      1.00000000    1.0000000000
#> 954_Kiaa1143           0.5731500      0.41581591    0.1395399686
#> 955_Bcar1              0.9937178      0.62588442    0.4014644565
#> 956_Gatad2b            0.4884807      0.57955504    0.5319516213
#> 957_Eif4ebp2           0.3045423      0.32541914    0.4358371316
#> 958_Fam76b             0.4884807      0.57340016    0.3458807846
#> 959_Camlg              0.8539081      0.67937672    0.3770190551
#> 960_LOC100754077       0.9996591      0.63814186    0.1343299522
#> 961_NA                 0.9996591      0.29247726    0.0237123566
#> 962_Epb41l2            0.8754165      0.99480848    0.4945708006
#> 963_Ankrd34a           1.0000000      1.00000000    1.0000000000
#> 964_Zc3h15             0.9996591      0.78620962    0.2610754681
#> 965_Def6               0.4335248      0.92434829    0.4034772289
#> 966_Srsf6              1.0000000      1.00000000    0.5502559641
#> 967_H671_4g11480       0.5975031      0.74796187    0.7613779935
#> 968_Top2b              1.0000000      1.00000000    0.9464121243
#> 969_LOC100769471       0.9996591      0.83239038    0.4889794954
#> 970_Raver1             0.4884807      0.67222253    0.4987772572
#> 971_Etv3               1.0000000      0.96707635    1.0000000000
#> 972_Psd                0.9996591      0.89797917    0.8042667469
#> 973_Usp15              1.0000000      1.00000000    0.7211170537
#> 974_Nol7               0.9159876      0.74938527    0.8085069881
#> 975_Stk38              0.9996591      0.58893964    0.1691476201
#> 976_Smim13             0.4937630      0.76934373    0.3812768780
#> 977_Etv3               0.4884807      0.10784466    0.1374908882
#> 978_Synm               0.4884807      0.39726695    0.6449466883
#> 979_Pwp1               0.9996591      0.75726680    0.4111144282
#> 980_Fus                0.9150685      0.89309712    0.9520677088
#> 981_Junb               0.8739042      0.90581890    0.7589135521
#> 982_Phf8               0.7482648      0.70094781    0.9164066845
#> 983_Nelfa              1.0000000      1.00000000    1.0000000000
#> 984_Prpf4b             1.0000000      1.00000000    1.0000000000
#> 985_Abraxas1           0.4470790      0.84244807    0.6840169429
#> 986_Prpf4b             0.9996591      1.00000000    1.0000000000
#> 987_Raver1             1.0000000      1.00000000    1.0000000000
#> 988_Caap1              0.9159876      0.95648901    0.9016849838
#> 989_Rpap3              0.9937178      0.94475915    0.9174607580
#> 990_Hsph1              0.9937178      0.74796187    0.4028573087
#> 991_LOC100750437       1.0000000      1.00000000    1.0000000000
#> 992_Mepce              0.8156605      0.98300959    0.6582681907
#> 993_Efs                0.9937178      0.68895932    0.8591498257
#> 994_Epb41l2            0.9996591      0.90756350    0.9628343032
#> 995_Abcf1              0.4884807      0.50159122    0.4823171066
#> 996_NA                 0.9085372      0.96609855    0.6724947438
#> 997_Eif4ebp2           1.0000000      1.00000000    1.0000000000
#> 998_Pfkfb3             0.6295160      0.91266136    0.8120905602
#> 999_Hnrnpc             0.7224935      0.44694320    0.3587160453
#> 1000_Psmd2             0.9307018      0.50159122    0.4967507613
#>                  120_vs_neighbors
#> 1_Top2b                0.99332811
#> 2_NA                   0.65273244
#> 3_Snip1                1.00000000
#> 4_Tomm34               0.86859820
#> 5_Pus3                 0.96317292
#> 6_Ints1                1.00000000
#> 7_Mlh1                 0.95822477
#> 8_LOC100750437         0.78837220
#> 9_Pabpc1               0.81550950
#> 10_Top2b               0.87684327
#> 11_Gorasp1             0.57624609
#> 12_Ints1               1.00000000
#> 13_Syvn1               0.98272300
#> 14_Znf280b             0.71378983
#> 15_Mrnip               0.87125683
#> 16_Rragc               0.87684327
#> 17_Gorasp1             0.74897796
#> 18_Tomm34              1.00000000
#> 19_LOC100757430        0.86650073
#> 20_Ubxn1               0.89691748
#> 21_H671_1g1131         1.00000000
#> 22_Luzp1               0.84433926
#> 23_Efs                 0.96309797
#> 24_Mta2                0.59299794
#> 25_Nedd1               0.53998241
#> 26_Gigyf1              0.95401692
#> 27_Myh9                0.96317292
#> 28_Caskin2             0.61804675
#> 29_Papolg              0.53998241
#> 30_Tfg                 1.00000000
#> 31_Rpl34               0.89691748
#> 32_Mideas              0.53647630
#> 33_Gys1                0.53998241
#> 34_Arhgef6             0.53998241
#> 35_Ctdspl2             0.99040703
#> 36_Ptpn14              0.76552820
#> 37_Raly                1.00000000
#> 38_Znhit3              0.82285588
#> 39_LOC113833392        1.00000000
#> 40_Luc7l3              0.99040703
#> 41_Rplp0               1.00000000
#> 42_Gys1                0.65792365
#> 43_Rpl22l1             0.61804675
#> 44_Eif3b               0.80741720
#> 45_Med26               0.71891081
#> 46_Mepce               0.71378983
#> 47_Pdcd11              0.97937571
#> 48_Twf1                1.00000000
#> 49_LOC100759640        1.00000000
#> 50_Wrnip1              0.76463145
#> 51_Poldip3             0.89691748
#> 52_Ampd2               0.64034789
#> 53_Mea1                0.81456911
#> 54_Dbn1                0.89691748
#> 55_Snip1               0.61127133
#> 56_Srsf6               0.84680256
#> 57_LOC113834282        0.11563729
#> 58_Map9                0.69586387
#> 59_Cdc42ep1            0.53998241
#> 60_Poldip3             0.96317292
#> 61_LOC100764225        0.54836325
#> 62_Epb41l2             1.00000000
#> 63_H671_4g11480        0.76463145
#> 64_Nbn                 1.00000000
#> 65_U2surp              1.00000000
#> 66_Gigyf1              0.59107576
#> 67_NA                  0.57624609
#> 68_Luc7l3              0.57624609
#> 69_LOC100752363        0.86650073
#> 70_Ampd2               0.53647630
#> 71_LOC100759640        0.53647630
#> 72_Stam                0.89691748
#> 73_Nsfl1c              0.71891081
#> 74_Pfkfb3              0.96317292
#> 75_Rad23a              0.61804675
#> 76_Elf2                0.97937571
#> 77_Crem                0.87125683
#> 78_Rragc               0.82801076
#> 79_Lrrfip2             0.89691748
#> 80_Zyx                 0.46858598
#> 81_Lrrfip2             1.00000000
#> 82_Gatad2b             0.75084189
#> 83_Bcar1               0.82531804
#> 84_Ehd1                0.61127133
#> 85_LOC113834282        0.71891081
#> 86_Tmem230             0.76552820
#> 87_Ncbp1               0.89691748
#> 88_Mllt1               1.00000000
#> 89_Stk17b              0.53647630
#> 90_Dlgap4              0.76552820
#> 91_Papolg              0.67058556
#> 92_Cyld                0.61127133
#> 93_Gigyf1              0.99790083
#> 94_Lrrfip2             0.95365528
#> 95_Lrrfip2             0.76463145
#> 96_Rlim                0.89977838
#> 97_Eif3b               0.95365528
#> 98_Mphosph10           0.65273244
#> 99_Gatad2b             0.96317292
#> 100_Srsf6              0.58004623
#> 101_Zyx                0.80007140
#> 102_Mphosph10          0.59107576
#> 103_Psip1              0.89691748
#> 104_Fbl                1.00000000
#> 105_H671_1g2680        0.98272300
#> 106_Sgtb               0.99040703
#> 107_Gnl3               0.96317292
#> 108_Eif3b              0.69265742
#> 109_Serpinb1           0.61127133
#> 110_N4bp1              0.79370312
#> 111_Snip1              0.82801076
#> 112_Psip1              0.86686270
#> 113_Mlh1               0.87684327
#> 114_Bsg                0.80519954
#> 115_Tnpo1              1.00000000
#> 116_H671_1g2680        0.53998241
#> 117_Cbx8               0.71378983
#> 118_Mideas             0.53998241
#> 119_Mideas             0.96317292
#> 120_Dcun1d3            0.96317292
#> 121_Dlg1               0.96317292
#> 122_Rad23a             0.89691748
#> 123_Srsf6              0.87125683
#> 124_Stx7               0.84433926
#> 125_Pdcd11             0.92281298
#> 126_Kiaa1958           0.55994835
#> 127_Pwp1               0.67058556
#> 128_Txlng              0.76552820
#> 129_Junb               0.86650073
#> 130_LOC100759640       0.89691748
#> 131_Dbn1               0.61127133
#> 132_Top2b              0.59107576
#> 133_Rusc2              0.80007140
#> 134_NA                 1.00000000
#> 135_LOC113837251       0.97937571
#> 136_Fam76b             0.98272300
#> 137_Ptpn14             0.76463145
#> 138_Chmp4b             0.98272300
#> 139_Prpf4b             0.61127133
#> 140_Eif3b              1.00000000
#> 141_Nsfl1c             0.87684327
#> 142_Pdlim7             0.95731924
#> 143_Rnf113a            0.80822413
#> 144_Epb41l2            0.71891081
#> 145_Hnrnpc             0.75084189
#> 146_LOC113834282       0.89691748
#> 147_Plekho2            0.86319749
#> 148_Med26              0.98272300
#> 149_Arhgef40           0.95013062
#> 150_NA                 0.53647630
#> 151_Phf8               0.78837220
#> 152_Minar1             1.00000000
#> 153_H671_21690         0.87435458
#> 154_Arhgef40           0.70198358
#> 155_Chaf1b             0.84519728
#> 156_Prpf4b             0.89691748
#> 157_Znf367             0.65273244
#> 158_Luzp1              0.87435458
#> 159_LOC113833882       0.81550950
#> 160_Hnrnpc             0.61127133
#> 161_Mepce              0.76552820
#> 162_Ubxn1              0.96317292
#> 163_Mllt1              0.71891081
#> 164_Chaf1b             1.00000000
#> 165_Raly               0.89691748
#> 166_Gas2l1             0.74161192
#> 167_Dlg1               0.46858598
#> 168_Hoxc10             1.00000000
#> 169_Gigyf1             1.00000000
#> 170_Luzp1              0.89691748
#> 171_Srp72              0.81284520
#> 172_LOC100771461       1.00000000
#> 173_Chaf1b             0.65273244
#> 174_C3H11orf58         0.73477355
#> 175_Pdcd11             0.99040703
#> 176_Psip1              0.81456911
#> 177_Prpf4b             0.53998241
#> 178_Rnf113a            0.61127133
#> 179_Irf3               0.87684327
#> 180_Smim13             0.89691748
#> 181_Gnl3               0.89691748
#> 182_Psma5              0.55097037
#> 183_Ptpn14             0.97937571
#> 184_Prpf4b             0.78837220
#> 185_Top2b              0.91498011
#> 186_Prpf38b            0.89691748
#> 187_Epb41l2            1.00000000
#> 188_Eif3b              0.95731924
#> 189_Hnrnpc             0.67058556
#> 190_LOC100758278       0.96425465
#> 191_Prpf4b             0.76463145
#> 192_Caskin2            0.31331155
#> 193_LOC100752363       0.96317292
#> 194_Septin6            0.61127133
#> 195_Max                0.59299794
#> 196_Mid1ip1            0.88019189
#> 197_NA                 0.76552820
#> 198_Hsph1              0.67058556
#> 199_Nol7               0.87684327
#> 200_Raly               0.96317292
#> 201_Smim13             0.85184664
#> 202_LOC100757535       0.87684327
#> 203_Net1               0.95831218
#> 204_LOC100754077       0.99821395
#> 205_Snip1              0.76552820
#> 206_Hnrnpc             0.01545919
#> 207_Ldlrap1            0.89691748
#> 208_Luzp1              0.61127133
#> 209_Rpl26              0.87684327
#> 210_Epb41l2            0.78427438
#> 211_Znf367             0.54346611
#> 212_Dlgap4             0.61127133
#> 213_Plekho2            1.00000000
#> 214_Zpr1               0.72214605
#> 215_Dlgap4             1.00000000
#> 216_Def6               0.76552820
#> 217_Eif4ebp2           0.61127133
#> 218_Eef1b2             0.89691748
#> 219_Rad23a             0.80007140
#> 220_Morf4l2            0.61127133
#> 221_Arhgef40           0.82069736
#> 222_NA                 0.87684327
#> 223_LOC100773565       0.81550950
#> 224_Dus2               0.95897789
#> 225_Pip4p2             0.92872331
#> 226_Top2b              0.87684327
#> 227_Znf280b            0.76463145
#> 228_Pdcd11             0.71378983
#> 229_Bckdk              0.76552820
#> 230_Arhgef40           0.53647630
#> 231_Mepce              0.53998241
#> 232_Ccnd3              0.89691748
#> 233_Phf8               0.87684327
#> 234_H671_1g2680        1.00000000
#> 235_Ell                0.96317292
#> 236_U2surp             0.53998241
#> 237_Rps10              0.86686270
#> 238_Ctdspl2            0.46858598
#> 239_Top2b              0.61947174
#> 240_Msantd3            0.87684327
#> 241_Fam76b             0.96317292
#> 242_Ppp4r3a            0.64518630
#> 243_Gpatch4            1.00000000
#> 244_Nudc               0.99040703
#> 245_Nol7               0.96425465
#> 246_Plekho2            0.86650073
#> 247_Prpf4b             0.53647630
#> 248_Mta2               1.00000000
#> 249_U2surp             0.79131251
#> 250_Ubxn1              0.71891081
#> 251_Rlim               0.97937571
#> 252_Atat1              0.76552820
#> 253_Ubxn1              0.96317292
#> 254_H671_1g2680        0.76552820
#> 255_eIF2aK2            0.68492351
#> 256_Skiv2l             0.71378983
#> 257_Rpl28              0.82260728
#> 258_LOC100759640       0.53998241
#> 259_Gatad2b            0.80822413
#> 260_NA                 0.89691748
#> 261_Gprasp1            0.82285588
#> 262_Luzp1              0.88372742
#> 263_Slc1a5             0.81550950
#> 264_LOC113834282       0.71378983
#> 265_Srsf6              0.98272300
#> 266_Cdc42ep1           0.53998241
#> 267_Net1               1.00000000
#> 268_Caskin2            0.55406551
#> 269_LOC100759640       0.82243218
#> 270_Mideas             0.91879652
#> 271_Luzp1              1.00000000
#> 272_Emd                0.91734510
#> 273_Plpp6              0.81907261
#> 274_LOC100759640       0.74897796
#> 275_Rps7               0.74897796
#> 276_Fkbp1a             0.97937571
#> 277_Gatad2b            0.53998241
#> 278_Znf385a            0.65792365
#> 279_Arhgef6            0.84433926
#> 280_Slirp              0.85184664
#> 281_Skiv2l             0.71378983
#> 282_H671_21690         0.89691748
#> 283_Kat8               0.95822477
#> 284_Nkap               0.89691748
#> 285_Gsk3b              0.95731924
#> 286_Ints1              0.74897796
#> 287_Gas2l1             1.00000000
#> 288_LOC100759640       0.81456911
#> 289_Top2b              0.65273244
#> 290_Kif20b             0.87684327
#> 291_Phf8               1.00000000
#> 292_Snip1              0.61127133
#> 293_Gsk3b              0.96317292
#> 294_Caskin2            0.96317292
#> 295_C3H11orf58         1.00000000
#> 296_Lrch4              0.85184664
#> 297_LOC113834282       0.80822413
#> 298_LOC100750407       0.87125683
#> 299_LOC113833392       0.89691748
#> 300_LOC113833882       0.76463145
#> 301_Ldlrap1            0.87684327
#> 302_Wee1               0.74897796
#> 303_Caap1              1.00000000
#> 304_Eif4ebp2           0.81284520
#> 305_Ripk2              0.86563801
#> 306_Srp72              0.59299794
#> 307_Taok2              1.00000000
#> 308_Nr2f6              0.82285588
#> 309_Arhgef40           0.93783336
#> 310_Gys1               0.76552820
#> 311_Dlg1               0.86650073
#> 312_Vapb               0.78483347
#> 313_LOC100757535       1.00000000
#> 314_Mkrn2              0.99040703
#> 315_Eif3b              0.81456911
#> 316_Isyna1             0.61804675
#> 317_Prpf4b             0.81912718
#> 318_LOC113833882       0.98272300
#> 319_Lrch4              0.96425465
#> 320_Dbn1               0.93345068
#> 321_Abcf1              0.71891081
#> 322_Ints1              0.76463145
#> 323_C3H11orf58         0.78837220
#> 324_Psma5              1.00000000
#> 325_Fundc1             0.98272300
#> 326_Papolg             0.78837220
#> 327_Mideas             0.61198963
#> 328_Ubxn1              0.89691748
#> 329_Synm               1.00000000
#> 330_Arhgef6            1.00000000
#> 331_Ptpn14             1.00000000
#> 332_Pgrmc1             1.00000000
#> 333_Myh9               0.89691748
#> 334_Etv3               0.81550950
#> 335_Ip6k1              0.59299794
#> 336_Luzp1              0.96425465
#> 337_Ptpn14             0.78837220
#> 338_Caskin2            0.96317292
#> 339_Chaf1b             0.61127133
#> 340_Ubxn1              1.00000000
#> 341_Ube2c              0.57624609
#> 342_Gins2              0.98272300
#> 343_Nlgn2              0.71425490
#> 344_Nf2                0.76463145
#> 345_Pip4p2             0.90249140
#> 346_Emd                0.89691748
#> 347_Top2b              0.85184664
#> 348_Trim35             0.99790083
#> 349_NA                 0.81550950
#> 350_NA                 1.00000000
#> 351_Mideas             0.61804675
#> 352_Gas2l1             0.96317292
#> 353_Ampd2              0.99040703
#> 354_Calu               0.85779284
#> 355_Fam76b             0.86650073
#> 356_Dlg1               1.00000000
#> 357_Srsf6              0.89691748
#> 358_Chaf1b             0.89691748
#> 359_Dbn1               0.96302284
#> 360_Tcf25              0.61804675
#> 361_Psip1              0.73733493
#> 362_Cnpy3              0.80822413
#> 363_LOC100759640       0.83332028
#> 364_Zyx                1.00000000
#> 365_Lrch4              0.86650073
#> 366_Bola1              1.00000000
#> 367_Znf385a            0.89691748
#> 368_Kif20b             0.87435458
#> 369_Ell                1.00000000
#> 370_Ell                1.00000000
#> 371_Srsf6              0.61127133
#> 372_Pwp1               0.87684327
#> 373_Def6               0.98272300
#> 374_Cbx8               1.00000000
#> 375_Ddx51              0.75084189
#> 376_Psip1              0.89691748
#> 377_Arhgef40           0.96317292
#> 378_Raly               0.53647630
#> 379_NA                 0.99040703
#> 380_Lrrfip2            0.70208246
#> 381_Gnl3               0.98272300
#> 382_Caskin2            0.81456911
#> 383_Rragc              0.89691748
#> 384_Caskin2            0.61127133
#> 385_Bcar1              0.98272300
#> 386_Homer3             0.95731924
#> 387_Luzp1              0.96317292
#> 388_N4bp1              1.00000000
#> 389_Ppp4r3a            0.91526116
#> 390_H671_1g2680        1.00000000
#> 391_Gnl3               0.61127133
#> 392_Top2b              0.71425490
#> 393_Oser1              0.98272300
#> 394_Snrk               0.57624609
#> 395_Kat8               1.00000000
#> 396_Raver1             1.00000000
#> 397_Pdcd11             0.95365528
#> 398_Rps20              0.76552820
#> 399_Bsg                1.00000000
#> 400_Raly               0.53998241
#> 401_Pdcd2              0.65273244
#> 402_Caskin2            0.61127133
#> 403_LOC100773571       0.53647630
#> 404_Papolg             0.89691748
#> 405_LOC100757535       0.84139140
#> 406_Caap1              0.53998241
#> 407_Psip1              0.73477355
#> 408_Dbn1               0.81456911
#> 409_Mta2               0.96317292
#> 410_Abcf1              1.00000000
#> 411_LOC100754108       0.97937571
#> 412_Slirp              0.61127133
#> 413_Nelfa              0.98891432
#> 414_Aggf1              0.99040703
#> 415_Bap1               0.86859820
#> 416_Luc7l3             0.82069736
#> 417_Rrp1               0.88019189
#> 418_Wrnip1             1.00000000
#> 419_NA                 0.87684327
#> 420_Abcf1              0.77944157
#> 421_Cluap1             0.89691748
#> 422_Hnrnpc             0.53647630
#> 423_Ptpn1              0.46858598
#> 424_Myh9               0.89691748
#> 425_U2surp             0.76552820
#> 426_NA                 0.98272300
#> 427_Arhgef40           0.86859820
#> 428_Chaf1b             0.94766133
#> 429_Prpf4b             1.00000000
#> 430_Epb41l2            1.00000000
#> 431_Eif3b              0.76552820
#> 432_Isyna1             0.55097037
#> 433_U2surp             0.99040703
#> 434_LOC100765020       0.53998241
#> 435_Arhgef6            0.53647630
#> 436_Ptpn1              0.70202424
#> 437_Prpf4b             0.99040703
#> 438_Rpl35a             1.00000000
#> 439_Prpf4b             0.99404010
#> 440_Zyx                0.89691748
#> 441_Dbn1               0.87684327
#> 442_Chaf1b             0.74897796
#> 443_LOC113834282       0.76463145
#> 444_Gpsm2              0.81550950
#> 445_LOC100757535       0.82254097
#> 446_Cfap410            0.98891432
#> 447_Epb41l2            0.53998241
#> 448_Ncbp1              0.76463145
#> 449_Pacsin1            0.53647630
#> 450_Cstf2              0.87684327
#> 451_LOC100769437       0.84745243
#> 452_eIF2aK2            0.76552820
#> 453_Kiaa1191           0.99040703
#> 454_Mepce              0.61127133
#> 455_Cbx8               0.97937571
#> 456_Eed                0.53998241
#> 457_Cdc42ep1           0.89691748
#> 458_Lrrfip2            0.57624609
#> 459_Pacsin1            0.70453687
#> 460_Gpatch4            0.78427438
#> 461_Plin4              0.89691748
#> 462_NA                 0.95731924
#> 463_Snip1              0.43271982
#> 464_Cyld               0.61127133
#> 465_Plin4              0.98272300
#> 466_Twf1               1.00000000
#> 467_LOC113834282       1.00000000
#> 468_Snip1              0.87684327
#> 469_Ppp4r3a            1.00000000
#> 470_Psip1              0.71378983
#> 471_Dnajc5             0.89691748
#> 472_Phf8               0.46858598
#> 473_Bola1              0.87684327
#> 474_Cdc42ep1           0.86650073
#> 475_Eif4ebp2           0.80605441
#> 476_Prpf38b            0.89691748
#> 477_Klhl26             0.67058556
#> 478_Hsph1              0.96317292
#> 479_Snip1              0.80007140
#> 480_Caskin2            0.80519954
#> 481_Plpp6              0.98272300
#> 482_NA                 1.00000000
#> 483_Mlh1               0.96317292
#> 484_Gys1               0.99790083
#> 485_Tfg                0.92779889
#> 486_Arhgef6            0.65273244
#> 487_Mphosph10          0.65273244
#> 488_Hoxc10             0.83634429
#> 489_LOC100759640       0.67058556
#> 490_Arhgef40           0.99040703
#> 491_Dnajc5             0.98272300
#> 492_Tbc1d23            1.00000000
#> 493_Ubxn1              0.59422876
#> 494_Rab1a              0.89691748
#> 495_Eif3b              0.87684327
#> 496_Tceal8             0.80350549
#> 497_Dlgap4             0.78837220
#> 498_Smim13             0.53647630
#> 499_NA                 1.00000000
#> 500_Lrch4              0.46858598
#> 501_Bola1              0.59107576
#> 502_NA                 0.68324247
#> 503_Ptpn14             0.97937571
#> 504_LOC100759640       0.61127133
#> 505_Rps10              0.99040703
#> 506_Top2b              0.99040703
#> 507_Ssr3               0.61127133
#> 508_Homer3             0.61127133
#> 509_Phf8               0.87684327
#> 510_LOC100767716       0.86617797
#> 511_Xpa                0.53998241
#> 512_H671_21690         0.76463145
#> 513_LOC100769471       0.74897796
#> 514_Gas2l1             0.57624609
#> 515_Luzp1              0.53998241
#> 516_Gpbp1              0.57624609
#> 517_Gatad2b            1.00000000
#> 518_Gys1               0.89691748
#> 519_Top2b              0.89391002
#> 520_LOC100757535       1.00000000
#> 521_Lpcat4             0.98891432
#> 522_Arhgef6            0.65273244
#> 523_Cavin3             0.61127133
#> 524_Gpatch4            0.99333603
#> 525_Prpf38b            0.76463145
#> 526_Timm8a             0.53998241
#> 527_Cavin3             0.82069736
#> 528_Mkrn2              0.76552820
#> 529_Oser1              0.87684327
#> 530_Gsk3b              0.53998241
#> 531_Eef1b2             1.00000000
#> 532_Ampd2              0.97937571
#> 533_Lrrfip2            0.86650073
#> 534_Ring1              0.79131251
#> 535_Rlim               1.00000000
#> 536_LOC100759640       0.76552820
#> 537_LOC100759640       0.74897796
#> 538_Atp5pf             0.94766133
#> 539_Max                0.76552820
#> 540_Bap1               0.96309797
#> 541_Nsfl1c             0.53998241
#> 542_Prpf4b             0.89691748
#> 543_LOC100757535       0.87684327
#> 544_Mtmr10             0.66402033
#> 545_Hoxc10             1.00000000
#> 546_Trim35             0.61127133
#> 547_Eif4ebp2           0.92779889
#> 548_Dlgap4             1.00000000
#> 549_Gys1               0.86686270
#> 550_Sgtb               0.97937571
#> 551_Eri2               0.95365528
#> 552_Ccnd3              1.00000000
#> 553_Smim13             0.64380397
#> 554_Snrk               0.82069736
#> 555_Caskin2            0.70453687
#> 556_Pdcd11             1.00000000
#> 557_Pgam5              0.89691748
#> 558_Mphosph10          0.86650073
#> 559_Mideas             0.84433926
#> 560_Top2b              0.76463145
#> 561_LOC100763014       0.61117857
#> 562_Snip1              0.53998241
#> 563_Ubxn1              0.86859820
#> 564_LOC100750407       1.00000000
#> 565_Morf4l2            1.00000000
#> 566_Ctdspl2            0.89691748
#> 567_Cwf19l1            0.89691748
#> 568_Eef1b2             1.00000000
#> 569_C1H12orf45         1.00000000
#> 570_Znf367             0.53647630
#> 571_Ankrd34a           0.69265742
#> 572_Mllt11             0.99616460
#> 573_LOC100774792       0.88253951
#> 574_NA                 0.64270521
#> 575_Cbx8               0.86319749
#> 576_Bckdk              0.99040703
#> 577_Snip1              1.00000000
#> 578_Nsfl1c             1.00000000
#> 579_Gas2l1             0.65273244
#> 580_Nudc               0.46858598
#> 581_Epb41l2            1.00000000
#> 582_Mtmr6              0.74897796
#> 583_Znf668             0.95731924
#> 584_Hsph1              0.70198358
#> 585_LOC113834282       0.92669903
#> 586_Ctdspl2            1.00000000
#> 587_Foxf1              0.61127133
#> 588_Luzp1              0.82243218
#> 589_Xpa                0.88019189
#> 590_Psip1              0.57624609
#> 591_Rbm7               0.82148606
#> 592_Mtrex              0.91740623
#> 593_Arhgef40           1.00000000
#> 594_Plekho2            0.95731924
#> 595_Bckdk              0.76463145
#> 596_Dut                0.87684327
#> 597_Abcf1              0.87684327
#> 598_Txnl1              0.99500130
#> 599_Nudc               0.77944157
#> 600_Sh3gl1             0.99040703
#> 601_Gatad2b            0.82254097
#> 602_Homer3             0.89691748
#> 603_Septin6            0.89691748
#> 604_Smim13             0.85184664
#> 605_Arhgef40           0.81456911
#> 606_Rpl32              1.00000000
#> 607_Tomm34             0.96317292
#> 608_Mlh1               0.46858598
#> 609_Tbcc               1.00000000
#> 610_Eif3d              0.89691748
#> 611_Snrk               0.65273244
#> 612_Bckdk              0.97937571
#> 613_Wdr3               1.00000000
#> 614_LOC100757535       0.81153692
#> 615_Dlg1               0.99040703
#> 616_LOC100767716       0.70453687
#> 617_Hnrnpc             0.98272300
#> 618_Mphosph10          0.95365528
#> 619_Eif3b              0.95831218
#> 620_Emd                0.94326412
#> 621_Txlng              1.00000000
#> 622_Prpf4b             0.76552820
#> 623_Rlim               0.71891081
#> 624_Eef1b2             0.81907261
#> 625_Def6               0.71378983
#> 626_LOC100765020       0.61127133
#> 627_U2surp             0.81153692
#> 628_Elf2               1.00000000
#> 629_Slc1a5             0.99040703
#> 630_NA                 0.96317292
#> 631_Tfg                0.81153692
#> 632_Top2b              0.81550950
#> 633_Pip4p2             0.82203853
#> 634_Cdc42ep1           0.66219472
#> 635_Hsph1              1.00000000
#> 636_Twf1               1.00000000
#> 637_Nbn                0.95365528
#> 638_Psmd4              0.98272300
#> 639_Bap1               0.89691748
#> 640_Mepce              0.73741749
#> 641_Mideas             1.00000000
#> 642_LOC100759640       0.81550950
#> 643_Epb41l2            0.89691748
#> 644_Sav1               0.89691748
#> 645_Prpf4b             0.59107576
#> 646_Gnas               0.59107576
#> 647_Mllt1              0.67058556
#> 648_Poldip3            0.97937571
#> 649_Aldoa              0.61127133
#> 650_Rbbp8              0.98272300
#> 651_LOC113834282       0.76463145
#> 652_Gys1               0.86650073
#> 653_Hnrnpc             0.95365528
#> 654_Vps35              0.80822413
#> 655_Miga2              0.83144533
#> 656_Epb41l2            0.87684327
#> 657_Tob2               1.00000000
#> 658_Lamtor1            0.89691748
#> 659_LOC100759640       1.00000000
#> 660_Epb41l2            1.00000000
#> 661_Rlim               0.70202424
#> 662_Gys1               0.53647630
#> 663_LOC100750437       0.69265742
#> 664_NA                 0.46858598
#> 665_Nbn                0.89691748
#> 666_Tyw3               0.71891081
#> 667_Gas2l1             0.76463145
#> 668_Fus                1.00000000
#> 669_Prpf38b            0.71891081
#> 670_Calu               1.00000000
#> 671_Rras2              0.53647630
#> 672_Prpf4b             0.80519954
#> 673_Nelfa              0.98272300
#> 674_LOC100754077       1.00000000
#> 675_Rbm28              0.97937571
#> 676_Nsfl1c             0.87684327
#> 677_Rnf126             0.70134699
#> 678_Eme1               0.98272300
#> 679_Nbn                0.11563729
#> 680_Eif4ebp2           0.53647630
#> 681_Wee1               0.80741720
#> 682_Prpf38b            0.89691748
#> 683_Luzp1              0.95365528
#> 684_Gas2l1             1.00000000
#> 685_Pdcd11             0.98272300
#> 686_Chaf1b             0.89691748
#> 687_Pycr1              0.85184664
#> 688_Phf8               0.61804675
#> 689_Raver1             0.98272300
#> 690_Dbn1               0.86650073
#> 691_Dut                0.65273244
#> 692_Prpf4b             0.78401202
#> 693_Prpf4b             0.86650073
#> 694_Efs                0.96317292
#> 695_NA                 0.82069736
#> 696_Ppp2r5b            0.98272300
#> 697_Caskin2            0.46858598
#> 698_Arhgef40           0.65273244
#> 699_Zyx                0.89691748
#> 700_Mphosph10          0.53998241
#> 701_LOC113833392       0.89793287
#> 702_Cdc42ep1           0.90249140
#> 703_Snrpa1             0.74897796
#> 704_Ncbp1              0.95822477
#> 705_Gas2l1             1.00000000
#> 706_Gas2l1             0.93504308
#> 707_Bap1               0.89691748
#> 708_LOC100759640       0.96012819
#> 709_Cherp              0.89691748
#> 710_Nbn                0.53647630
#> 711_LOC100759640       1.00000000
#> 712_NA                 0.53647630
#> 713_Eif3b              1.00000000
#> 714_Miga2              0.86650073
#> 715_Prpf4b             0.61127133
#> 716_Dbn1               0.53998241
#> 717_Ppp2r5b            0.66402033
#> 718_Exosc9             0.75851223
#> 719_Eif3b              0.89691748
#> 720_Ripk2              0.81153692
#> 721_Dlg1               0.91740623
#> 722_N4bp1              0.97848956
#> 723_Nudc               0.61127133
#> 724_Znf367             1.00000000
#> 725_Ring1              1.00000000
#> 726_Snrpa1             0.65273244
#> 727_U2surp             0.80822413
#> 728_LOC100764225       0.59107576
#> 729_Cdc42ep1           0.59107576
#> 730_Znf385a            0.61127133
#> 731_Ints1              0.61127133
#> 732_LOC113833392       0.61127133
#> 733_Lrch4              0.83634429
#> 734_Ctdspl2            0.53647630
#> 735_Prpf4b             0.89691748
#> 736_Luzp1              0.84433926
#> 737_Eif3b              0.98272300
#> 738_Ptpn14             0.76552820
#> 739_Rrp1               0.89691748
#> 740_Lrrfip2            0.46858598
#> 741_Nsfl1c             0.89691748
#> 742_Ddx51              0.82801076
#> 743_Prpf38b            0.55406551
#> 744_Eef1b2             0.89691748
#> 745_Znf385a            0.95365528
#> 746_Map9               0.89691748
#> 747_Rflnb              0.89691748
#> 748_NA                 0.96425465
#> 749_C1H12orf45         0.96425465
#> 750_U2surp             0.70134699
#> 751_Caskin2            0.61127133
#> 752_Eri1               0.80519954
#> 753_Gsk3b              0.81456911
#> 754_LOC100766946       0.98272300
#> 755_Cnpy3              0.99123043
#> 756_Hnrnpc             1.00000000
#> 757_Ptpn14             1.00000000
#> 758_Slc7a11            0.96425465
#> 759_Hnrnpc             0.94766133
#> 760_Cdc37l1            0.76463145
#> 761_LOC100768405       0.82069736
#> 762_Rragc              0.82069736
#> 763_LOC113834282       0.71200897
#> 764_Fus                0.87684327
#> 765_Ubxn1              0.94365155
#> 766_Mmut               0.99821395
#> 767_Pdcd11             0.99040703
#> 768_LOC100757535       0.75851223
#> 769_Eif3b              0.53998241
#> 770_Rnf113a            0.96317292
#> 771_Sytl4              0.74825817
#> 772_Tlnrd1             0.90249140
#> 773_H671_1g1131        0.80007140
#> 774_Neurl1             0.86563801
#> 775_Zyx                1.00000000
#> 776_Ctdspl2            0.64380397
#> 777_Chaf1b             0.53647630
#> 778_Rragc              0.65421849
#> 779_Srfbp1             0.46858598
#> 780_Gys1               0.92300444
#> 781_Usp15              0.98272300
#> 782_Arhgef40           0.87684327
#> 783_Gigyf1             0.95365528
#> 784_Minar1             0.84433926
#> 785_Dus2               0.76552820
#> 786_Gatad2b            0.97937571
#> 787_Eif5               0.89691748
#> 788_Epb41l2            1.00000000
#> 789_Arl6ip4            0.87684327
#> 790_Plin4              0.82069736
#> 791_Elf2               0.85184664
#> 792_Plin4              0.54346611
#> 793_Snip1              0.81456911
#> 794_Txlng              0.62865083
#> 795_LOC100769437       0.53647630
#> 796_Caskin2            0.65273244
#> 797_NA                 0.65273244
#> 798_Synm               0.81550950
#> 799_Synm               1.00000000
#> 800_Ube2c              1.00000000
#> 801_Sgtb               0.92872331
#> 802_Prpf4b             0.61127133
#> 803_Epb41l2            0.71378983
#> 804_Mllt1              0.43271982
#> 805_LOC100759640       1.00000000
#> 806_Epb41l2            0.95731924
#> 807_Znf280b            0.87684327
#> 808_Kiaa1143           1.00000000
#> 809_Gas2l1             0.71891081
#> 810_Srp72              0.87684327
#> 811_Tomm22             0.87684327
#> 812_Psip1              0.92267626
#> 813_Arhgef37           0.53647630
#> 814_Bckdk              0.97937571
#> 815_Strip1             0.71378983
#> 816_Usp15              0.78218056
#> 817_Ssr3               0.82973183
#> 818_Strip1             0.89691748
#> 819_Eif3b              0.61804675
#> 820_U2surp             0.81456911
#> 821_Bend3              0.89691748
#> 822_Rps10              1.00000000
#> 823_Rpl23a             0.78605465
#> 824_Nbn                0.70208246
#> 825_Rpap3              0.96317292
#> 826_LOC100759640       0.53998241
#> 827_Ric8a              0.58004623
#> 828_Hsph1              0.84433926
#> 829_LOC100759640       0.96317292
#> 830_LOC100757535       0.91418289
#> 831_Gigyf1             0.71378983
#> 832_Dbn1               1.00000000
#> 833_Snrk               1.00000000
#> 834_Prpf38b            0.99040703
#> 835_LOC100766868       1.00000000
#> 836_LOC100766868       0.53998241
#> 837_Wbp11              0.77944157
#> 838_Rusc2              0.61127133
#> 839_Eif3b              0.98969498
#> 840_Ptpn14             0.96425465
#> 841_Rlim               1.00000000
#> 842_Ints1              0.71378983
#> 843_Chaf1b             0.46858598
#> 844_Dlg1               1.00000000
#> 845_Lamtor1            0.97937571
#> 846_Tab1               0.99040703
#> 847_Dbn1               0.53998241
#> 848_Psip1              0.53998241
#> 849_Dbn1               0.99790083
#> 850_Pabpc1             0.99040703
#> 851_Hnrnpc             0.89691748
#> 852_Emd                0.76463145
#> 853_LOC100764225       0.61127133
#> 854_Nup50              0.76463145
#> 855_Ctcf               0.71891081
#> 856_Raly               0.96317292
#> 857_Bard1              0.87684327
#> 858_Ptpn14             0.87125683
#> 859_LOC100757535       0.43271982
#> 860_Psmd2              1.00000000
#> 861_Junb               0.53998241
#> 862_C1qbp              0.97937571
#> 863_Lrch4              0.98891432
#> 864_CUNH14orf93        0.82069736
#> 865_U2surp             0.66402033
#> 866_Raly               0.98272300
#> 867_LOC100774417       0.80605441
#> 868_Srp72              1.00000000
#> 869_LOC100764225       1.00000000
#> 870_Morf4l2            0.68492351
#> 871_CUNH9orf40         0.71378983
#> 872_Gas2l1             0.57624609
#> 873_Atp5pf             0.87684327
#> 874_Lrrfip2            0.71378983
#> 875_Prpf4b             0.62405509
#> 876_Top2b              0.53998241
#> 877_Mepce              0.84433926
#> 878_Ptpn14             0.91734510
#> 879_Dnajc25            0.46858598
#> 880_Cbx8               0.95731924
#> 881_Synm               0.82203853
#> 882_Def6               1.00000000
#> 883_Gys1               0.96317292
#> 884_Luzp1              0.76463145
#> 885_Synm               0.87184148
#> 886_Snip1              0.65273244
#> 887_Top2b              0.80741720
#> 888_NA                 1.00000000
#> 889_Trim35             1.00000000
#> 890_Znf385a            0.96317292
#> 891_Chaf1b             0.71891081
#> 892_Abcf1              1.00000000
#> 893_Pdcd11             0.95365528
#> 894_Dlg1               0.79131251
#> 895_Dbn1               0.69476960
#> 896_LOC100752363       0.92669903
#> 897_Ppp4r3a            0.59107576
#> 898_Gas2l1             0.76552820
#> 899_Mtmr10             1.00000000
#> 900_Cyld               0.61127133
#> 901_NA                 0.76552820
#> 902_Rnf113a            0.53647630
#> 903_Nelfa              0.76463145
#> 904_Zkscan1            0.81550950
#> 905_Chaf1b             0.65273244
#> 906_Eif3b              1.00000000
#> 907_Top2b              0.65273244
#> 908_Chaf1b             0.98272300
#> 909_Epb41l2            0.99616460
#> 910_C3H11orf58         0.76463145
#> 911_Top2b              0.99040703
#> 912_Wee1               0.99040703
#> 913_Raly               0.98891432
#> 914_H671_1g2680        0.91418289
#> 915_Eef1b2             0.65273244
#> 916_Gas2l1             0.46858598
#> 917_Epb41l2            0.61127133
#> 918_Rpl23a             0.97937571
#> 919_Chmp2b             0.98272300
#> 920_Lrrfip2            0.46858598
#> 921_Aldoa              0.61127133
#> 922_Cby1               0.76552820
#> 923_LOC100759640       0.96317292
#> 924_Rbm28              0.99616460
#> 925_Skiv2l             0.76552820
#> 926_Ints1              0.95731924
#> 927_Ehd1               1.00000000
#> 928_Nr2f6              0.84433926
#> 929_Top2b              0.53998241
#> 930_Lrrfip2            0.96317292
#> 931_Pip4p2             0.92300444
#> 932_Srp72              0.76463145
#> 933_Mtmr9              0.65273244
#> 934_Gigyf1             0.98380907
#> 935_Rbm7               0.87684327
#> 936_LOC100773565       0.74897796
#> 937_Trim35             0.86650073
#> 938_Cbx8               0.81456911
#> 939_Rplp0              1.00000000
#> 940_Aldoa              0.74897796
#> 941_NA                 0.94365155
#> 942_Zyx                0.76463145
#> 943_Psip1              0.95897789
#> 944_Slc7a11            0.96317292
#> 945_Miga2              0.61127133
#> 946_Arhgef6            0.53998241
#> 947_Dlgap4             0.95365528
#> 948_Ampd2              1.00000000
#> 949_Luzp1              0.91740623
#> 950_Camlg              0.81456911
#> 951_Pfkfb3             0.73477355
#> 952_NA                 1.00000000
#> 953_Raly               1.00000000
#> 954_Kiaa1143           0.43271982
#> 955_Bcar1              0.76552820
#> 956_Gatad2b            0.74897796
#> 957_Eif4ebp2           0.98272300
#> 958_Fam76b             0.57736706
#> 959_Camlg              0.82069736
#> 960_LOC100754077       0.31331155
#> 961_NA                 0.43271982
#> 962_Epb41l2            0.99500130
#> 963_Ankrd34a           1.00000000
#> 964_Zc3h15             0.53647630
#> 965_Def6               0.76552820
#> 966_Srsf6              0.69265742
#> 967_H671_4g11480       0.76463145
#> 968_Top2b              0.88019189
#> 969_LOC100769471       0.87684327
#> 970_Raver1             0.81550950
#> 971_Etv3               1.00000000
#> 972_Psd                0.69476960
#> 973_Usp15              0.99500130
#> 974_Nol7               0.89691748
#> 975_Stk38              0.89691748
#> 976_Smim13             0.98272300
#> 977_Etv3               0.80741720
#> 978_Synm               0.96317292
#> 979_Pwp1               0.76463145
#> 980_Fus                0.87684327
#> 981_Junb               0.89691748
#> 982_Phf8               0.66699560
#> 983_Nelfa              0.98272300
#> 984_Prpf4b             0.87435458
#> 985_Abraxas1           0.87684327
#> 986_Prpf4b             1.00000000
#> 987_Raver1             1.00000000
#> 988_Caap1              0.94365155
#> 989_Rpap3              0.89691748
#> 990_Hsph1              0.69265742
#> 991_LOC100750437       1.00000000
#> 992_Mepce              0.75851223
#> 993_Efs                0.95731924
#> 994_Epb41l2            0.92593443
#> 995_Abcf1              0.80822413
#> 996_NA                 0.86319749
#> 997_Eif4ebp2           1.00000000
#> 998_Pfkfb3             0.98272300
#> 999_Hnrnpc             0.89691748
#> 1000_Psmd2             0.99500130
#> 
#> $Exponential$pvc_pattern_summary
#>   -60 15 60 90 120 240
#> p   0  0  0  1   1   0
#> v   0  1  0  2   0   0
#> b   0  0  1  0   0   0
#> t   0  0  0  0   0   0
#> 
#> 
#> $Stationary
#> $Stationary$pvc_adj_pvals
#>                  15_vs_neighbors 60_vs_neighbors 90_vs_neighbors
#> 1_Top2b               0.47623282     0.010069485       0.2528074
#> 2_NA                  1.00000000     1.000000000       1.0000000
#> 3_Snip1               1.00000000     1.000000000       0.9914096
#> 4_Tomm34              0.97399727     0.730884695       0.8594514
#> 5_Pus3                0.93408439     0.612667440       1.0000000
#> 6_Ints1               0.90515009     0.860650480       0.9790592
#> 7_Mlh1                0.93408439     0.949911982       0.7816496
#> 8_LOC100750437        1.00000000     1.000000000       0.2528074
#> 9_Pabpc1              0.54502481     1.000000000       0.7774182
#> 10_Top2b              0.99508555     0.848200886       0.7605719
#> 11_Gorasp1            0.53591151     0.212803322       0.8452398
#> 12_Ints1              1.00000000     1.000000000       1.0000000
#> 13_Syvn1              0.60118421     0.457565692       0.8538114
#> 14_Znf280b            0.63289895     0.536829699       0.7656499
#> 15_Mrnip              0.79319446     0.616363428       0.9998830
#> 16_Rragc              0.90515009     0.839109126       0.7012143
#> 17_Gorasp1            0.71724978     0.119293142       0.7305634
#> 18_Tomm34             0.83839993     0.481902036       0.6563396
#> 19_LOC100757430       0.62388094     0.127441307       0.6923074
#> 20_Ubxn1              0.88292911     0.368330199       0.8049503
#> 21_H671_1g1131        0.54502481     0.400510484       0.2456534
#> 22_Luzp1              0.58745270     0.502705746       0.9298273
#> 23_Efs                0.54502481     0.048783400       0.3971096
#> 24_Mta2               0.94427730     0.523907271       0.7685783
#> 25_Nedd1              0.71607132     0.127428892       0.3648974
#> 26_Gigyf1             0.54325567     0.130931101       0.7416436
#> 27_Myh9               0.90515009     0.093865695       0.2456534
#> 28_Caskin2            0.64073199     0.028561245       0.2528074
#> 29_Papolg             0.54325567     0.310570702       0.3648974
#> 30_Tfg                0.65406828     0.834107317       0.6170691
#> 31_Rpl34              0.52376249     0.113354095       0.5515318
#> 32_Mideas             0.90515009     0.769933446       0.6464593
#> 33_Gys1               0.29488064     0.004123866       0.1865649
#> 34_Arhgef6            0.85192356     0.170774257       0.2792817
#> 35_Ctdspl2            1.00000000     1.000000000       1.0000000
#> 36_Ptpn14             1.00000000     1.000000000       0.6451309
#> 37_Raly               1.00000000     0.996871421       0.9965477
#> 38_Znhit3             0.63851869     0.082558457       0.3226632
#> 39_LOC113833392       0.38373964     0.024970798       1.0000000
#> 40_Luc7l3             0.61241133     0.118235739       0.4873530
#> 41_Rplp0              1.00000000     1.000000000       1.0000000
#> 42_Gys1               0.54502481     0.507654023       0.9790592
#> 43_Rpl22l1            0.49554493     0.299407772       0.7888888
#> 44_Eif3b              0.90515009     0.818415544       0.7073961
#> 45_Med26              0.69624464     0.040469441       0.5214596
#> 46_Mepce              0.97194962     0.874141663       0.7728969
#> 47_Pdcd11             0.63993993     0.307212373       0.8613130
#> 48_Twf1               0.53917994     0.217102161       0.2528074
#> 49_LOC100759640       0.73891031     1.000000000       1.0000000
#> 50_Wrnip1             0.74053880     0.329225607       0.7605719
#> 51_Poldip3            0.52376249     0.047911897       0.4629072
#> 52_Ampd2              0.88137715     0.554525769       0.7038252
#> 53_Mea1               0.79677029     0.949911982       0.9905692
#> 54_Dbn1               0.66796513     0.469833614       0.8939102
#> 55_Snip1              0.97194962     0.981399784       0.9368268
#> 56_Srsf6              0.52063115     0.011705607       0.3975070
#> 57_LOC113834282       1.00000000     1.000000000       1.0000000
#> 58_Map9               0.63385599     0.085099598       0.5004048
#> 59_Cdc42ep1           0.66005075     0.761503234       0.9790592
#> 60_Poldip3            1.00000000     1.000000000       1.0000000
#> 61_LOC100764225       1.00000000     1.000000000       1.0000000
#> 62_Epb41l2            1.00000000     1.000000000       1.0000000
#> 63_H671_4g11480       0.76257173     0.958954018       0.9757581
#> 64_Nbn                1.00000000     1.000000000       1.0000000
#> 65_U2surp             1.00000000     1.000000000       1.0000000
#> 66_Gigyf1             0.54325567     0.073928051       0.3779698
#> 67_NA                 0.92364395     0.554712047       0.7121590
#> 68_Luc7l3             1.00000000     1.000000000       1.0000000
#> 69_LOC100752363       0.72613118     0.101906724       0.4771061
#> 70_Ampd2              1.00000000     1.000000000       1.0000000
#> 71_LOC100759640       1.00000000     1.000000000       1.0000000
#> 72_Stam               0.92291089     0.355297643       0.3614870
#> 73_Nsfl1c             0.54502481     0.678847456       0.4980875
#> 74_Pfkfb3             0.93408439     0.208936896       0.2788533
#> 75_Rad23a             0.92291089     0.041040655       0.2528074
#> 76_Elf2               0.62288573     0.009592916       0.2456534
#> 77_Crem               0.90515009     0.872113131       0.9316856
#> 78_Rragc              0.57155299     0.277582354       0.4412340
#> 79_Lrrfip2            0.93408439     0.934269330       0.9998830
#> 80_Zyx                0.90515009     0.383464190       0.6505823
#> 81_Lrrfip2            0.46779026     0.085215420       0.4586578
#> 82_Gatad2b            0.89679055     0.115365078       0.5098348
#> 83_Bcar1              0.92291089     0.709226993       0.9546590
#> 84_Ehd1               0.54325567     0.025047244       0.3226632
#> 85_LOC113834282       0.53917994     0.048783400       0.4046589
#> 86_Tmem230            0.97878789     0.218109162       0.7605719
#> 87_Ncbp1              0.62945080     0.044152701       0.3966787
#> 88_Mllt1              1.00000000     1.000000000       1.0000000
#> 89_Stk17b             0.73942586     0.099296146       0.3756990
#> 90_Dlgap4             0.97878789     0.218109162       0.7605719
#> 91_Papolg             0.84802386     0.093254560       0.2456534
#> 92_Cyld               0.80885508     0.218109162       0.7605719
#> 93_Gigyf1             0.64073199     0.854188318       0.9164056
#> 94_Lrrfip2            1.00000000     1.000000000       1.0000000
#> 95_Lrrfip2            0.97194962     0.100570871       0.3779698
#> 96_Rlim               0.70303031     0.990436120       0.7969119
#> 97_Eif3b              0.71607132     0.980519065       0.7857770
#> 98_Mphosph10          0.64552746     0.375556171       1.0000000
#> 99_Gatad2b            0.39686202     0.161989284       0.7686332
#> 100_Srsf6             0.60118421     0.062796911       0.5474512
#> 101_Zyx               0.88292911     0.386375070       0.7774182
#> 102_Mphosph10         1.00000000     1.000000000       1.0000000
#> 103_Psip1             0.81804469     0.834107317       0.9542486
#> 104_Fbl               0.99508555     0.835850865       0.9316856
#> 105_H671_1g2680       0.82649872     0.818669157       0.4959790
#> 106_Sgtb              0.62288573     0.881643420       0.9905692
#> 107_Gnl3              0.93408439     0.218109162       0.6170691
#> 108_Eif3b             0.63385599     0.492326464       0.9090431
#> 109_Serpinb1          1.00000000     1.000000000       1.0000000
#> 110_N4bp1             1.00000000     1.000000000       0.7991360
#> 111_Snip1             1.00000000     1.000000000       1.0000000
#> 112_Psip1             1.00000000     1.000000000       1.0000000
#> 113_Mlh1              0.90515009     0.473245549       0.2792817
#> 114_Bsg               0.97194962     0.646473799       0.9606716
#> 115_Tnpo1             0.72880676     0.809159384       0.8049503
#> 116_H671_1g2680       1.00000000     1.000000000       1.0000000
#> 117_Cbx8              1.00000000     1.000000000       1.0000000
#> 118_Mideas            1.00000000     1.000000000       1.0000000
#> 119_Mideas            0.97194962     0.118412423       0.3678174
#> 120_Dcun1d3           0.97194962     0.118412423       0.3678174
#> 121_Dlg1              1.00000000     1.000000000       1.0000000
#> 122_Rad23a            0.54325567     0.024200859       0.5683353
#> 123_Srsf6             0.72613118     0.231705463       0.9790592
#> 124_Stx7              1.00000000     1.000000000       1.0000000
#> 125_Pdcd11            0.88292911     0.726453304       0.6138471
#> 126_Kiaa1958          0.90515009     0.799298832       0.9790592
#> 127_Pwp1              0.89639358     0.822406373       0.8487388
#> 128_Txlng             0.71607132     0.790034699       0.9998830
#> 129_Junb              1.00000000     1.000000000       1.0000000
#> 130_LOC100759640      1.00000000     1.000000000       1.0000000
#> 131_Dbn1              0.45454337     0.048653042       0.7305634
#> 132_Top2b             0.62288573     0.461331485       0.8184153
#> 133_Rusc2             0.93408439     0.399784993       0.6138471
#> 134_NA                0.94606430     1.000000000       1.0000000
#> 135_LOC113837251      0.88563466     0.588310352       0.9857687
#> 136_Fam76b            1.00000000     1.000000000       1.0000000
#> 137_Ptpn14            0.75827918     0.098525135       0.3779698
#> 138_Chmp4b            0.54325567     0.028456138       0.8258043
#> 139_Prpf4b            1.00000000     1.000000000       1.0000000
#> 140_Eif3b             1.00000000     1.000000000       1.0000000
#> 141_Nsfl1c            1.00000000     1.000000000       1.0000000
#> 142_Pdlim7            1.00000000     1.000000000       1.0000000
#> 143_Rnf113a           1.00000000     1.000000000       1.0000000
#> 144_Epb41l2           0.54094413     0.222577072       0.9905692
#> 145_Hnrnpc            1.00000000     1.000000000       1.0000000
#> 146_LOC113834282      0.54325567     0.044152701       0.6170691
#> 147_Plekho2           0.88292911     0.322882788       0.9724011
#> 148_Med26             0.93408439     0.936900048       0.6675377
#> 149_Arhgef40          0.90515009     0.796500636       0.7321951
#> 150_NA                1.00000000     1.000000000       0.8566426
#> 151_Phf8              0.46779026     0.010131878       0.2211880
#> 152_Minar1            0.79319446     1.000000000       1.0000000
#> 153_H671_21690        0.88292911     0.553896530       0.5683353
#> 154_Arhgef40          0.51979551     0.018650257       0.4959790
#> 155_Chaf1b            0.99983034     0.589597408       0.5683353
#> 156_Prpf4b            0.62479340     0.429382386       0.8786434
#> 157_Znf367            0.64828195     0.838455031       0.9740321
#> 158_Luzp1             0.70908163     0.804669079       0.6138471
#> 159_LOC113833882      0.57155299     0.234312222       0.7605719
#> 160_Hnrnpc            1.00000000     1.000000000       1.0000000
#> 161_Mepce             0.53591151     0.095752128       0.3678174
#> 162_Ubxn1             0.93408439     0.947257152       0.6260825
#> 163_Mllt1             0.05820201     0.024970798       0.6749909
#> 164_Chaf1b            0.54094413     0.268296359       0.9905692
#> 165_Raly              0.67549601     0.183136137       0.6669812
#> 166_Gas2l1            0.66502052     0.368796409       0.7605719
#> 167_Dlg1              0.88292911     0.804749173       0.8303413
#> 168_Hoxc10            1.00000000     1.000000000       0.6995470
#> 169_Gigyf1            0.67190462     1.000000000       1.0000000
#> 170_Luzp1             0.71607132     0.495954955       0.7660779
#> 171_Srp72             0.64073199     0.152487189       0.7658717
#> 172_LOC100771461      0.54325567     0.342023429       1.0000000
#> 173_Chaf1b            0.71794780     0.568194019       0.9316856
#> 174_C3H11orf58        0.78137419     0.799298832       1.0000000
#> 175_Pdcd11            1.00000000     1.000000000       1.0000000
#> 176_Psip1             0.62288573     0.822406373       0.2528074
#> 177_Prpf4b            0.97194962     0.231298837       0.3614870
#> 178_Rnf113a           0.57155299     0.044152701       0.4980875
#> 179_Irf3              0.94427730     0.752570747       0.5683353
#> 180_Smim13            0.88292911     0.834107317       0.5369183
#> 181_Gnl3              0.94012802     0.589597408       0.4441837
#> 182_Psma5             0.93408439     0.728733616       0.6223072
#> 183_Ptpn14            0.57155299     0.044152701       0.1973920
#> 184_Prpf4b            1.00000000     0.132571176       0.3917706
#> 185_Top2b             1.00000000     1.000000000       1.0000000
#> 186_Prpf38b           0.65406828     0.072705811       0.2456534
#> 187_Epb41l2           0.93408439     0.864692664       0.6995470
#> 188_Eif3b             1.00000000     1.000000000       1.0000000
#> 189_Hnrnpc            0.85192356     0.962032310       0.7423490
#> 190_LOC100758278      0.89679055     0.260257950       0.3913066
#> 191_Prpf4b            0.92291089     0.818669157       0.5734605
#> 192_Caskin2           0.58496599     0.477961599       0.4010473
#> 193_LOC100752363      0.54325567     0.700740104       0.8184153
#> 194_Septin6           0.93408439     0.818415544       0.6532590
#> 195_Max               0.62388094     0.329225607       0.9688836
#> 196_Mid1ip1           0.61241133     0.093652409       0.5515318
#> 197_NA                0.94792345     0.399784993       0.8391381
#> 198_Hsph1             1.00000000     1.000000000       1.0000000
#> 199_Nol7              0.32526276     0.015678224       0.2528074
#> 200_Raly              0.52376249     0.405200601       0.9149144
#> 201_Smim13            0.70908163     0.944038423       0.5382447
#> 202_LOC100757535      0.65406828     0.093254560       0.2535009
#> 203_Net1              0.89679055     0.072788586       0.2456534
#> 204_LOC100754077      1.00000000     1.000000000       1.0000000
#> 205_Snip1             0.46779026     0.249506635       0.7309207
#> 206_Hnrnpc            0.62288573     0.312784861       0.9905692
#> 207_Ldlrap1           0.85192356     0.868830731       0.8613130
#> 208_Luzp1             0.93408439     0.118412423       0.2792817
#> 209_Rpl26             0.93408439     0.847488627       0.3779698
#> 210_Epb41l2           0.89639358     0.262881922       0.6532590
#> 211_Znf367            0.81448384     0.945572261       0.9585265
#> 212_Dlgap4            1.00000000     1.000000000       1.0000000
#> 213_Plekho2           0.76257173     1.000000000       1.0000000
#> 214_Zpr1              0.92291089     0.804069909       0.8109574
#> 215_Dlgap4            1.00000000     1.000000000       1.0000000
#> 216_Def6              0.94427730     0.754718745       0.9740321
#> 217_Eif4ebp2          0.62479340     0.412411042       0.7605719
#> 218_Eef1b2            0.57311992     0.690767762       0.8594514
#> 219_Rad23a            0.76257173     0.616238304       0.6995470
#> 220_Morf4l2           0.72613118     0.818415544       0.9790592
#> 221_Arhgef40          0.72613118     0.447603231       0.8711279
#> 222_NA                0.78223484     0.585388668       0.4144947
#> 223_LOC100773565      0.77895142     0.799298832       0.7305634
#> 224_Dus2              0.78315628     0.088940513       0.6563396
#> 225_Pip4p2            0.45454337     0.141221624       0.7564088
#> 226_Top2b             0.61035300     0.934269330       0.9894765
#> 227_Znf280b           0.72251712     0.306129267       0.9790592
#> 228_Pdcd11            0.89867219     0.368796409       0.7495548
#> 229_Bckdk             0.79319446     0.937388320       1.0000000
#> 230_Arhgef40          1.00000000     1.000000000       1.0000000
#> 231_Mepce             1.00000000     1.000000000       1.0000000
#> 232_Ccnd3             0.98509912     0.548421073       0.6995470
#> 233_Phf8              1.00000000     1.000000000       1.0000000
#> 234_H671_1g2680       0.99983034     0.813074192       0.9790592
#> 235_Ell               0.93408439     0.991723987       0.8594514
#> 236_U2surp            0.52376249     0.422107012       0.8487388
#> 237_Rps10             0.76257173     0.809465395       0.9905692
#> 238_Ctdspl2           0.54502481     0.144987365       0.9316856
#> 239_Top2b             1.00000000     1.000000000       1.0000000
#> 240_Msantd3           0.62945080     0.938325214       0.3614870
#> 241_Fam76b            0.97399727     0.754718745       0.9790592
#> 242_Ppp4r3a           1.00000000     1.000000000       1.0000000
#> 243_Gpatch4           1.00000000     1.000000000       0.7660779
#> 244_Nudc              0.97194962     0.804749173       0.7245342
#> 245_Nol7              0.82667600     0.460414349       0.5528622
#> 246_Plekho2           0.89639358     0.049525675       0.3614870
#> 247_Prpf4b            0.71794780     0.116295046       0.7073961
#> 248_Mta2              0.64073199     0.285687994       0.8487388
#> 249_U2surp            0.62388094     0.082341425       0.9740321
#> 250_Ubxn1             1.00000000     1.000000000       1.0000000
#> 251_Rlim              0.97194962     0.337140363       0.2817729
#> 252_Atat1             0.54094413     0.237143782       0.7321951
#> 253_Ubxn1             0.98602325     0.798374776       0.8757057
#> 254_H671_1g2680       1.00000000     0.619146987       0.5098348
#> 255_eIF2aK2           0.93975237     0.100881143       0.2792817
#> 256_Skiv2l            0.63851869     0.110015956       1.0000000
#> 257_Rpl28             0.54325567     0.645812581       1.0000000
#> 258_LOC100759640      1.00000000     1.000000000       1.0000000
#> 259_Gatad2b           0.97878789     0.834331354       0.8613130
#> 260_NA                0.62945080     0.383464190       0.9998830
#> 261_Gprasp1           0.93975237     0.589927693       0.8883499
#> 262_Luzp1             1.00000000     1.000000000       1.0000000
#> 263_Slc1a5            1.00000000     1.000000000       1.0000000
#> 264_LOC113834282      0.86639744     0.834107317       0.6110255
#> 265_Srsf6             0.65406828     0.914920513       0.7605719
#> 266_Cdc42ep1          0.72880676     0.378206060       0.6923074
#> 267_Net1              1.00000000     1.000000000       0.5048485
#> 268_Caskin2           0.24835879     0.082341425       0.7964405
#> 269_LOC100759640      0.88292911     0.046953528       0.2762895
#> 270_Mideas            0.70908163     0.859444878       0.4959790
#> 271_Luzp1             0.72613118     0.763870844       0.7746328
#> 272_Emd               1.00000000     1.000000000       1.0000000
#> 273_Plpp6             0.53917994     0.044182198       0.7685783
#> 274_LOC100759640      1.00000000     1.000000000       1.0000000
#> 275_Rps7              0.71607132     0.158835844       0.5528622
#> 276_Fkbp1a            0.54094413     0.554525769       0.6138471
#> 277_Gatad2b           0.97878789     0.307629261       0.3913066
#> 278_Znf385a           0.45454337     0.056503461       0.6494172
#> 279_Arhgef6           0.93408439     0.427929351       0.9905692
#> 280_Slirp             0.94276860     0.935348601       0.8580075
#> 281_Skiv2l            0.75621416     0.985476771       0.8315399
#> 282_H671_21690        0.89639358     0.368796409       0.5098348
#> 283_Kat8              0.54325567     0.078777292       0.6923074
#> 284_Nkap              0.82649872     0.352472407       0.8013220
#> 285_Gsk3b             0.88292911     0.114460768       0.2841615
#> 286_Ints1             0.54325567     0.137187744       0.6563396
#> 287_Gas2l1            0.97935078     1.000000000       1.0000000
#> 288_LOC100759640      1.00000000     1.000000000       1.0000000
#> 289_Top2b             0.72880676     0.024200859       0.3022223
#> 290_Kif20b            0.56867904     0.256751492       0.5853013
#> 291_Phf8              1.00000000     1.000000000       1.0000000
#> 292_Snip1             1.00000000     0.164093919       0.2528074
#> 293_Gsk3b             1.00000000     1.000000000       1.0000000
#> 294_Caskin2           1.00000000     1.000000000       1.0000000
#> 295_C3H11orf58        0.53917994     0.251777985       0.4258505
#> 296_Lrch4             0.71724978     0.249506635       0.7788113
#> 297_LOC113834282      0.54502481     0.218109162       0.6140788
#> 298_LOC100750407      1.00000000     1.000000000       1.0000000
#> 299_LOC113833392      1.00000000     1.000000000       1.0000000
#> 300_LOC113833882      0.89410612     0.256944312       0.4771061
#> 301_Ldlrap1           0.54325567     0.145603610       0.3175328
#> 302_Wee1              0.88292911     0.761386308       0.9790592
#> 303_Caap1             1.00000000     1.000000000       1.0000000
#> 304_Eif4ebp2          0.99508555     0.937520812       0.9790592
#> 305_Ripk2             0.97301371     0.447195064       0.3975070
#> 306_Srp72             0.76257173     0.388957672       0.5683353
#> 307_Taok2             0.85518808     0.941668344       0.6604538
#> 308_Nr2f6             0.45454337     0.010428941       0.2528074
#> 309_Arhgef40          0.27742871     0.079178318       0.5098348
#> 310_Gys1              1.00000000     1.000000000       1.0000000
#> 311_Dlg1              1.00000000     1.000000000       1.0000000
#> 312_Vapb              0.79319446     0.839336901       0.5683353
#> 313_LOC100757535      0.53917994     0.010428941       0.2456534
#> 314_Mkrn2             0.90515009     0.913964935       0.7434850
#> 315_Eif3b             0.87609268     0.225209194       0.5515318
#> 316_Isyna1            0.94543672     1.000000000       1.0000000
#> 317_Prpf4b            0.71724978     0.044182198       0.3175328
#> 318_LOC113833882      0.24835879     0.015678224       0.7073961
#> 319_Lrch4             0.56592593     0.181894212       0.8585766
#> 320_Dbn1              0.69475214     0.140667781       0.5528622
#> 321_Abcf1             0.52376249     0.044152701       0.3175328
#> 322_Ints1             0.32526276     0.011705607       0.2528074
#> 323_C3H11orf58        1.00000000     1.000000000       1.0000000
#> 324_Psma5             0.71524930     0.302100917       0.5677689
#> 325_Fundc1            0.98602325     0.778244792       0.7605719
#> 326_Papolg            0.83839993     1.000000000       1.0000000
#> 327_Mideas            0.93408439     0.368194560       0.8038231
#> 328_Ubxn1             0.71633877     0.546856546       0.8210004
#> 329_Synm              0.85192356     0.774335787       0.9790592
#> 330_Arhgef6           1.00000000     1.000000000       1.0000000
#> 331_Ptpn14            1.00000000     0.343423518       0.7073961
#> 332_Pgrmc1            1.00000000     0.343423518       0.7073961
#> 333_Myh9              0.88292911     0.915527348       0.6563396
#> 334_Etv3              1.00000000     1.000000000       1.0000000
#> 335_Ip6k1             1.00000000     1.000000000       1.0000000
#> 336_Luzp1             1.00000000     1.000000000       1.0000000
#> 337_Ptpn14            1.00000000     1.000000000       1.0000000
#> 338_Caskin2           1.00000000     1.000000000       1.0000000
#> 339_Chaf1b            1.00000000     1.000000000       0.5956977
#> 340_Ubxn1             1.00000000     1.000000000       1.0000000
#> 341_Ube2c             1.00000000     1.000000000       1.0000000
#> 342_Gins2             0.70908163     0.902538619       0.9032349
#> 343_Nlgn2             0.62238996     0.676683041       0.9724011
#> 344_Nf2               0.54325567     0.105583422       0.5683353
#> 345_Pip4p2            0.54502481     0.285227044       0.7073961
#> 346_Emd               1.00000000     1.000000000       1.0000000
#> 347_Top2b             0.81864075     0.393617907       0.7660779
#> 348_Trim35            0.93408439     0.285818903       0.2528074
#> 349_NA                0.66796513     0.635527827       0.3226632
#> 350_NA                0.75827918     0.137491563       1.0000000
#> 351_Mideas            1.00000000     1.000000000       1.0000000
#> 352_Gas2l1            1.00000000     0.124564801       0.2841615
#> 353_Ampd2             0.88292911     0.453321923       0.9332729
#> 354_Calu              1.00000000     1.000000000       1.0000000
#> 355_Fam76b            1.00000000     1.000000000       1.0000000
#> 356_Dlg1              0.65300837     0.079545642       0.4629072
#> 357_Srsf6             0.94792345     0.120554949       0.4555868
#> 358_Chaf1b            0.94792345     0.120554949       0.4555868
#> 359_Dbn1              0.66796513     0.337140363       0.3913066
#> 360_Tcf25             0.93975237     0.834107317       0.8452398
#> 361_Psip1             0.12918526     0.009592916       0.4629072
#> 362_Cnpy3             0.60253962     0.224001282       0.7882100
#> 363_LOC100759640      0.61241133     0.453321923       0.3975070
#> 364_Zyx               1.00000000     1.000000000       0.6771960
#> 365_Lrch4             0.62479340     0.912339550       0.6995470
#> 366_Bola1             1.00000000     1.000000000       1.0000000
#> 367_Znf385a           1.00000000     0.640415997       0.6260825
#> 368_Kif20b            0.98096368     0.457525130       0.6260825
#> 369_Ell               0.54502481     0.462802123       0.6433664
#> 370_Ell               0.77895142     0.316343659       0.7611982
#> 371_Srsf6             0.57186421     0.915557817       0.6715757
#> 372_Pwp1              0.88292911     0.033535898       0.2456534
#> 373_Def6              0.97878789     0.090761275       0.3226632
#> 374_Cbx8              0.66930492     0.752570747       0.8480630
#> 375_Ddx51             0.79319446     0.140667781       0.3447812
#> 376_Psip1             0.93408439     0.295813395       0.7415917
#> 377_Arhgef40          0.63454585     0.713568546       0.7858065
#> 378_Raly              0.89639358     0.818415544       0.9790592
#> 379_NA                0.54325567     0.226989403       0.6767307
#> 380_Lrrfip2           0.93408439     0.833342900       0.5034136
#> 381_Gnl3              0.58496599     0.133601735       0.5474512
#> 382_Caskin2           0.71633877     0.130931101       0.5515318
#> 383_Rragc             0.98100329     0.852754600       0.7305634
#> 384_Caskin2           0.97194962     0.950493002       0.9905692
#> 385_Bcar1             0.54094413     1.000000000       1.0000000
#> 386_Homer3            0.65640867     0.009592916       0.2456534
#> 387_Luzp1             0.78913563     0.048653042       0.6715757
#> 388_N4bp1             0.92291089     0.460414349       1.0000000
#> 389_Ppp4r3a           0.71607132     0.265196832       0.6171612
#> 390_H671_1g2680       0.94427730     0.848200886       0.5853013
#> 391_Gnl3              0.72459850     0.089170899       0.6464664
#> 392_Top2b             0.52376249     0.115348038       0.8613130
#> 393_Oser1             1.00000000     1.000000000       1.0000000
#> 394_Snrk              0.97301371     0.949911982       0.9657031
#> 395_Kat8              0.97763126     0.085358109       0.2456534
#> 396_Raver1            0.65300837     0.438379999       0.8613130
#> 397_Pdcd11            0.52063115     0.020193592       0.2841615
#> 398_Rps20             0.65406828     0.800283545       0.6494172
#> 399_Bsg               0.60253962     0.484162075       0.9332729
#> 400_Raly              0.76257173     0.340197749       0.7660779
#> 401_Pdcd2             0.97194962     0.876503074       0.5683353
#> 402_Caskin2           0.73942586     0.494439642       0.1973920
#> 403_LOC100773571      0.93408439     0.174869256       0.6995470
#> 404_Papolg            1.00000000     1.000000000       1.0000000
#> 405_LOC100757535      0.97878789     0.098525135       0.2528074
#> 406_Caap1             0.89688873     0.788507330       0.6178613
#> 407_Psip1             0.65406828     0.796500636       0.9790592
#> 408_Dbn1              0.71607132     0.267648475       0.8993151
#> 409_Mta2              0.46779026     0.086647449       0.8184153
#> 410_Abcf1             0.46779026     1.000000000       1.0000000
#> 411_LOC100754108      0.84802386     0.378690117       0.5515318
#> 412_Slirp             0.16916527     0.009592916       0.7714479
#> 413_Nelfa             0.85691664     0.120554949       0.3614870
#> 414_Aggf1             0.86772434     0.834107317       0.5567138
#> 415_Bap1              0.69625731     0.015227573       0.2456534
#> 416_Luc7l3            0.54325567     0.062795263       0.4017848
#> 417_Rrp1              0.45454337     0.011705607       0.3631106
#> 418_Wrnip1            1.00000000     1.000000000       0.7073961
#> 419_NA                0.62479340     0.103781228       0.3226632
#> 420_Abcf1             0.62288573     0.272113625       0.7888888
#> 421_Cluap1            0.76670578     0.240172943       0.6464593
#> 422_Hnrnpc            0.87608828     0.942618066       0.5643064
#> 423_Ptpn1             0.97194962     0.806223317       0.9790592
#> 424_Myh9              0.52376249     0.035324687       0.1973920
#> 425_U2surp            0.88009809     0.589597408       0.7251536
#> 426_NA                0.63564757     0.067817361       0.3320750
#> 427_Arhgef40          1.00000000     1.000000000       1.0000000
#> 428_Chaf1b            0.88292911     0.312784861       0.5378301
#> 429_Prpf4b            0.54094413     0.326203498       0.6972062
#> 430_Epb41l2           0.53917994     0.072843266       0.6138471
#> 431_Eif3b             1.00000000     1.000000000       1.0000000
#> 432_Isyna1            0.53917994     0.650300774       0.3226632
#> 433_U2surp            0.92291089     0.123733874       0.2514204
#> 434_LOC100765020      0.63385599     0.798306192       0.9032349
#> 435_Arhgef6           0.52376249     0.919930357       0.8566426
#> 436_Ptpn1             0.93279831     0.231705463       0.8786434
#> 437_Prpf4b            0.93490147     0.368796409       0.8639562
#> 438_Rpl35a            0.64073199     0.272113625       0.7920882
#> 439_Prpf4b            1.00000000     1.000000000       1.0000000
#> 440_Zyx               0.93408439     0.957142372       0.8487388
#> 441_Dbn1              0.86772434     0.279979469       0.2456534
#> 442_Chaf1b            1.00000000     1.000000000       0.7685783
#> 443_LOC113834282      0.89431253     0.158835844       0.5853013
#> 444_Gpsm2             0.45454337     0.011705607       0.1351724
#> 445_LOC100757535      1.00000000     1.000000000       1.0000000
#> 446_Cfap410           0.51985473     0.051078092       0.5853013
#> 447_Epb41l2           0.66796513     0.343423518       0.7605719
#> 448_Ncbp1             1.00000000     1.000000000       1.0000000
#> 449_Pacsin1           1.00000000     1.000000000       1.0000000
#> 450_Cstf2             0.61888944     0.093865695       0.7251536
#> 451_LOC100769437      0.76257173     0.981399784       0.9790592
#> 452_eIF2aK2           0.92291089     0.872113131       0.9894765
#> 453_Kiaa1191          0.97194962     0.799298832       0.9998830
#> 454_Mepce             0.54325567     0.678847456       0.8594514
#> 455_Cbx8              0.34788533     0.011705607       0.4629072
#> 456_Eed               0.62288573     0.494439642       0.9790592
#> 457_Cdc42ep1          0.88292911     0.227624172       0.5515318
#> 458_Lrrfip2           1.00000000     1.000000000       1.0000000
#> 459_Pacsin1           0.61241133     0.536829699       0.9790592
#> 460_Gpatch4           0.73967958     0.167212385       0.8303413
#> 461_Plin4             1.00000000     1.000000000       1.0000000
#> 462_NA                0.66687948     0.044182198       0.6715757
#> 463_Snip1             0.94543672     0.329225607       0.9546590
#> 464_Cyld              0.97194962     0.312011245       0.8184153
#> 465_Plin4             0.99508555     0.796160277       0.7895709
#> 466_Twf1              0.86772434     1.000000000       1.0000000
#> 467_LOC113834282      0.86772434     1.000000000       1.0000000
#> 468_Snip1             0.72880676     0.521352863       0.7580342
#> 469_Ppp4r3a           0.90515009     0.833342900       0.9790592
#> 470_Psip1             0.97399727     0.455298716       0.6747728
#> 471_Dnajc5            1.00000000     1.000000000       1.0000000
#> 472_Phf8              0.93408439     0.738183804       0.9998830
#> 473_Bola1             0.53591151     0.269847805       0.9047093
#> 474_Cdc42ep1          1.00000000     0.839266628       0.8690306
#> 475_Eif4ebp2          0.63385599     0.282421060       0.9201656
#> 476_Prpf38b           0.99508555     0.330329679       0.4711825
#> 477_Klhl26            0.83989471     0.065943776       0.2456534
#> 478_Hsph1             1.00000000     1.000000000       1.0000000
#> 479_Snip1             0.97194962     0.938325214       0.9164056
#> 480_Caskin2           0.54325567     0.065565791       0.8480630
#> 481_Plpp6             0.57311992     0.249506635       0.7857770
#> 482_NA                0.61241133     0.061202699       0.2456534
#> 483_Mlh1              0.79604118     0.181894212       0.7038252
#> 484_Gys1              0.57155299     0.009592916       0.2456534
#> 485_Tfg               0.16916527     0.009592916       0.5477890
#> 486_Arhgef6           1.00000000     1.000000000       1.0000000
#> 487_Mphosph10         1.00000000     1.000000000       1.0000000
#> 488_Hoxc10            0.92291089     0.322882788       0.8476647
#> 489_LOC100759640      1.00000000     1.000000000       1.0000000
#> 490_Arhgef40          0.93975237     0.355297643       0.6260825
#> 491_Dnajc5            1.00000000     1.000000000       1.0000000
#> 492_Tbc1d23           1.00000000     1.000000000       1.0000000
#> 493_Ubxn1             0.99251976     0.864692664       0.9047093
#> 494_Rab1a             1.00000000     1.000000000       1.0000000
#> 495_Eif3b             1.00000000     1.000000000       1.0000000
#> 496_Tceal8            0.57155299     0.093254560       0.5515318
#> 497_Dlgap4            0.76257173     0.457565692       0.9316856
#> 498_Smim13            0.92364395     0.853078223       0.7605719
#> 499_NA                1.00000000     1.000000000       1.0000000
#> 500_Lrch4             1.00000000     1.000000000       1.0000000
#> 501_Bola1             1.00000000     1.000000000       1.0000000
#> 502_NA                0.87609268     0.452259090       0.5382447
#> 503_Ptpn14            0.84177776     0.492332140       0.7434850
#> 504_LOC100759640      0.78328468     0.606386718       0.7855578
#> 505_Rps10             0.98509912     0.713568546       0.9509986
#> 506_Top2b             1.00000000     1.000000000       1.0000000
#> 507_Ssr3              0.86774666     0.554812277       0.8181409
#> 508_Homer3            1.00000000     1.000000000       1.0000000
#> 509_Phf8              0.58745270     0.492326464       0.6592642
#> 510_LOC100767716      0.58496599     0.046804463       0.5055045
#> 511_Xpa               0.86863314     0.962032310       0.8711279
#> 512_H671_21690        0.54502481     0.108920594       0.3975070
#> 513_LOC100769471      0.90515009     0.487247401       0.8303413
#> 514_Gas2l1            0.54325567     0.359604268       0.8487388
#> 515_Luzp1             0.71607132     0.043582798       0.6260825
#> 516_Gpbp1             1.00000000     0.285818903       0.9905692
#> 517_Gatad2b           1.00000000     1.000000000       1.0000000
#> 518_Gys1              0.70880247     0.598849644       0.5435710
#> 519_Top2b             0.72880676     0.251648545       0.4164747
#> 520_LOC100757535      1.00000000     0.329225607       0.5535905
#> 521_Lpcat4            0.88292911     0.118235739       0.5341971
#> 522_Arhgef6           1.00000000     1.000000000       1.0000000
#> 523_Cavin3            0.93975237     0.378690117       0.7605719
#> 524_Gpatch4           0.54502481     0.915557817       0.4960784
#> 525_Prpf38b           0.60502460     0.178156017       0.5683353
#> 526_Timm8a            0.64073199     0.254215304       0.4959790
#> 527_Cavin3            0.52376249     0.036591030       0.3922119
#> 528_Mkrn2             0.64828195     0.161989284       0.5098348
#> 529_Oser1             0.78999981     0.211749285       0.6171612
#> 530_Gsk3b             0.64073199     0.254215304       0.4959790
#> 531_Eef1b2            1.00000000     1.000000000       1.0000000
#> 532_Ampd2             0.54094413     0.251648545       0.6532590
#> 533_Lrrfip2           0.47623282     0.598849644       0.8402942
#> 534_Ring1             0.59557209     0.018987223       0.2528074
#> 535_Rlim              0.47623282     0.027421699       0.2528074
#> 536_LOC100759640      0.89679055     0.101906724       0.5546358
#> 537_LOC100759640      0.99508555     0.441780823       0.4428999
#> 538_Atp5pf            0.69624464     0.390006870       0.7321951
#> 539_Max               0.79319446     0.165335918       0.9164056
#> 540_Bap1              0.78328468     1.000000000       1.0000000
#> 541_Nsfl1c            0.54820576     0.536829699       0.4959790
#> 542_Prpf4b            1.00000000     0.936428449       0.8184153
#> 543_LOC100757535      0.93408439     0.712149728       0.6068931
#> 544_Mtmr10            0.64828195     0.354124539       0.4010473
#> 545_Hoxc10            0.47504094     0.095752128       0.7658717
#> 546_Trim35            1.00000000     1.000000000       1.0000000
#> 547_Eif4ebp2          0.77929811     0.145603610       0.5683353
#> 548_Dlgap4            0.72880676     0.507654023       1.0000000
#> 549_Gys1              0.54502481     0.318408177       0.8594514
#> 550_Sgtb              0.71607132     0.223869493       0.9332729
#> 551_Eri2              1.00000000     1.000000000       0.8639562
#> 552_Ccnd3             0.93408439     1.000000000       1.0000000
#> 553_Smim13            0.61241133     0.825474804       0.9047093
#> 554_Snrk              0.54325567     0.195989817       0.4555868
#> 555_Caskin2           1.00000000     1.000000000       1.0000000
#> 556_Pdcd11            0.47504094     0.095752128       0.7658717
#> 557_Pgam5             0.88292911     0.319757459       0.6995470
#> 558_Mphosph10         1.00000000     1.000000000       1.0000000
#> 559_Mideas            1.00000000     1.000000000       1.0000000
#> 560_Top2b             0.60118421     0.190393736       0.9790592
#> 561_LOC100763014      0.92291089     0.181894212       0.3397790
#> 562_Snip1             1.00000000     1.000000000       1.0000000
#> 563_Ubxn1             0.88292911     0.755595789       0.9740321
#> 564_LOC100750407      0.60580980     0.015227573       0.3141525
#> 565_Morf4l2           1.00000000     1.000000000       1.0000000
#> 566_Ctdspl2           1.00000000     1.000000000       1.0000000
#> 567_Cwf19l1           1.00000000     1.000000000       1.0000000
#> 568_Eef1b2            1.00000000     1.000000000       1.0000000
#> 569_C1H12orf45        1.00000000     1.000000000       1.0000000
#> 570_Znf367            1.00000000     0.589597408       0.7154080
#> 571_Ankrd34a          0.54325567     0.118235739       0.7038252
#> 572_Mllt11            0.72613118     0.921066131       0.7858065
#> 573_LOC100774792      0.67424968     0.399639137       0.7309207
#> 574_NA                0.94427730     0.108920594       0.2515247
#> 575_Cbx8              0.70908163     0.598849644       0.6138471
#> 576_Bckdk             0.89639358     0.452259090       0.9790592
#> 577_Snip1             0.46033394     1.000000000       1.0000000
#> 578_Nsfl1c            0.97194962     0.839336901       0.8356342
#> 579_Gas2l1            1.00000000     1.000000000       1.0000000
#> 580_Nudc              0.90515009     0.936428449       0.9790592
#> 581_Epb41l2           0.54094413     0.058615300       0.2456534
#> 582_Mtmr6             0.70731481     1.000000000       1.0000000
#> 583_Znf668            0.53591151     0.047572498       0.6256015
#> 584_Hsph1             1.00000000     1.000000000       0.9905692
#> 585_LOC113834282      0.71607132     0.770566581       0.8594514
#> 586_Ctdspl2           0.54325567     0.072843266       0.7605719
#> 587_Foxf1             0.87608828     0.351779941       0.8480630
#> 588_Luzp1             0.99508555     0.272113625       0.7321951
#> 589_Xpa               0.63993993     0.799298832       1.0000000
#> 590_Psip1             0.51979551     0.030851148       0.5721803
#> 591_Rbm7              0.65406828     0.062795263       0.5341971
#> 592_Mtrex             0.97194962     0.625117006       0.6923074
#> 593_Arhgef40          0.54502481     0.231705463       0.6715757
#> 594_Plekho2           0.89679055     0.803553113       0.5965493
#> 595_Bckdk             0.64073199     0.660718983       0.9870435
#> 596_Dut               1.00000000     1.000000000       1.0000000
#> 597_Abcf1             1.00000000     1.000000000       1.0000000
#> 598_Txnl1             1.00000000     1.000000000       1.0000000
#> 599_Nudc              1.00000000     1.000000000       1.0000000
#> 600_Sh3gl1            0.47623282     0.028743673       0.3614870
#> 601_Gatad2b           0.60253962     0.279979469       0.6171612
#> 602_Homer3            0.66796513     0.438379999       0.8391381
#> 603_Septin6           0.47623282     0.100570871       0.6835401
#> 604_Smim13            0.99176908     0.834107317       0.6464593
#> 605_Arhgef40          0.95686779     0.880579860       0.7416436
#> 606_Rpl32             0.97878789     0.209461038       0.5956977
#> 607_Tomm34            0.71607132     0.359986277       0.8303413
#> 608_Mlh1              0.92364395     0.756572162       0.9790592
#> 609_Tbcc              0.71724978     0.990436120       0.5474512
#> 610_Eif3d             0.54325567     0.809465395       0.5148762
#> 611_Snrk              0.62288573     0.685883250       0.7858065
#> 612_Bckdk             1.00000000     1.000000000       1.0000000
#> 613_Wdr3              0.86427681     0.388957672       0.6464664
#> 614_LOC100757535      1.00000000     1.000000000       1.0000000
#> 615_Dlg1              0.86772434     0.532690647       0.5956977
#> 616_LOC100767716      0.62479340     0.322882788       0.6995470
#> 617_Hnrnpc            0.54502481     0.452259090       0.9998830
#> 618_Mphosph10         0.05820201     0.103875602       0.8594514
#> 619_Eif3b             0.93490147     0.707124470       0.8367759
#> 620_Emd               1.00000000     1.000000000       1.0000000
#> 621_Txlng             0.71607132     0.860650480       0.6138471
#> 622_Prpf4b            1.00000000     1.000000000       1.0000000
#> 623_Rlim              0.79319446     1.000000000       1.0000000
#> 624_Eef1b2            1.00000000     1.000000000       0.9998830
#> 625_Def6              0.52376249     0.364712626       0.9869263
#> 626_LOC100765020      0.64552746     0.623579485       0.9998830
#> 627_U2surp            0.75854480     0.995949366       0.7901560
#> 628_Elf2              0.54094413     0.011705607       0.2456534
#> 629_Slc1a5            0.70908163     0.124564801       0.5098348
#> 630_NA                0.45904970     0.367092213       0.9058986
#> 631_Tfg               1.00000000     1.000000000       1.0000000
#> 632_Top2b             0.93408439     0.254215304       0.3756990
#> 633_Pip4p2            0.99176908     0.585388668       0.4892722
#> 634_Cdc42ep1          1.00000000     1.000000000       1.0000000
#> 635_Hsph1             1.00000000     1.000000000       1.0000000
#> 636_Twf1              1.00000000     1.000000000       1.0000000
#> 637_Nbn               0.71466746     0.317281054       0.6464664
#> 638_Psmd4             0.68196409     0.756572162       0.9795073
#> 639_Bap1              0.89639358     0.188600072       0.6860277
#> 640_Mepce             0.93408439     0.368796409       0.4398391
#> 641_Mideas            0.89047058     0.093583162       0.2456534
#> 642_LOC100759640      0.54325567     0.081632755       0.5214596
#> 643_Epb41l2           1.00000000     1.000000000       1.0000000
#> 644_Sav1              0.75827918     0.926824194       0.7888888
#> 645_Prpf4b            0.84802386     0.219937022       0.7305634
#> 646_Gnas              0.24835879     0.048653042       0.6110255
#> 647_Mllt1             1.00000000     1.000000000       1.0000000
#> 648_Poldip3           0.47623282     0.044152701       0.3975070
#> 649_Aldoa             0.14642828     0.009592916       0.2762895
#> 650_Rbbp8             0.46779026     0.018987223       0.2788533
#> 651_LOC113834282      0.71607132     0.072705811       0.2528074
#> 652_Gys1              0.94427730     0.752570747       0.4959790
#> 653_Hnrnpc            0.54502481     0.099296146       0.7423490
#> 654_Vps35             0.90515009     0.612667440       0.8594514
#> 655_Miga2             0.82187420     0.339218607       0.6093150
#> 656_Epb41l2           0.85192356     0.329225607       0.7321951
#> 657_Tob2              1.00000000     1.000000000       1.0000000
#> 658_Lamtor1           1.00000000     1.000000000       1.0000000
#> 659_LOC100759640      1.00000000     1.000000000       1.0000000
#> 660_Epb41l2           1.00000000     0.023674400       0.3175328
#> 661_Rlim              0.49554493     0.011705607       0.2657991
#> 662_Gys1              0.26292836     1.000000000       1.0000000
#> 663_LOC100750437      0.77895142     0.224001282       0.2762895
#> 664_NA                0.99983034     0.170774257       0.9998830
#> 665_Nbn               0.98602325     0.058324968       0.4555868
#> 666_Tyw3              1.00000000     1.000000000       1.0000000
#> 667_Gas2l1            0.92291089     0.752570747       0.8391381
#> 668_Fus               1.00000000     1.000000000       1.0000000
#> 669_Prpf38b           0.53917994     0.487247401       0.5683353
#> 670_Calu              1.00000000     1.000000000       0.6138471
#> 671_Rras2             0.54502481     0.272113625       0.9130044
#> 672_Prpf4b            0.81804469     0.860650480       0.9661887
#> 673_Nelfa             0.54325567     0.009592916       0.1865649
#> 674_LOC100754077      0.54094413     0.009592916       0.1865649
#> 675_Rbm28             0.98100329     0.874141663       0.6464593
#> 676_Nsfl1c            0.66930492     0.030810786       0.3141525
#> 677_Rnf126            0.89047058     0.161989284       0.6260825
#> 678_Eme1              0.64073199     0.616363428       0.5956977
#> 679_Nbn               0.82932424     0.044182198       0.2456534
#> 680_Eif4ebp2          0.62288573     0.044152701       0.2792817
#> 681_Wee1              1.00000000     1.000000000       1.0000000
#> 682_Prpf38b           1.00000000     0.569001074       0.9407237
#> 683_Luzp1             0.60118421     0.251648545       0.8480630
#> 684_Gas2l1            1.00000000     1.000000000       1.0000000
#> 685_Pdcd11            1.00000000     1.000000000       1.0000000
#> 686_Chaf1b            0.66796513     0.041040655       0.4306813
#> 687_Pycr1             0.47623282     0.033239223       0.4586578
#> 688_Phf8              1.00000000     0.021259362       0.1973920
#> 689_Raver1            1.00000000     1.000000000       1.0000000
#> 690_Dbn1              0.62388094     0.015227573       0.4959790
#> 691_Dut               0.46779026     0.078754217       0.3397691
#> 692_Prpf4b            0.97878789     0.589597408       1.0000000
#> 693_Prpf4b            0.90515009     0.885437511       0.9672425
#> 694_Efs               0.99508555     0.219937022       0.6532590
#> 695_NA                0.82561145     0.113354095       0.3071868
#> 696_Ppp2r5b           1.00000000     1.000000000       1.0000000
#> 697_Caskin2           0.99508555     0.727132902       0.9790592
#> 698_Arhgef40          0.72399299     0.769933446       0.9998830
#> 699_Zyx               1.00000000     1.000000000       1.0000000
#> 700_Mphosph10         0.64828195     0.532733162       0.6532590
#> 701_LOC113833392      0.49486143     0.080714181       0.5382447
#> 702_Cdc42ep1          0.52376249     0.004123866       0.1351724
#> 703_Snrpa1            0.67549601     1.000000000       1.0000000
#> 704_Ncbp1             0.88292911     0.990436120       0.9790592
#> 705_Gas2l1            1.00000000     1.000000000       1.0000000
#> 706_Gas2l1            0.14642828     0.004123866       0.2456534
#> 707_Bap1              0.94276860     0.343423518       0.6464664
#> 708_LOC100759640      0.64828195     0.100570871       0.6068931
#> 709_Cherp             1.00000000     1.000000000       1.0000000
#> 710_Nbn               0.64828195     0.080948689       0.2792817
#> 711_LOC100759640      1.00000000     1.000000000       1.0000000
#> 712_NA                1.00000000     1.000000000       1.0000000
#> 713_Eif3b             1.00000000     1.000000000       1.0000000
#> 714_Miga2             0.87581193     0.616363428       0.6563396
#> 715_Prpf4b            0.52376249     0.616363428       0.8939102
#> 716_Dbn1              0.71724978     0.635079569       0.9090431
#> 717_Ppp2r5b           0.47623282     0.007058261       0.2456534
#> 718_Exosc9            1.00000000     1.000000000       1.0000000
#> 719_Eif3b             0.62388094     0.114460768       0.9790592
#> 720_Ripk2             1.00000000     1.000000000       1.0000000
#> 721_Dlg1              0.65406828     0.158394831       0.6464664
#> 722_N4bp1             0.83654835     0.319858659       0.8487388
#> 723_Nudc              0.94792345     0.100570871       0.4010473
#> 724_Znf367            0.62288573     0.015227573       0.6464664
#> 725_Ring1             0.62388094     0.839109126       0.5214596
#> 726_Snrpa1            0.62945080     0.737458633       0.4714871
#> 727_U2surp            0.45454337     0.080268242       0.4555868
#> 728_LOC100764225      0.29488064     0.018737895       0.5683353
#> 729_Cdc42ep1          1.00000000     1.000000000       1.0000000
#> 730_Znf385a           1.00000000     1.000000000       1.0000000
#> 731_Ints1             1.00000000     1.000000000       1.0000000
#> 732_LOC113833392      1.00000000     1.000000000       1.0000000
#> 733_Lrch4             0.86639744     0.756572162       0.9796476
#> 734_Ctdspl2           0.54094413     0.268296359       0.7027710
#> 735_Prpf4b            0.97399727     0.834331354       0.8480630
#> 736_Luzp1             1.00000000     1.000000000       1.0000000
#> 737_Eif3b             1.00000000     1.000000000       1.0000000
#> 738_Ptpn14            1.00000000     1.000000000       1.0000000
#> 739_Rrp1              0.69625731     0.319757459       0.7305634
#> 740_Lrrfip2           0.66796513     1.000000000       1.0000000
#> 741_Nsfl1c            1.00000000     1.000000000       1.0000000
#> 742_Ddx51             0.97194962     0.291608060       0.7245342
#> 743_Prpf38b           0.71607132     0.589597408       0.8184153
#> 744_Eef1b2            0.97878789     0.885437511       0.5535905
#> 745_Znf385a           1.00000000     1.000000000       1.0000000
#> 746_Map9              0.53917994     0.081632755       0.8129386
#> 747_Rflnb             0.53917994     0.081632755       0.8129386
#> 748_NA                0.65406828     0.588411064       0.3975070
#> 749_C1H12orf45        0.65406828     0.588411064       0.3975070
#> 750_U2surp            0.82667600     0.180007676       0.4629072
#> 751_Caskin2           0.47623282     0.009592916       0.2528074
#> 752_Eri1              0.71633877     0.874141663       0.7044591
#> 753_Gsk3b             0.72880676     0.848200886       0.7660779
#> 754_LOC100766946      0.52376249     0.028456138       0.6530086
#> 755_Cnpy3             1.00000000     1.000000000       0.5214596
#> 756_Hnrnpc            1.00000000     1.000000000       0.7321951
#> 757_Ptpn14            1.00000000     1.000000000       0.7321951
#> 758_Slc7a11           0.63385599     0.018987223       0.2528074
#> 759_Hnrnpc            0.71633877     0.065157196       0.3648974
#> 760_Cdc37l1           0.97194962     0.663905553       0.8303413
#> 761_LOC100768405      0.93408439     0.554812277       0.7305634
#> 762_Rragc             0.93408439     0.554812277       0.7305634
#> 763_LOC113834282      0.46779026     0.030685499       0.3779698
#> 764_Fus               0.62288573     0.799298832       0.6171612
#> 765_Ubxn1             0.45454337     0.011705607       0.3975070
#> 766_Mmut              0.93408439     0.209461038       0.4144947
#> 767_Pdcd11            0.88924209     0.499536019       0.5554079
#> 768_LOC100757535      1.00000000     1.000000000       1.0000000
#> 769_Eif3b             0.84802386     0.721771211       0.9790592
#> 770_Rnf113a           1.00000000     1.000000000       1.0000000
#> 771_Sytl4             1.00000000     1.000000000       1.0000000
#> 772_Tlnrd1            0.54502481     0.417597821       0.7321951
#> 773_H671_1g1131       0.92496501     0.809738231       0.9657031
#> 774_Neurl1            0.88292911     0.838455031       0.9755035
#> 775_Zyx               0.52376249     0.108821927       0.2528074
#> 776_Ctdspl2           0.70908163     0.158835844       0.7201570
#> 777_Chaf1b            0.89867219     0.350582445       0.5148762
#> 778_Rragc             1.00000000     1.000000000       1.0000000
#> 779_Srfbp1            0.65406828     0.815429174       0.9316856
#> 780_Gys1              1.00000000     0.822406373       0.6974470
#> 781_Usp15             0.93408439     0.901620154       0.8452398
#> 782_Arhgef40          0.97935078     0.452259090       0.4873530
#> 783_Gigyf1            1.00000000     0.368796409       0.7787875
#> 784_Minar1            0.57451925     0.124140958       0.2792817
#> 785_Dus2              1.00000000     1.000000000       1.0000000
#> 786_Gatad2b           0.72613118     0.813074192       0.7564088
#> 787_Eif5              0.94427730     0.985476771       0.9740321
#> 788_Epb41l2           0.57448002     0.255538068       0.9790592
#> 789_Arl6ip4           0.87608828     0.721593970       0.9998830
#> 790_Plin4             0.88292911     0.368796409       0.9918682
#> 791_Elf2              0.99275666     0.401225862       0.5535905
#> 792_Plin4             0.97399727     0.342023429       1.0000000
#> 793_Snip1             1.00000000     1.000000000       1.0000000
#> 794_Txlng             1.00000000     1.000000000       1.0000000
#> 795_LOC100769437      1.00000000     0.996871421       0.9788357
#> 796_Caskin2           0.67729723     0.322882788       0.6434996
#> 797_NA                0.24835879     0.025997938       0.3346013
#> 798_Synm              1.00000000     1.000000000       1.0000000
#> 799_Synm              0.97194962     0.295813395       0.8184153
#> 800_Ube2c             1.00000000     0.609625186       0.9905692
#> 801_Sgtb              0.66796513     0.082558457       0.5430705
#> 802_Prpf4b            0.58745270     0.797826368       0.9585265
#> 803_Epb41l2           0.14642828     0.161679006       0.9740321
#> 804_Mllt1             0.84802386     0.834107317       0.7776935
#> 805_LOC100759640      1.00000000     1.000000000       1.0000000
#> 806_Epb41l2           0.63851869     0.633144523       0.8786434
#> 807_Znf280b           1.00000000     0.388957672       0.8452398
#> 808_Kiaa1143          1.00000000     1.000000000       1.0000000
#> 809_Gas2l1            1.00000000     1.000000000       1.0000000
#> 810_Srp72             0.64828195     0.359667514       0.8939102
#> 811_Tomm22            0.93408439     0.422286817       0.4549154
#> 812_Psip1             0.93975237     0.071271034       0.3036211
#> 813_Arhgef37          0.92291089     0.632610498       0.9905692
#> 814_Bckdk             0.93975237     0.804749173       0.9449472
#> 815_Strip1            0.47504094     0.009592916       0.2762895
#> 816_Usp15             1.00000000     1.000000000       1.0000000
#> 817_Ssr3              0.47623282     0.616363428       0.5341971
#> 818_Strip1            0.54094413     0.141221624       0.5515318
#> 819_Eif3b             0.32526276     0.008305263       0.2514204
#> 820_U2surp            0.93408439     0.076614751       0.2528074
#> 821_Bend3             0.64073199     0.018987223       0.3320750
#> 822_Rps10             0.82649872     0.634193211       0.5956977
#> 823_Rpl23a            0.71633877     0.700740104       0.7400550
#> 824_Nbn               0.97194962     0.848200886       0.5515318
#> 825_Rpap3             0.72613118     0.316008364       0.6563396
#> 826_LOC100759640      1.00000000     1.000000000       1.0000000
#> 827_Ric8a             0.94792345     0.690767762       0.7073961
#> 828_Hsph1             0.67729723     0.551322191       0.8487388
#> 829_LOC100759640      0.66796513     0.137187744       0.3175328
#> 830_LOC100757535      0.89639358     0.838455031       0.7656499
#> 831_Gigyf1            0.94427730     0.195989817       0.9790592
#> 832_Dbn1              0.71724978     0.105583422       0.6260825
#> 833_Snrk              0.88292911     0.119427976       0.6715757
#> 834_Prpf38b           0.54094413     0.246083797       0.9728359
#> 835_LOC100766868      0.05820201     0.004123866       0.6171612
#> 836_LOC100766868      0.87530023     0.492052392       0.8993151
#> 837_Wbp11             0.93408439     0.752570747       0.7888888
#> 838_Rusc2             0.90515009     0.044182198       0.1316251
#> 839_Eif3b             0.62288573     0.378206060       0.4909478
#> 840_Ptpn14            0.72880676     0.492052392       0.9606716
#> 841_Rlim              0.83839993     0.941385992       0.9332729
#> 842_Ints1             0.64828195     1.000000000       1.0000000
#> 843_Chaf1b            0.54502481     0.598849644       0.9790592
#> 844_Dlg1              1.00000000     1.000000000       1.0000000
#> 845_Lamtor1           0.56985618     0.028561245       0.2528074
#> 846_Tab1              0.65194228     0.044152701       0.3226632
#> 847_Dbn1              1.00000000     1.000000000       1.0000000
#> 848_Psip1             1.00000000     1.000000000       1.0000000
#> 849_Dbn1              0.66796513     0.541728497       0.9332729
#> 850_Pabpc1            0.76257173     0.339218607       0.7714479
#> 851_Hnrnpc            0.91576521     0.158835844       0.2528074
#> 852_Emd               0.54325567     0.025997938       0.1865649
#> 853_LOC100764225      0.54325567     0.083810939       0.4412340
#> 854_Nup50             0.90515009     0.322882788       0.6013451
#> 855_Ctcf              0.71724978     0.124065223       0.6138471
#> 856_Raly              0.94792345     0.246900570       0.3966787
#> 857_Bard1             0.45454337     0.022880728       0.4586578
#> 858_Ptpn14            0.86863314     0.460414349       0.8258043
#> 859_LOC100757535      1.00000000     0.279979469       0.8049503
#> 860_Psmd2             0.88292911     0.640060350       1.0000000
#> 861_Junb              0.97399727     0.532733162       0.5515318
#> 862_C1qbp             0.99983034     0.804069909       0.9790592
#> 863_Lrch4             0.97878789     0.045386395       0.2456534
#> 864_CUNH14orf93       0.61241133     0.591600668       1.0000000
#> 865_U2surp            0.76257173     0.788096638       0.5535905
#> 866_Raly              1.00000000     1.000000000       1.0000000
#> 867_LOC100774417      0.93408439     0.187920561       0.3226632
#> 868_Srp72             0.98602325     0.860650480       0.8402942
#> 869_LOC100764225      1.00000000     1.000000000       1.0000000
#> 870_Morf4l2           0.88137715     0.809465395       0.6485884
#> 871_CUNH9orf40        0.66005075     0.079545642       0.4808161
#> 872_Gas2l1            0.16916527     0.004123866       0.1316251
#> 873_Atp5pf            0.72251712     0.099723851       0.2841615
#> 874_Lrrfip2           0.64073199     0.116295046       0.5683353
#> 875_Prpf4b            1.00000000     1.000000000       1.0000000
#> 876_Top2b             1.00000000     1.000000000       1.0000000
#> 877_Mepce             0.93975237     0.958954018       0.8711279
#> 878_Ptpn14            0.72880676     0.145603610       0.9332729
#> 879_Dnajc25           0.97194962     0.998833004       0.7743909
#> 880_Cbx8              0.62288573     0.373980374       0.9998830
#> 881_Synm              0.65406828     0.220944799       0.6171612
#> 882_Def6              0.53793240     0.241416749       0.6260825
#> 883_Gys1              0.89639358     0.462476170       0.8303413
#> 884_Luzp1             0.99508555     0.532733162       0.7012143
#> 885_Synm              0.79677029     0.566071310       0.9740321
#> 886_Snip1             1.00000000     1.000000000       1.0000000
#> 887_Top2b             0.47623282     0.005603759       0.1351724
#> 888_NA                0.54094413     0.113354095       0.8594514
#> 889_Trim35            0.72613118     0.981399784       0.8391381
#> 890_Znf385a           0.54325567     0.093254560       0.7337186
#> 891_Chaf1b            0.93408439     0.756572162       0.5683353
#> 892_Abcf1             1.00000000     0.269441152       0.6260825
#> 893_Pdcd11            1.00000000     1.000000000       1.0000000
#> 894_Dlg1              0.70908163     0.634718999       0.7321951
#> 895_Dbn1              1.00000000     1.000000000       1.0000000
#> 896_LOC100752363      0.79677029     0.113354095       0.5341971
#> 897_Ppp4r3a           1.00000000     1.000000000       1.0000000
#> 898_Gas2l1            0.97194962     0.227624172       0.4555868
#> 899_Mtmr10            1.00000000     1.000000000       1.0000000
#> 900_Cyld              0.54502481     0.087519180       0.4808161
#> 901_NA                0.54325567     0.048653042       0.4157832
#> 902_Rnf113a           1.00000000     1.000000000       1.0000000
#> 903_Nelfa             1.00000000     1.000000000       1.0000000
#> 904_Zkscan1           1.00000000     1.000000000       1.0000000
#> 905_Chaf1b            0.94013089     0.859444878       0.7685402
#> 906_Eif3b             0.64073199     0.848200886       0.6592642
#> 907_Top2b             0.52376249     0.252674563       0.9998830
#> 908_Chaf1b            1.00000000     1.000000000       0.7415917
#> 909_Epb41l2           1.00000000     1.000000000       1.0000000
#> 910_C3H11orf58        1.00000000     0.722042209       0.8402942
#> 911_Top2b             0.64828195     0.017767918       0.1865649
#> 912_Wee1              0.99176908     0.399784993       0.4258505
#> 913_Raly              0.64828195     0.756572162       0.7479895
#> 914_H671_1g2680       0.87530023     0.074667937       0.1973920
#> 915_Eef1b2            1.00000000     1.000000000       1.0000000
#> 916_Gas2l1            0.54325567     1.000000000       1.0000000
#> 917_Epb41l2           0.71607132     0.460414349       0.8639562
#> 918_Rpl23a            0.71607132     0.011705607       0.1316251
#> 919_Chmp2b            0.86772434     0.336096117       0.3779698
#> 920_Lrrfip2           0.95524306     0.553896530       0.7564088
#> 921_Aldoa             0.63454585     0.786737186       0.6260825
#> 922_Cby1              0.97935078     0.462802123       0.6485884
#> 923_LOC100759640      0.88292911     0.187885817       0.2528074
#> 924_Rbm28             0.47623282     0.337140363       0.7073961
#> 925_Skiv2l            0.62288573     0.460926236       0.9914096
#> 926_Ints1             0.97194962     0.243179869       0.3779698
#> 927_Ehd1              0.62388094     0.913964935       0.9905692
#> 928_Nr2f6             0.54325567     0.011493459       0.1772224
#> 929_Top2b             0.86062850     0.240009924       0.4711825
#> 930_Lrrfip2           0.70908163     0.099723851       0.3975070
#> 931_Pip4p2            0.72880676     0.250783169       0.9332729
#> 932_Srp72             0.66687948     0.495954955       0.8049503
#> 933_Mtmr9             0.53917994     0.067817361       0.2515247
#> 934_Gigyf1            0.97878789     0.796500636       0.7079214
#> 935_Rbm7              0.86772434     0.205373417       0.7305634
#> 936_LOC100773565      0.52376249     0.080268242       0.4469415
#> 937_Trim35            0.72613118     0.835850865       0.7166115
#> 938_Cbx8              0.79319446     0.583764133       0.3320750
#> 939_Rplp0             0.52376249     0.219937022       0.5602351
#> 940_Aldoa             0.05820201     0.144987365       0.9905692
#> 941_NA                0.71724978     0.011493459       0.1865649
#> 942_Zyx               0.54502481     0.024200859       0.1316251
#> 943_Psip1             0.70026237     0.025997938       0.2528074
#> 944_Slc7a11           0.52376249     0.011705607       0.1973920
#> 945_Miga2             0.86772434     0.023674400       0.1316251
#> 946_Arhgef6           1.00000000     1.000000000       1.0000000
#> 947_Dlgap4            0.65384212     0.152487189       0.6995470
#> 948_Ampd2             1.00000000     1.000000000       1.0000000
#> 949_Luzp1             0.60118421     0.756572162       0.8095756
#> 950_Camlg             0.92291089     0.368330199       0.6923074
#> 951_Pfkfb3            1.00000000     1.000000000       1.0000000
#> 952_NA                0.75854480     0.296230394       0.8129386
#> 953_Raly              0.93975237     0.603016843       0.8939102
#> 954_Kiaa1143          0.78315628     0.429382386       0.7079214
#> 955_Bcar1             0.46779026     0.011705607       0.3779698
#> 956_Gatad2b           1.00000000     1.000000000       1.0000000
#> 957_Eif4ebp2          0.90069073     0.676683041       0.6563396
#> 958_Fam76b            1.00000000     1.000000000       1.0000000
#> 959_Camlg             0.24684158     0.018982808       0.5956977
#> 960_LOC100754077      1.00000000     1.000000000       0.9998830
#> 961_NA                0.86427681     0.769933446       0.5853013
#> 962_Epb41l2           0.52376249     0.605252603       0.8711279
#> 963_Ankrd34a          0.62479340     0.072705811       0.2792817
#> 964_Zc3h15            0.93975237     0.553896530       0.7660779
#> 965_Def6              0.93975237     0.247766927       0.7816251
#> 966_Srsf6             0.59313749     0.553896530       1.0000000
#> 967_H671_4g11480      0.89639358     0.249506635       0.5341971
#> 968_Top2b             0.97194962     0.839266628       0.8997262
#> 969_LOC100769471      0.46779026     0.268195707       0.9661887
#> 970_Raver1            0.62288573     0.589597408       0.8049503
#> 971_Etv3              0.93408439     0.020193592       0.2528074
#> 972_Psd               0.89679055     0.521352863       0.4950697
#> 973_Usp15             1.00000000     1.000000000       0.8480630
#> 974_Nol7              1.00000000     1.000000000       0.9740321
#> 975_Stk38             0.14868675     0.025997938       0.5757619
#> 976_Smim13            1.00000000     1.000000000       1.0000000
#> 977_Etv3              0.29488064     0.009592916       0.2792817
#> 978_Synm              0.76257173     0.158835844       0.5214596
#> 979_Pwp1              0.98684650     0.650300774       0.8613130
#> 980_Fus               0.89639358     0.834331354       0.8129386
#> 981_Junb              0.92291089     0.839266628       0.7605719
#> 982_Phf8              0.89410612     0.913964935       1.0000000
#> 983_Nelfa             1.00000000     1.000000000       1.0000000
#> 984_Prpf4b            1.00000000     1.000000000       1.0000000
#> 985_Abraxas1          0.47623282     0.009592916       0.2515247
#> 986_Prpf4b            0.84802386     0.352472407       0.8258043
#> 987_Raver1            0.81948788     0.589597408       0.8402942
#> 988_Caap1             0.72815581     0.937520812       0.8184153
#> 989_Rpap3             1.00000000     1.000000000       1.0000000
#> 990_Hsph1             0.89639358     0.813074192       0.3175328
#> 991_LOC100750437      0.97194962     0.399784993       1.0000000
#> 992_Mepce             0.58980148     0.233136373       0.9790592
#> 993_Efs               0.89639358     0.452259090       0.6995470
#> 994_Epb41l2           0.13952113     0.009592916       0.2528074
#> 995_Abcf1             0.99508555     0.438379999       0.9790592
#> 996_NA                0.73942586     0.269847805       0.9790592
#> 997_Eif4ebp2          0.79604118     0.590957303       0.9790592
#> 998_Pfkfb3            0.71607132     0.319858659       0.6236418
#> 999_Hnrnpc            0.78328468     0.913964935       0.9816002
#> 1000_Psmd2            0.85570913     0.752570747       0.8997262
#>                  120_vs_neighbors
#> 1_Top2b                 0.9974120
#> 2_NA                    1.0000000
#> 3_Snip1                 0.9974120
#> 4_Tomm34                0.9974120
#> 5_Pus3                  1.0000000
#> 6_Ints1                 0.9974120
#> 7_Mlh1                  0.9974120
#> 8_LOC100750437          0.9974120
#> 9_Pabpc1                0.9974120
#> 10_Top2b                0.9974120
#> 11_Gorasp1              0.9974120
#> 12_Ints1                1.0000000
#> 13_Syvn1                0.9974120
#> 14_Znf280b              0.9974120
#> 15_Mrnip                0.9974120
#> 16_Rragc                0.9974120
#> 17_Gorasp1              0.9974120
#> 18_Tomm34               0.9974120
#> 19_LOC100757430         0.9974120
#> 20_Ubxn1                0.9974120
#> 21_H671_1g1131          0.9974120
#> 22_Luzp1                0.9974120
#> 23_Efs                  0.9974120
#> 24_Mta2                 0.9974120
#> 25_Nedd1                0.9974120
#> 26_Gigyf1               0.9974120
#> 27_Myh9                 0.9974120
#> 28_Caskin2              0.9974120
#> 29_Papolg               0.9974120
#> 30_Tfg                  0.9974120
#> 31_Rpl34                0.9974120
#> 32_Mideas               0.9974120
#> 33_Gys1                 0.9974120
#> 34_Arhgef6              0.9974120
#> 35_Ctdspl2              1.0000000
#> 36_Ptpn14               0.9974120
#> 37_Raly                 0.9974120
#> 38_Znhit3               0.9974120
#> 39_LOC113833392         1.0000000
#> 40_Luc7l3               0.9974120
#> 41_Rplp0                1.0000000
#> 42_Gys1                 0.9974120
#> 43_Rpl22l1              0.9974120
#> 44_Eif3b                0.9974120
#> 45_Med26                0.9974120
#> 46_Mepce                0.9974120
#> 47_Pdcd11               0.9974120
#> 48_Twf1                 0.9974120
#> 49_LOC100759640         1.0000000
#> 50_Wrnip1               0.9974120
#> 51_Poldip3              0.9974120
#> 52_Ampd2                0.9974120
#> 53_Mea1                 0.9974120
#> 54_Dbn1                 0.9974120
#> 55_Snip1                0.9974120
#> 56_Srsf6                0.9974120
#> 57_LOC113834282         0.9974120
#> 58_Map9                 0.9974120
#> 59_Cdc42ep1             0.9974120
#> 60_Poldip3              1.0000000
#> 61_LOC100764225         1.0000000
#> 62_Epb41l2              0.9974120
#> 63_H671_4g11480         0.9974120
#> 64_Nbn                  0.9974120
#> 65_U2surp               0.9974120
#> 66_Gigyf1               0.9974120
#> 67_NA                   0.9974120
#> 68_Luc7l3               1.0000000
#> 69_LOC100752363         0.9974120
#> 70_Ampd2                1.0000000
#> 71_LOC100759640         1.0000000
#> 72_Stam                 0.9974120
#> 73_Nsfl1c               0.9974120
#> 74_Pfkfb3               0.9974120
#> 75_Rad23a               0.9974120
#> 76_Elf2                 0.9974120
#> 77_Crem                 0.9974120
#> 78_Rragc                0.9974120
#> 79_Lrrfip2              0.9974120
#> 80_Zyx                  0.9974120
#> 81_Lrrfip2              0.9974120
#> 82_Gatad2b              0.9995772
#> 83_Bcar1                0.9974120
#> 84_Ehd1                 0.9974120
#> 85_LOC113834282         0.9974120
#> 86_Tmem230              0.9974120
#> 87_Ncbp1                0.9974120
#> 88_Mllt1                1.0000000
#> 89_Stk17b               0.9974120
#> 90_Dlgap4               0.9974120
#> 91_Papolg               0.9974120
#> 92_Cyld                 0.9974120
#> 93_Gigyf1               0.9974120
#> 94_Lrrfip2              1.0000000
#> 95_Lrrfip2              0.9974120
#> 96_Rlim                 0.9974120
#> 97_Eif3b                0.9974120
#> 98_Mphosph10            1.0000000
#> 99_Gatad2b              0.9974120
#> 100_Srsf6               0.9974120
#> 101_Zyx                 0.9974120
#> 102_Mphosph10           1.0000000
#> 103_Psip1               0.9974120
#> 104_Fbl                 0.9974120
#> 105_H671_1g2680         0.9974120
#> 106_Sgtb                0.9974120
#> 107_Gnl3                0.9974120
#> 108_Eif3b               0.9974120
#> 109_Serpinb1            1.0000000
#> 110_N4bp1               0.9974120
#> 111_Snip1               1.0000000
#> 112_Psip1               1.0000000
#> 113_Mlh1                0.9974120
#> 114_Bsg                 0.9974120
#> 115_Tnpo1               0.9974120
#> 116_H671_1g2680         1.0000000
#> 117_Cbx8                1.0000000
#> 118_Mideas              1.0000000
#> 119_Mideas              0.9974120
#> 120_Dcun1d3             0.9974120
#> 121_Dlg1                1.0000000
#> 122_Rad23a              0.9974120
#> 123_Srsf6               0.9974120
#> 124_Stx7                0.9995772
#> 125_Pdcd11              0.9974120
#> 126_Kiaa1958            0.9995772
#> 127_Pwp1                0.9974120
#> 128_Txlng               0.9974120
#> 129_Junb                1.0000000
#> 130_LOC100759640        1.0000000
#> 131_Dbn1                0.9974120
#> 132_Top2b               0.9974120
#> 133_Rusc2               0.9974120
#> 134_NA                  1.0000000
#> 135_LOC113837251        0.9974120
#> 136_Fam76b              1.0000000
#> 137_Ptpn14              0.9974120
#> 138_Chmp4b              0.9974120
#> 139_Prpf4b              1.0000000
#> 140_Eif3b               1.0000000
#> 141_Nsfl1c              1.0000000
#> 142_Pdlim7              1.0000000
#> 143_Rnf113a             1.0000000
#> 144_Epb41l2             0.9974120
#> 145_Hnrnpc              1.0000000
#> 146_LOC113834282        0.9974120
#> 147_Plekho2             0.9974120
#> 148_Med26               0.9974120
#> 149_Arhgef40            0.9974120
#> 150_NA                  0.9974120
#> 151_Phf8                0.9974120
#> 152_Minar1              1.0000000
#> 153_H671_21690          0.9974120
#> 154_Arhgef40            0.9974120
#> 155_Chaf1b              0.9974120
#> 156_Prpf4b              0.9974120
#> 157_Znf367              0.9974120
#> 158_Luzp1               0.9974120
#> 159_LOC113833882        0.9974120
#> 160_Hnrnpc              0.9974120
#> 161_Mepce               0.9974120
#> 162_Ubxn1               0.9974120
#> 163_Mllt1               0.9974120
#> 164_Chaf1b              0.9974120
#> 165_Raly                0.9974120
#> 166_Gas2l1              0.9974120
#> 167_Dlg1                0.9974120
#> 168_Hoxc10              0.9974120
#> 169_Gigyf1              1.0000000
#> 170_Luzp1               0.9974120
#> 171_Srp72               0.9974120
#> 172_LOC100771461        1.0000000
#> 173_Chaf1b              0.9974120
#> 174_C3H11orf58          1.0000000
#> 175_Pdcd11              0.9974120
#> 176_Psip1               0.9974120
#> 177_Prpf4b              0.9974120
#> 178_Rnf113a             0.9974120
#> 179_Irf3                0.9974120
#> 180_Smim13              0.9974120
#> 181_Gnl3                0.9974120
#> 182_Psma5               0.9974120
#> 183_Ptpn14              0.9974120
#> 184_Prpf4b              1.0000000
#> 185_Top2b               1.0000000
#> 186_Prpf38b             0.9974120
#> 187_Epb41l2             0.9974120
#> 188_Eif3b               1.0000000
#> 189_Hnrnpc              0.9974120
#> 190_LOC100758278        0.9974120
#> 191_Prpf4b              0.9974120
#> 192_Caskin2             0.9974120
#> 193_LOC100752363        0.9974120
#> 194_Septin6             0.9974120
#> 195_Max                 0.9974120
#> 196_Mid1ip1             0.9974120
#> 197_NA                  0.9974120
#> 198_Hsph1               1.0000000
#> 199_Nol7                0.9974120
#> 200_Raly                0.9974120
#> 201_Smim13              0.9974120
#> 202_LOC100757535        0.9974120
#> 203_Net1                0.9974120
#> 204_LOC100754077        1.0000000
#> 205_Snip1               0.9974120
#> 206_Hnrnpc              0.9974120
#> 207_Ldlrap1             0.9974120
#> 208_Luzp1               0.9974120
#> 209_Rpl26               0.9974120
#> 210_Epb41l2             0.9974120
#> 211_Znf367              0.9974120
#> 212_Dlgap4              1.0000000
#> 213_Plekho2             1.0000000
#> 214_Zpr1                0.9974120
#> 215_Dlgap4              1.0000000
#> 216_Def6                0.9974120
#> 217_Eif4ebp2            0.9974120
#> 218_Eef1b2              0.9974120
#> 219_Rad23a              0.9974120
#> 220_Morf4l2             0.9974120
#> 221_Arhgef40            0.9974120
#> 222_NA                  0.9974120
#> 223_LOC100773565        0.9974120
#> 224_Dus2                0.9974120
#> 225_Pip4p2              0.9974120
#> 226_Top2b               0.9974120
#> 227_Znf280b             0.9974120
#> 228_Pdcd11              0.9974120
#> 229_Bckdk               1.0000000
#> 230_Arhgef40            1.0000000
#> 231_Mepce               1.0000000
#> 232_Ccnd3               0.9974120
#> 233_Phf8                0.9974120
#> 234_H671_1g2680         0.9974120
#> 235_Ell                 0.9974120
#> 236_U2surp              0.9974120
#> 237_Rps10               0.9974120
#> 238_Ctdspl2             0.9974120
#> 239_Top2b               0.9974120
#> 240_Msantd3             0.9974120
#> 241_Fam76b              0.9974120
#> 242_Ppp4r3a             0.9974120
#> 243_Gpatch4             0.9974120
#> 244_Nudc                0.9974120
#> 245_Nol7                0.9974120
#> 246_Plekho2             0.9974120
#> 247_Prpf4b              0.9974120
#> 248_Mta2                0.9974120
#> 249_U2surp              0.9974120
#> 250_Ubxn1               1.0000000
#> 251_Rlim                0.9974120
#> 252_Atat1               0.9974120
#> 253_Ubxn1               0.9995772
#> 254_H671_1g2680         0.9974120
#> 255_eIF2aK2             0.9974120
#> 256_Skiv2l              1.0000000
#> 257_Rpl28               1.0000000
#> 258_LOC100759640        1.0000000
#> 259_Gatad2b             0.9974120
#> 260_NA                  0.9974120
#> 261_Gprasp1             0.9974120
#> 262_Luzp1               1.0000000
#> 263_Slc1a5              1.0000000
#> 264_LOC113834282        0.9974120
#> 265_Srsf6               0.9974120
#> 266_Cdc42ep1            0.9974120
#> 267_Net1                0.9974120
#> 268_Caskin2             0.9974120
#> 269_LOC100759640        0.9974120
#> 270_Mideas              0.9974120
#> 271_Luzp1               0.9974120
#> 272_Emd                 1.0000000
#> 273_Plpp6               0.9974120
#> 274_LOC100759640        1.0000000
#> 275_Rps7                0.9974120
#> 276_Fkbp1a              0.9974120
#> 277_Gatad2b             0.9974120
#> 278_Znf385a             0.9974120
#> 279_Arhgef6             0.9974120
#> 280_Slirp               0.9974120
#> 281_Skiv2l              0.9974120
#> 282_H671_21690          0.9974120
#> 283_Kat8                0.9974120
#> 284_Nkap                0.9974120
#> 285_Gsk3b               0.9974120
#> 286_Ints1               0.9974120
#> 287_Gas2l1              1.0000000
#> 288_LOC100759640        1.0000000
#> 289_Top2b               0.9974120
#> 290_Kif20b              0.9974120
#> 291_Phf8                1.0000000
#> 292_Snip1               0.9974120
#> 293_Gsk3b               1.0000000
#> 294_Caskin2             1.0000000
#> 295_C3H11orf58          0.9974120
#> 296_Lrch4               0.9974120
#> 297_LOC113834282        0.9974120
#> 298_LOC100750407        1.0000000
#> 299_LOC113833392        1.0000000
#> 300_LOC113833882        0.9974120
#> 301_Ldlrap1             0.9974120
#> 302_Wee1                0.9974120
#> 303_Caap1               1.0000000
#> 304_Eif4ebp2            0.9974120
#> 305_Ripk2               0.9974120
#> 306_Srp72               0.9974120
#> 307_Taok2               0.9974120
#> 308_Nr2f6               0.9974120
#> 309_Arhgef40            0.9974120
#> 310_Gys1                1.0000000
#> 311_Dlg1                1.0000000
#> 312_Vapb                0.9974120
#> 313_LOC100757535        0.9974120
#> 314_Mkrn2               0.9974120
#> 315_Eif3b               0.9974120
#> 316_Isyna1              1.0000000
#> 317_Prpf4b              0.9974120
#> 318_LOC113833882        0.9974120
#> 319_Lrch4               0.9974120
#> 320_Dbn1                0.9974120
#> 321_Abcf1               0.9974120
#> 322_Ints1               0.9974120
#> 323_C3H11orf58          1.0000000
#> 324_Psma5               0.9974120
#> 325_Fundc1              0.9974120
#> 326_Papolg              1.0000000
#> 327_Mideas              0.9974120
#> 328_Ubxn1               0.9974120
#> 329_Synm                0.9974120
#> 330_Arhgef6             1.0000000
#> 331_Ptpn14              0.9974120
#> 332_Pgrmc1              0.9974120
#> 333_Myh9                0.9974120
#> 334_Etv3                1.0000000
#> 335_Ip6k1               1.0000000
#> 336_Luzp1               1.0000000
#> 337_Ptpn14              1.0000000
#> 338_Caskin2             1.0000000
#> 339_Chaf1b              0.9974120
#> 340_Ubxn1               1.0000000
#> 341_Ube2c               1.0000000
#> 342_Gins2               0.9974120
#> 343_Nlgn2               0.9974120
#> 344_Nf2                 0.9974120
#> 345_Pip4p2              0.9974120
#> 346_Emd                 0.9974120
#> 347_Top2b               0.9974120
#> 348_Trim35              0.9974120
#> 349_NA                  0.9974120
#> 350_NA                  1.0000000
#> 351_Mideas              0.9974120
#> 352_Gas2l1              0.9974120
#> 353_Ampd2               0.9974120
#> 354_Calu                1.0000000
#> 355_Fam76b              1.0000000
#> 356_Dlg1                0.9995772
#> 357_Srsf6               0.9974120
#> 358_Chaf1b              0.9974120
#> 359_Dbn1                0.9974120
#> 360_Tcf25               0.9974120
#> 361_Psip1               0.9974120
#> 362_Cnpy3               0.9974120
#> 363_LOC100759640        0.9974120
#> 364_Zyx                 0.9974120
#> 365_Lrch4               0.9974120
#> 366_Bola1               1.0000000
#> 367_Znf385a             0.9974120
#> 368_Kif20b              0.9974120
#> 369_Ell                 0.9974120
#> 370_Ell                 0.9974120
#> 371_Srsf6               0.9974120
#> 372_Pwp1                0.9974120
#> 373_Def6                0.9974120
#> 374_Cbx8                1.0000000
#> 375_Ddx51               0.9974120
#> 376_Psip1               0.9974120
#> 377_Arhgef40            0.9974120
#> 378_Raly                0.9974120
#> 379_NA                  0.9974120
#> 380_Lrrfip2             0.9974120
#> 381_Gnl3                0.9974120
#> 382_Caskin2             0.9974120
#> 383_Rragc               0.9974120
#> 384_Caskin2             0.9974120
#> 385_Bcar1               1.0000000
#> 386_Homer3              0.9974120
#> 387_Luzp1               0.9974120
#> 388_N4bp1               1.0000000
#> 389_Ppp4r3a             0.9974120
#> 390_H671_1g2680         0.9974120
#> 391_Gnl3                0.9974120
#> 392_Top2b               0.9974120
#> 393_Oser1               1.0000000
#> 394_Snrk                0.9974120
#> 395_Kat8                0.9974120
#> 396_Raver1              1.0000000
#> 397_Pdcd11              0.9974120
#> 398_Rps20               0.9974120
#> 399_Bsg                 0.9974120
#> 400_Raly                0.9974120
#> 401_Pdcd2               0.9974120
#> 402_Caskin2             0.9974120
#> 403_LOC100773571        0.9974120
#> 404_Papolg              1.0000000
#> 405_LOC100757535        0.9974120
#> 406_Caap1               0.9974120
#> 407_Psip1               0.9974120
#> 408_Dbn1                0.9974120
#> 409_Mta2                0.9995772
#> 410_Abcf1               1.0000000
#> 411_LOC100754108        0.9974120
#> 412_Slirp               0.9974120
#> 413_Nelfa               0.9974120
#> 414_Aggf1               0.9974120
#> 415_Bap1                0.9974120
#> 416_Luc7l3              0.9974120
#> 417_Rrp1                0.9974120
#> 418_Wrnip1              0.9974120
#> 419_NA                  0.9974120
#> 420_Abcf1               0.9974120
#> 421_Cluap1              0.9974120
#> 422_Hnrnpc              0.9974120
#> 423_Ptpn1               0.9974120
#> 424_Myh9                0.9974120
#> 425_U2surp              0.9974120
#> 426_NA                  0.9974120
#> 427_Arhgef40            1.0000000
#> 428_Chaf1b              0.9974120
#> 429_Prpf4b              0.9974120
#> 430_Epb41l2             0.9974120
#> 431_Eif3b               1.0000000
#> 432_Isyna1              0.9974120
#> 433_U2surp              0.9974120
#> 434_LOC100765020        0.9974120
#> 435_Arhgef6             0.9974120
#> 436_Ptpn1               0.9974120
#> 437_Prpf4b              0.9974120
#> 438_Rpl35a              0.9974120
#> 439_Prpf4b              1.0000000
#> 440_Zyx                 0.9974120
#> 441_Dbn1                0.9974120
#> 442_Chaf1b              1.0000000
#> 443_LOC113834282        0.9974120
#> 444_Gpsm2               0.9974120
#> 445_LOC100757535        1.0000000
#> 446_Cfap410             0.9974120
#> 447_Epb41l2             0.9974120
#> 448_Ncbp1               1.0000000
#> 449_Pacsin1             1.0000000
#> 450_Cstf2               0.9974120
#> 451_LOC100769437        0.9974120
#> 452_eIF2aK2             0.9974120
#> 453_Kiaa1191            0.9974120
#> 454_Mepce               0.9974120
#> 455_Cbx8                0.9974120
#> 456_Eed                 0.9974120
#> 457_Cdc42ep1            0.9974120
#> 458_Lrrfip2             1.0000000
#> 459_Pacsin1             0.9974120
#> 460_Gpatch4             0.9974120
#> 461_Plin4               1.0000000
#> 462_NA                  0.9974120
#> 463_Snip1               0.9974120
#> 464_Cyld                0.9974120
#> 465_Plin4               0.9974120
#> 466_Twf1                1.0000000
#> 467_LOC113834282        1.0000000
#> 468_Snip1               0.9974120
#> 469_Ppp4r3a             0.9974120
#> 470_Psip1               0.9974120
#> 471_Dnajc5              1.0000000
#> 472_Phf8                0.9974120
#> 473_Bola1               0.9974120
#> 474_Cdc42ep1            0.9974120
#> 475_Eif4ebp2            0.9974120
#> 476_Prpf38b             0.9974120
#> 477_Klhl26              0.9974120
#> 478_Hsph1               1.0000000
#> 479_Snip1               0.9974120
#> 480_Caskin2             0.9974120
#> 481_Plpp6               0.9974120
#> 482_NA                  0.9974120
#> 483_Mlh1                0.9974120
#> 484_Gys1                0.9995772
#> 485_Tfg                 0.9974120
#> 486_Arhgef6             1.0000000
#> 487_Mphosph10           1.0000000
#> 488_Hoxc10              0.9974120
#> 489_LOC100759640        1.0000000
#> 490_Arhgef40            0.9974120
#> 491_Dnajc5              1.0000000
#> 492_Tbc1d23             1.0000000
#> 493_Ubxn1               0.9974120
#> 494_Rab1a               1.0000000
#> 495_Eif3b               1.0000000
#> 496_Tceal8              0.9974120
#> 497_Dlgap4              0.9974120
#> 498_Smim13              0.9974120
#> 499_NA                  1.0000000
#> 500_Lrch4               1.0000000
#> 501_Bola1               1.0000000
#> 502_NA                  0.9974120
#> 503_Ptpn14              0.9974120
#> 504_LOC100759640        0.9974120
#> 505_Rps10               0.9974120
#> 506_Top2b               1.0000000
#> 507_Ssr3                0.9974120
#> 508_Homer3              1.0000000
#> 509_Phf8                0.9974120
#> 510_LOC100767716        0.9974120
#> 511_Xpa                 0.9974120
#> 512_H671_21690          0.9974120
#> 513_LOC100769471        0.9974120
#> 514_Gas2l1              0.9974120
#> 515_Luzp1               0.9974120
#> 516_Gpbp1               0.9974120
#> 517_Gatad2b             1.0000000
#> 518_Gys1                0.9974120
#> 519_Top2b               0.9974120
#> 520_LOC100757535        0.9974120
#> 521_Lpcat4              0.9974120
#> 522_Arhgef6             1.0000000
#> 523_Cavin3              0.9974120
#> 524_Gpatch4             0.9974120
#> 525_Prpf38b             0.9974120
#> 526_Timm8a              0.9974120
#> 527_Cavin3              0.9974120
#> 528_Mkrn2               0.9974120
#> 529_Oser1               0.9974120
#> 530_Gsk3b               0.9974120
#> 531_Eef1b2              1.0000000
#> 532_Ampd2               0.9974120
#> 533_Lrrfip2             0.9974120
#> 534_Ring1               0.9974120
#> 535_Rlim                0.9974120
#> 536_LOC100759640        0.9974120
#> 537_LOC100759640        0.9974120
#> 538_Atp5pf              0.9974120
#> 539_Max                 0.9974120
#> 540_Bap1                1.0000000
#> 541_Nsfl1c              0.9974120
#> 542_Prpf4b              1.0000000
#> 543_LOC100757535        0.9974120
#> 544_Mtmr10              0.9974120
#> 545_Hoxc10              0.9974120
#> 546_Trim35              1.0000000
#> 547_Eif4ebp2            0.9974120
#> 548_Dlgap4              1.0000000
#> 549_Gys1                0.9974120
#> 550_Sgtb                0.9974120
#> 551_Eri2                0.9974120
#> 552_Ccnd3               1.0000000
#> 553_Smim13              0.9974120
#> 554_Snrk                0.9974120
#> 555_Caskin2             1.0000000
#> 556_Pdcd11              0.9974120
#> 557_Pgam5               0.9974120
#> 558_Mphosph10           1.0000000
#> 559_Mideas              1.0000000
#> 560_Top2b               0.9974120
#> 561_LOC100763014        0.9974120
#> 562_Snip1               1.0000000
#> 563_Ubxn1               0.9974120
#> 564_LOC100750407        0.9974120
#> 565_Morf4l2             1.0000000
#> 566_Ctdspl2             1.0000000
#> 567_Cwf19l1             1.0000000
#> 568_Eef1b2              1.0000000
#> 569_C1H12orf45          1.0000000
#> 570_Znf367              0.9974120
#> 571_Ankrd34a            0.9974120
#> 572_Mllt11              0.9974120
#> 573_LOC100774792        0.9974120
#> 574_NA                  0.9974120
#> 575_Cbx8                0.9974120
#> 576_Bckdk               0.9974120
#> 577_Snip1               1.0000000
#> 578_Nsfl1c              0.9974120
#> 579_Gas2l1              1.0000000
#> 580_Nudc                0.9974120
#> 581_Epb41l2             0.9974120
#> 582_Mtmr6               1.0000000
#> 583_Znf668              0.9974120
#> 584_Hsph1               0.9974120
#> 585_LOC113834282        0.9974120
#> 586_Ctdspl2             0.9974120
#> 587_Foxf1               0.9974120
#> 588_Luzp1               0.9974120
#> 589_Xpa                 1.0000000
#> 590_Psip1               0.9974120
#> 591_Rbm7                0.9974120
#> 592_Mtrex               0.9974120
#> 593_Arhgef40            0.9974120
#> 594_Plekho2             0.9974120
#> 595_Bckdk               0.9974120
#> 596_Dut                 1.0000000
#> 597_Abcf1               1.0000000
#> 598_Txnl1               1.0000000
#> 599_Nudc                1.0000000
#> 600_Sh3gl1              0.9974120
#> 601_Gatad2b             0.9974120
#> 602_Homer3              0.9974120
#> 603_Septin6             0.9974120
#> 604_Smim13              0.9974120
#> 605_Arhgef40            0.9974120
#> 606_Rpl32               1.0000000
#> 607_Tomm34              0.9974120
#> 608_Mlh1                0.9974120
#> 609_Tbcc                1.0000000
#> 610_Eif3d               0.9974120
#> 611_Snrk                0.9974120
#> 612_Bckdk               1.0000000
#> 613_Wdr3                0.9974120
#> 614_LOC100757535        1.0000000
#> 615_Dlg1                0.9974120
#> 616_LOC100767716        0.9974120
#> 617_Hnrnpc              0.9974120
#> 618_Mphosph10           0.9974120
#> 619_Eif3b               0.9974120
#> 620_Emd                 1.0000000
#> 621_Txlng               0.9974120
#> 622_Prpf4b              1.0000000
#> 623_Rlim                1.0000000
#> 624_Eef1b2              0.9974120
#> 625_Def6                0.9974120
#> 626_LOC100765020        0.9974120
#> 627_U2surp              0.9974120
#> 628_Elf2                0.9974120
#> 629_Slc1a5              0.9974120
#> 630_NA                  0.9974120
#> 631_Tfg                 1.0000000
#> 632_Top2b               0.9974120
#> 633_Pip4p2              0.9974120
#> 634_Cdc42ep1            1.0000000
#> 635_Hsph1               1.0000000
#> 636_Twf1                1.0000000
#> 637_Nbn                 0.9974120
#> 638_Psmd4               0.9974120
#> 639_Bap1                0.9974120
#> 640_Mepce               0.9974120
#> 641_Mideas              0.9974120
#> 642_LOC100759640        0.9974120
#> 643_Epb41l2             1.0000000
#> 644_Sav1                0.9974120
#> 645_Prpf4b              0.9974120
#> 646_Gnas                0.9974120
#> 647_Mllt1               1.0000000
#> 648_Poldip3             0.9974120
#> 649_Aldoa               0.9974120
#> 650_Rbbp8               0.9974120
#> 651_LOC113834282        0.9974120
#> 652_Gys1                0.9974120
#> 653_Hnrnpc              0.9974120
#> 654_Vps35               0.9974120
#> 655_Miga2               0.9974120
#> 656_Epb41l2             0.9974120
#> 657_Tob2                1.0000000
#> 658_Lamtor1             1.0000000
#> 659_LOC100759640        1.0000000
#> 660_Epb41l2             0.9974120
#> 661_Rlim                0.9974120
#> 662_Gys1                1.0000000
#> 663_LOC100750437        0.9974120
#> 664_NA                  0.9974120
#> 665_Nbn                 0.9974120
#> 666_Tyw3                1.0000000
#> 667_Gas2l1              0.9974120
#> 668_Fus                 1.0000000
#> 669_Prpf38b             0.9974120
#> 670_Calu                1.0000000
#> 671_Rras2               0.9974120
#> 672_Prpf4b              0.9974120
#> 673_Nelfa               0.9974120
#> 674_LOC100754077        0.9974120
#> 675_Rbm28               0.9974120
#> 676_Nsfl1c              0.9974120
#> 677_Rnf126              0.9974120
#> 678_Eme1                0.9974120
#> 679_Nbn                 0.9974120
#> 680_Eif4ebp2            0.9974120
#> 681_Wee1                1.0000000
#> 682_Prpf38b             1.0000000
#> 683_Luzp1               1.0000000
#> 684_Gas2l1              1.0000000
#> 685_Pdcd11              0.9974120
#> 686_Chaf1b              0.9974120
#> 687_Pycr1               0.9974120
#> 688_Phf8                0.9974120
#> 689_Raver1              1.0000000
#> 690_Dbn1                0.9974120
#> 691_Dut                 0.9974120
#> 692_Prpf4b              1.0000000
#> 693_Prpf4b              0.9974120
#> 694_Efs                 0.9974120
#> 695_NA                  0.9974120
#> 696_Ppp2r5b             1.0000000
#> 697_Caskin2             0.9974120
#> 698_Arhgef40            0.9974120
#> 699_Zyx                 1.0000000
#> 700_Mphosph10           0.9974120
#> 701_LOC113833392        0.9974120
#> 702_Cdc42ep1            0.9974120
#> 703_Snrpa1              1.0000000
#> 704_Ncbp1               0.9974120
#> 705_Gas2l1              1.0000000
#> 706_Gas2l1              0.9974120
#> 707_Bap1                0.9974120
#> 708_LOC100759640        0.9974120
#> 709_Cherp               1.0000000
#> 710_Nbn                 0.9974120
#> 711_LOC100759640        1.0000000
#> 712_NA                  1.0000000
#> 713_Eif3b               0.9974120
#> 714_Miga2               1.0000000
#> 715_Prpf4b              0.9974120
#> 716_Dbn1                0.9974120
#> 717_Ppp2r5b             0.9974120
#> 718_Exosc9              1.0000000
#> 719_Eif3b               0.9974120
#> 720_Ripk2               1.0000000
#> 721_Dlg1                0.9974120
#> 722_N4bp1               0.9974120
#> 723_Nudc                0.9974120
#> 724_Znf367              0.9974120
#> 725_Ring1               1.0000000
#> 726_Snrpa1              0.9974120
#> 727_U2surp              0.9974120
#> 728_LOC100764225        0.9974120
#> 729_Cdc42ep1            1.0000000
#> 730_Znf385a             1.0000000
#> 731_Ints1               1.0000000
#> 732_LOC113833392        1.0000000
#> 733_Lrch4               0.9974120
#> 734_Ctdspl2             0.9974120
#> 735_Prpf4b              0.9974120
#> 736_Luzp1               0.9974120
#> 737_Eif3b               1.0000000
#> 738_Ptpn14              1.0000000
#> 739_Rrp1                0.9974120
#> 740_Lrrfip2             1.0000000
#> 741_Nsfl1c              1.0000000
#> 742_Ddx51               0.9974120
#> 743_Prpf38b             0.9974120
#> 744_Eef1b2              0.9974120
#> 745_Znf385a             1.0000000
#> 746_Map9                1.0000000
#> 747_Rflnb               1.0000000
#> 748_NA                  0.9974120
#> 749_C1H12orf45          0.9974120
#> 750_U2surp              0.9974120
#> 751_Caskin2             0.9974120
#> 752_Eri1                0.9974120
#> 753_Gsk3b               0.9974120
#> 754_LOC100766946        0.9974120
#> 755_Cnpy3               0.9974120
#> 756_Hnrnpc              0.9974120
#> 757_Ptpn14              0.9974120
#> 758_Slc7a11             0.9974120
#> 759_Hnrnpc              0.9974120
#> 760_Cdc37l1             0.9974120
#> 761_LOC100768405        0.9974120
#> 762_Rragc               0.9974120
#> 763_LOC113834282        0.9974120
#> 764_Fus                 0.9974120
#> 765_Ubxn1               0.9974120
#> 766_Mmut                0.9974120
#> 767_Pdcd11              0.9974120
#> 768_LOC100757535        0.9974120
#> 769_Eif3b               0.9974120
#> 770_Rnf113a             1.0000000
#> 771_Sytl4               1.0000000
#> 772_Tlnrd1              0.9974120
#> 773_H671_1g1131         0.9974120
#> 774_Neurl1              0.9974120
#> 775_Zyx                 0.9974120
#> 776_Ctdspl2             0.9974120
#> 777_Chaf1b              0.9995772
#> 778_Rragc               1.0000000
#> 779_Srfbp1              0.9974120
#> 780_Gys1                0.9974120
#> 781_Usp15               0.9974120
#> 782_Arhgef40            0.9974120
#> 783_Gigyf1              0.9974120
#> 784_Minar1              0.9974120
#> 785_Dus2                1.0000000
#> 786_Gatad2b             0.9974120
#> 787_Eif5                0.9974120
#> 788_Epb41l2             0.9974120
#> 789_Arl6ip4             0.9974120
#> 790_Plin4               1.0000000
#> 791_Elf2                0.9974120
#> 792_Plin4               1.0000000
#> 793_Snip1               1.0000000
#> 794_Txlng               1.0000000
#> 795_LOC100769437        0.9974120
#> 796_Caskin2             0.9974120
#> 797_NA                  0.9974120
#> 798_Synm                1.0000000
#> 799_Synm                0.9974120
#> 800_Ube2c               0.9974120
#> 801_Sgtb                0.9974120
#> 802_Prpf4b              0.9974120
#> 803_Epb41l2             0.9974120
#> 804_Mllt1               0.9974120
#> 805_LOC100759640        1.0000000
#> 806_Epb41l2             0.9974120
#> 807_Znf280b             0.9974120
#> 808_Kiaa1143            1.0000000
#> 809_Gas2l1              1.0000000
#> 810_Srp72               0.9974120
#> 811_Tomm22              0.9974120
#> 812_Psip1               0.9974120
#> 813_Arhgef37            0.9974120
#> 814_Bckdk               0.9974120
#> 815_Strip1              0.9974120
#> 816_Usp15               1.0000000
#> 817_Ssr3                0.9974120
#> 818_Strip1              0.9974120
#> 819_Eif3b               0.9974120
#> 820_U2surp              0.9974120
#> 821_Bend3               0.9974120
#> 822_Rps10               0.9974120
#> 823_Rpl23a              0.9974120
#> 824_Nbn                 0.9974120
#> 825_Rpap3               0.9974120
#> 826_LOC100759640        1.0000000
#> 827_Ric8a               0.9974120
#> 828_Hsph1               0.9974120
#> 829_LOC100759640        0.9974120
#> 830_LOC100757535        0.9974120
#> 831_Gigyf1              0.9974120
#> 832_Dbn1                0.9974120
#> 833_Snrk                1.0000000
#> 834_Prpf38b             0.9974120
#> 835_LOC100766868        0.9974120
#> 836_LOC100766868        0.9974120
#> 837_Wbp11               0.9974120
#> 838_Rusc2               0.6510735
#> 839_Eif3b               0.9974120
#> 840_Ptpn14              1.0000000
#> 841_Rlim                0.9974120
#> 842_Ints1               1.0000000
#> 843_Chaf1b              0.9974120
#> 844_Dlg1                1.0000000
#> 845_Lamtor1             0.9974120
#> 846_Tab1                0.9974120
#> 847_Dbn1                0.9974120
#> 848_Psip1               0.9974120
#> 849_Dbn1                0.9974120
#> 850_Pabpc1              0.9974120
#> 851_Hnrnpc              0.9974120
#> 852_Emd                 0.9974120
#> 853_LOC100764225        0.9974120
#> 854_Nup50               1.0000000
#> 855_Ctcf                0.9974120
#> 856_Raly                0.9974120
#> 857_Bard1               0.9974120
#> 858_Ptpn14              0.9974120
#> 859_LOC100757535        0.9974120
#> 860_Psmd2               1.0000000
#> 861_Junb                0.9974120
#> 862_C1qbp               0.9974120
#> 863_Lrch4               0.9974120
#> 864_CUNH14orf93         1.0000000
#> 865_U2surp              0.9974120
#> 866_Raly                1.0000000
#> 867_LOC100774417        0.9974120
#> 868_Srp72               0.9974120
#> 869_LOC100764225        1.0000000
#> 870_Morf4l2             0.9974120
#> 871_CUNH9orf40          0.9974120
#> 872_Gas2l1              0.9974120
#> 873_Atp5pf              0.9974120
#> 874_Lrrfip2             0.9974120
#> 875_Prpf4b              0.9974120
#> 876_Top2b               1.0000000
#> 877_Mepce               0.9974120
#> 878_Ptpn14              0.9974120
#> 879_Dnajc25             0.9974120
#> 880_Cbx8                0.9974120
#> 881_Synm                0.9974120
#> 882_Def6                0.9974120
#> 883_Gys1                0.9974120
#> 884_Luzp1               0.9974120
#> 885_Synm                0.9974120
#> 886_Snip1               1.0000000
#> 887_Top2b               0.9974120
#> 888_NA                  0.9974120
#> 889_Trim35              0.9974120
#> 890_Znf385a             0.9974120
#> 891_Chaf1b              0.9974120
#> 892_Abcf1               0.9974120
#> 893_Pdcd11              1.0000000
#> 894_Dlg1                0.9974120
#> 895_Dbn1                1.0000000
#> 896_LOC100752363        0.9974120
#> 897_Ppp4r3a             1.0000000
#> 898_Gas2l1              0.9974120
#> 899_Mtmr10              1.0000000
#> 900_Cyld                0.9974120
#> 901_NA                  0.9974120
#> 902_Rnf113a             1.0000000
#> 903_Nelfa               1.0000000
#> 904_Zkscan1             1.0000000
#> 905_Chaf1b              0.9974120
#> 906_Eif3b               0.9974120
#> 907_Top2b               0.9995772
#> 908_Chaf1b              0.9974120
#> 909_Epb41l2             0.9974120
#> 910_C3H11orf58          0.9974120
#> 911_Top2b               0.9974120
#> 912_Wee1                0.9974120
#> 913_Raly                0.9974120
#> 914_H671_1g2680         0.9974120
#> 915_Eef1b2              1.0000000
#> 916_Gas2l1              1.0000000
#> 917_Epb41l2             0.9974120
#> 918_Rpl23a              0.9974120
#> 919_Chmp2b              0.9974120
#> 920_Lrrfip2             0.9974120
#> 921_Aldoa               0.9974120
#> 922_Cby1                0.9974120
#> 923_LOC100759640        0.9974120
#> 924_Rbm28               0.9974120
#> 925_Skiv2l              0.9974120
#> 926_Ints1               0.9974120
#> 927_Ehd1                0.9974120
#> 928_Nr2f6               0.9974120
#> 929_Top2b               0.9974120
#> 930_Lrrfip2             0.9974120
#> 931_Pip4p2              0.9974120
#> 932_Srp72               0.9974120
#> 933_Mtmr9               0.9974120
#> 934_Gigyf1              1.0000000
#> 935_Rbm7                0.9974120
#> 936_LOC100773565        0.9974120
#> 937_Trim35              0.9974120
#> 938_Cbx8                0.9974120
#> 939_Rplp0               0.9974120
#> 940_Aldoa               0.9974120
#> 941_NA                  0.9974120
#> 942_Zyx                 0.6510735
#> 943_Psip1               0.9974120
#> 944_Slc7a11             0.9974120
#> 945_Miga2               0.9974120
#> 946_Arhgef6             1.0000000
#> 947_Dlgap4              0.9974120
#> 948_Ampd2               1.0000000
#> 949_Luzp1               0.9974120
#> 950_Camlg               0.9974120
#> 951_Pfkfb3              1.0000000
#> 952_NA                  0.9974120
#> 953_Raly                0.9974120
#> 954_Kiaa1143            0.9974120
#> 955_Bcar1               0.9974120
#> 956_Gatad2b             0.9974120
#> 957_Eif4ebp2            0.9974120
#> 958_Fam76b              1.0000000
#> 959_Camlg               0.9974120
#> 960_LOC100754077        0.9974120
#> 961_NA                  0.9974120
#> 962_Epb41l2             0.9974120
#> 963_Ankrd34a            0.9995772
#> 964_Zc3h15              0.9974120
#> 965_Def6                0.9974120
#> 966_Srsf6               1.0000000
#> 967_H671_4g11480        0.9974120
#> 968_Top2b               0.9974120
#> 969_LOC100769471        0.9974120
#> 970_Raver1              0.9974120
#> 971_Etv3                0.9974120
#> 972_Psd                 0.9974120
#> 973_Usp15               1.0000000
#> 974_Nol7                0.9974120
#> 975_Stk38               0.9974120
#> 976_Smim13              0.9974120
#> 977_Etv3                0.9974120
#> 978_Synm                0.9974120
#> 979_Pwp1                0.9974120
#> 980_Fus                 0.9974120
#> 981_Junb                0.9974120
#> 982_Phf8                1.0000000
#> 983_Nelfa               0.9974120
#> 984_Prpf4b              0.9974120
#> 985_Abraxas1            0.9974120
#> 986_Prpf4b              0.9974120
#> 987_Raver1              0.9974120
#> 988_Caap1               0.9974120
#> 989_Rpap3               1.0000000
#> 990_Hsph1               0.9974120
#> 991_LOC100750437        1.0000000
#> 992_Mepce               0.9974120
#> 993_Efs                 0.9974120
#> 994_Epb41l2             0.9974120
#> 995_Abcf1               0.9974120
#> 996_NA                  0.9974120
#> 997_Eif4ebp2            0.9974120
#> 998_Pfkfb3              0.9974120
#> 999_Hnrnpc              0.9974120
#> 1000_Psmd2              0.9974120
#> 
#> $Stationary$pvc_pattern_summary
#>   -60 15 60 90 120 240
#> p   0  0 40  0   0   0
#> v   0  0 22  0   0   0
#> b   0  0  1  0   0   0
#> t   0  0  1  0   0   0
```

## Create the plots and HTML report

``` r
plot_info <- list(
    y_axis_label = "log2 intensity",
    time_unit = "min",
    treatment_labels = list(
        Exponential = "temp shift",
        Stationary = "temp shift"
    ),
    treatment_timepoints = list(
        Exponential = 146,
        Stationary = 146
    )
)

SplineOmics::create_pvc_report(
    splineomics = splineomics,
    pvc_results = pvc_results,
    plot_info = plot_info,
    report_dir = file.path(tempdir(), "splineomics_pvc_report")
) 
#> $Exponential
#> $Exponential$pvc_adj_pvals
#>                  15_vs_neighbors 60_vs_neighbors 90_vs_neighbors
#> 1_Top2b                0.9996591      0.89228318    0.9529197807
#> 2_NA                   1.0000000      1.00000000    0.6654892405
#> 3_Snip1                0.7976077      1.00000000    1.0000000000
#> 4_Tomm34               0.5975031      0.95301168    0.3458807846
#> 5_Pus3                 0.8241813      0.99844582    0.5103974918
#> 6_Ints1                0.7391776      1.00000000    1.0000000000
#> 7_Mlh1                 0.9996591      0.57340016    0.4231947052
#> 8_LOC100750437         0.9996591      0.91476163    0.4987772572
#> 9_Pabpc1               0.5731500      0.46668044    0.3669755656
#> 10_Top2b               0.9136090      0.89453090    0.4358371316
#> 11_Gorasp1             0.7710625      0.76519496    0.3357324904
#> 12_Ints1               0.7239466      0.57955504    0.5486871521
#> 13_Syvn1               0.9996591      0.63814186    0.8053523228
#> 14_Znf280b             0.7528198      0.98300959    0.9701805996
#> 15_Mrnip               0.8241813      0.57340016    0.8085069881
#> 16_Rragc               0.8915063      0.77220390    0.4872552464
#> 17_Gorasp1             0.9957487      0.77660746    0.3154653432
#> 18_Tomm34              0.9173402      0.92922311    0.6987798252
#> 19_LOC100757430        0.6485597      0.89453090    0.6803315820
#> 20_Ubxn1               0.9173402      0.78620962    0.9177365750
#> 21_H671_1g1131         0.6904936      0.40322706    1.0000000000
#> 22_Luzp1               0.9428566      0.94155079    0.8085069881
#> 23_Efs                 0.9996591      0.96609855    0.9271815649
#> 24_Mta2                0.9918588      0.96609855    0.5089919729
#> 25_Nedd1               0.8241813      0.61641970    0.1343299522
#> 26_Gigyf1              0.8539081      0.84244807    0.9562143923
#> 27_Myh9                0.8711111      0.64588748    0.6368831503
#> 28_Caskin2             0.5975031      0.89453090    0.3635388719
#> 29_Papolg              0.5975031      0.89453090    0.6353934056
#> 30_Tfg                 0.9996591      0.91274220    0.6992386870
#> 31_Rpl34               0.4335248      0.39726695    0.3277644801
#> 32_Mideas              0.9996591      0.74871175    0.1583039647
#> 33_Gys1                0.6481524      0.57340016    0.6460689701
#> 34_Arhgef6             0.9085372      0.89453090    0.2796179498
#> 35_Ctdspl2             0.9406578      0.96679510    0.9174607580
#> 36_Ptpn14              0.9206541      0.91266136    0.6803315820
#> 37_Raly                1.0000000      1.00000000    1.0000000000
#> 38_Znhit3              0.8775637      0.99844582    0.8179717078
#> 39_LOC113833392        1.0000000      1.00000000    1.0000000000
#> 40_Luc7l3              1.0000000      1.00000000    0.8930375367
#> 41_Rplp0               1.0000000      1.00000000    1.0000000000
#> 42_Gys1                0.9862244      0.96609855    0.2738360002
#> 43_Rpl22l1             0.9206541      0.65495110    0.2090589517
#> 44_Eif3b               0.8478867      0.89228318    0.5837544074
#> 45_Med26               0.8478867      0.97469270    0.2491244798
#> 46_Mepce               0.9996591      0.96707635    0.5299158386
#> 47_Pdcd11              0.8818738      0.91266136    0.7448404619
#> 48_Twf1                0.7976077      1.00000000    1.0000000000
#> 49_LOC100759640        1.0000000      1.00000000    1.0000000000
#> 50_Wrnip1              0.9937178      0.84503288    0.4987772572
#> 51_Poldip3             0.8340945      0.54140969    0.5103974918
#> 52_Ampd2               0.9409124      0.90654204    0.4184815759
#> 53_Mea1                0.9937178      0.60337996    0.3770190551
#> 54_Dbn1                0.3045423      0.32541914    0.3380338135
#> 55_Snip1               0.7140663      0.62588442    0.4391933166
#> 56_Srsf6               0.9307018      0.89453090    0.8591498257
#> 57_LOC113834282        0.7528198      0.97469270    0.1374908882
#> 58_Map9                0.8702312      0.89453090    0.9464121243
#> 59_Cdc42ep1            0.7528198      0.84244807    0.4084783215
#> 60_Poldip3             0.7619059      0.99844582    0.9187143328
#> 61_LOC100764225        1.0000000      1.00000000    1.0000000000
#> 62_Epb41l2             1.0000000      1.00000000    1.0000000000
#> 63_H671_4g11480        0.4884807      0.41581591    0.3669729209
#> 64_Nbn                 1.0000000      1.00000000    1.0000000000
#> 65_U2surp              1.0000000      1.00000000    1.0000000000
#> 66_Gigyf1              0.5975031      0.87597326    0.4240325835
#> 67_NA                  0.4335248      0.58893964    0.4007519987
#> 68_Luc7l3              1.0000000      0.76571006    0.4034772289
#> 69_LOC100752363        1.0000000      0.15911620    0.3352199052
#> 70_Ampd2               0.9996591      0.39726695    0.1163090667
#> 71_LOC100759640        0.9996591      0.39726695    0.1163090667
#> 72_Stam                0.9996591      0.84244807    0.7566755605
#> 73_Nsfl1c              0.9865255      0.98666147    0.9164066845
#> 74_Pfkfb3              0.9996591      0.89979241    0.9164066845
#> 75_Rad23a              0.9409124      0.89453090    0.4034772289
#> 76_Elf2                0.9996591      0.89453090    0.8085069881
#> 77_Crem                0.8540623      0.99844582    0.6437184324
#> 78_Rragc               0.9996591      0.84518166    0.6516309777
#> 79_Lrrfip2             0.9996591      0.96609855    0.9797810656
#> 80_Zyx                 0.7976077      0.39726695    0.1163090667
#> 81_Lrrfip2             1.0000000      1.00000000    1.0000000000
#> 82_Gatad2b             0.6473547      0.97469270    0.3587160453
#> 83_Bcar1               0.9426037      0.89148887    0.9993185527
#> 84_Ehd1                0.8608189      0.67222253    0.6468176079
#> 85_LOC113834282        0.8241813      0.84503288    0.4945708006
#> 86_Tmem230             0.4884807      0.75726680    0.4263010375
#> 87_Ncbp1               0.8711111      0.91266136    0.7814977065
#> 88_Mllt1               1.0000000      1.00000000    1.0000000000
#> 89_Stk17b              0.9085372      0.51929656    0.1163090667
#> 90_Dlgap4              0.4884807      0.75726680    0.4263010375
#> 91_Papolg              0.6933908      0.60337996    0.3186191377
#> 92_Cyld                0.5971522      0.58893964    0.3321813438
#> 93_Gigyf1              0.9937178      0.89453090    0.7448404619
#> 94_Lrrfip2             0.7027612      0.76571006    0.9798872443
#> 95_Lrrfip2             0.9996591      0.95179760    0.6629729761
#> 96_Rlim                0.7911358      0.92096884    0.9914761578
#> 97_Eif3b               0.8751329      0.94478306    0.9456073392
#> 98_Mphosph10           1.0000000      0.76511779    0.6044440309
#> 99_Gatad2b             0.8539081      0.57340016    0.4945708006
#> 100_Srsf6              0.5731500      0.40573021    0.1542353843
#> 101_Zyx                0.9865255      0.97469270    0.9271815649
#> 102_Mphosph10          0.7765806      0.39726695    0.1175948545
#> 103_Psip1              0.9996591      0.41581591    0.3669755656
#> 104_Fbl                1.0000000      1.00000000    1.0000000000
#> 105_H671_1g2680        0.9996591      0.76519496    0.3442993309
#> 106_Sgtb               0.8350348      0.97580669    0.6199397072
#> 107_Gnl3               0.7619059      0.99478987    0.8364333175
#> 108_Eif3b              0.7619059      0.65495110    0.8053523228
#> 109_Serpinb1           0.7528198      0.60337996    0.3444986722
#> 110_N4bp1              0.9409124      0.98261381    0.5103974918
#> 111_Snip1              0.8540623      0.98444277    0.6647202777
#> 112_Psip1              1.0000000      1.00000000    1.0000000000
#> 113_Mlh1               1.0000000      1.00000000    1.0000000000
#> 114_Bsg                0.9996591      0.39726695    0.2700018411
#> 115_Tnpo1              1.0000000      1.00000000    1.0000000000
#> 116_H671_1g2680        0.9530303      0.44163335    0.1163090667
#> 117_Cbx8               0.9362984      0.57340016    0.2646284031
#> 118_Mideas             0.9530303      0.44163335    0.1163090667
#> 119_Mideas             0.5975031      0.89453090    0.9628343032
#> 120_Dcun1d3            0.5975031      0.89453090    0.9628343032
#> 121_Dlg1               0.8478867      0.57340016    0.2221736132
#> 122_Rad23a             0.7708242      0.80301683    0.5641930786
#> 123_Srsf6              0.3045423      0.59567667    0.8287783448
#> 124_Stx7               0.9206541      0.89827722    0.8267719777
#> 125_Pdcd11             0.5601500      0.76519496    0.9973033256
#> 126_Kiaa1958           0.6657483      0.57955504    0.2543712380
#> 127_Pwp1               0.9377457      0.97864177    0.5103974918
#> 128_Txlng              0.5975031      0.68059529    0.6987798252
#> 129_Junb               0.8241813      0.92070000    0.7874312574
#> 130_LOC100759640       0.7224935      0.99844582    0.8318410441
#> 131_Dbn1               0.7391776      0.84244807    0.8474656873
#> 132_Top2b              0.9307018      0.72336026    0.1975058450
#> 133_Rusc2              0.8443625      0.89453090    0.8140711759
#> 134_NA                 1.0000000      1.00000000    1.0000000000
#> 135_LOC113837251       0.9996591      0.95561081    0.7126478180
#> 136_Fam76b             0.9173402      0.95301168    0.7566755605
#> 137_Ptpn14             0.9206541      0.99844582    0.4987772572
#> 138_Chmp4b             0.7391776      0.93310639    0.6987798252
#> 139_Prpf4b             0.5975031      0.94478306    0.6629729761
#> 140_Eif3b              0.4884807      1.00000000    1.0000000000
#> 141_Nsfl1c             0.8617286      0.84244807    0.8053523228
#> 142_Pdlim7             0.4884807      0.94478306    0.5103974918
#> 143_Rnf113a            0.7619059      0.58902875    0.4014644565
#> 144_Epb41l2            0.8447753      0.63814186    0.3587160453
#> 145_Hnrnpc             0.9409124      0.79747244    0.2480080014
#> 146_LOC113834282       0.5601500      0.97469270    0.8026677103
#> 147_Plekho2            0.9937178      0.76519496    0.2738360002
#> 148_Med26              1.0000000      0.84244807    0.8140711759
#> 149_Arhgef40           0.8327078      0.88113882    0.4945708006
#> 150_NA                 0.9937178      0.78620962    0.6368831503
#> 151_Phf8               0.8241813      0.86446114    0.9464121243
#> 152_Minar1             1.0000000      1.00000000    1.0000000000
#> 153_H671_21690         0.9361973      0.99844582    0.8473517850
#> 154_Arhgef40           0.9865255      0.83239038    0.3458807846
#> 155_Chaf1b             0.9805267      0.91266136    0.7100903489
#> 156_Prpf4b             0.9937178      0.89453090    0.9405569286
#> 157_Znf367             0.8447753      0.41581591    0.2738360002
#> 158_Luzp1              0.9996224      0.99844582    0.9164066845
#> 159_LOC113833882       0.7765806      0.84244807    0.8474656873
#> 160_Hnrnpc             0.4884807      0.59208712    0.6987798252
#> 161_Mepce              0.7976077      0.88610890    0.4888053927
#> 162_Ubxn1              0.8443625      0.87450350    0.7448404619
#> 163_Mllt1              0.9996591      0.81902350    0.3075246846
#> 164_Chaf1b             1.0000000      1.00000000    1.0000000000
#> 165_Raly               1.0000000      1.00000000    1.0000000000
#> 166_Gas2l1             0.7239466      0.41853713    0.2038598949
#> 167_Dlg1               0.9996591      0.58893964    0.1163090667
#> 168_Hoxc10             1.0000000      1.00000000    1.0000000000
#> 169_Gigyf1             1.0000000      1.00000000    1.0000000000
#> 170_Luzp1              0.9937178      0.89453090    0.8474656873
#> 171_Srp72              0.7976077      0.84244807    0.3989353266
#> 172_LOC100771461       0.7619059      0.57340016    0.4052192524
#> 173_Chaf1b             0.8818738      0.78620962    0.9271815649
#> 174_C3H11orf58         0.9996591      0.77083169    0.8053523228
#> 175_Pdcd11             0.8078319      0.65561346    0.7977575454
#> 176_Psip1              0.8540623      0.96609855    0.9848227291
#> 177_Prpf4b             0.9996591      0.89771461    0.3669755656
#> 178_Rnf113a            0.6473547      0.98840974    0.2729399625
#> 179_Irf3               1.0000000      0.64588748    0.5103974918
#> 180_Smim13             0.5682360      0.81938063    0.7074776148
#> 181_Gnl3               0.4884807      0.84244807    0.5981207601
#> 182_Psma5              0.8499250      0.67222253    0.2203807883
#> 183_Ptpn14             0.8818738      0.75726680    0.7809050456
#> 184_Prpf4b             0.4335248      0.39726695    0.5356628179
#> 185_Top2b              0.7528198      0.89228318    0.9271815649
#> 186_Prpf38b            0.4763454      0.91266136    0.7448404619
#> 187_Epb41l2            0.9426037      0.61874884    1.0000000000
#> 188_Eif3b              0.9996591      0.32541914    0.2562757161
#> 189_Hnrnpc             0.8327078      0.73152943    0.2562114509
#> 190_LOC100758278       0.4884807      0.89453090    0.5016274831
#> 191_Prpf4b             0.8241813      0.57955504    0.3183586939
#> 192_Caskin2            0.5731500      0.01144135    0.0001024892
#> 193_LOC100752363       0.9996591      0.84503288    0.6352867542
#> 194_Septin6            0.8540623      0.71646090    0.3154653432
#> 195_Max                0.8789016      0.57340016    0.1975058450
#> 196_Mid1ip1            0.7391776      0.94044689    0.6352867542
#> 197_NA                 0.8540623      0.94208942    0.9201386737
#> 198_Hsph1              1.0000000      1.00000000    0.7068411866
#> 199_Nol7               0.5682360      0.57340016    0.6353934056
#> 200_Raly               0.6473547      0.74503582    0.8094385315
#> 201_Smim13             0.7528198      0.90047417    0.5457472072
#> 202_LOC100757535       0.4884807      0.85102646    0.2611499299
#> 203_Net1               0.4884807      0.93743843    0.3937197943
#> 204_LOC100754077       1.0000000      1.00000000    1.0000000000
#> 205_Snip1              0.5601500      0.96609855    0.2116134333
#> 206_Hnrnpc             1.0000000      1.00000000    0.0237123566
#> 207_Ldlrap1            0.9085372      0.99844582    0.9906710175
#> 208_Luzp1              0.8540623      0.81938063    0.3154653432
#> 209_Rpl26              0.7987719      0.57340016    0.5243063397
#> 210_Epb41l2            0.8179315      0.89771461    0.2738360002
#> 211_Znf367             0.7239466      0.92434829    0.2992400844
#> 212_Dlgap4             0.6213394      0.97469270    0.3390104519
#> 213_Plekho2            1.0000000      1.00000000    1.0000000000
#> 214_Zpr1               0.8443625      0.86446114    0.3169935031
#> 215_Dlgap4             0.9996591      1.00000000    1.0000000000
#> 216_Def6               0.9865255      0.87450350    0.3243591041
#> 217_Eif4ebp2           0.9996224      0.63814186    0.3498488866
#> 218_Eef1b2             0.6485597      0.72336026    0.2728600300
#> 219_Rad23a             0.9362126      0.63814186    0.3974824450
#> 220_Morf4l2            0.6317540      0.63814186    0.8873517778
#> 221_Arhgef40           0.7911358      0.85102646    0.9973033256
#> 222_NA                 0.8443625      0.78296882    0.8179717078
#> 223_LOC100773565       0.5731500      0.62629746    0.4945708006
#> 224_Dus2               0.7585314      0.57340016    0.2682605900
#> 225_Pip4p2             0.9996591      0.89228318    0.3183586939
#> 226_Top2b              0.4937630      0.32541914    0.1492654263
#> 227_Znf280b            0.7365385      0.89453090    0.9690539503
#> 228_Pdcd11             0.5975031      0.62629746    0.7229251756
#> 229_Bckdk              0.7528198      0.81778149    0.2895665722
#> 230_Arhgef40           1.0000000      1.00000000    1.0000000000
#> 231_Mepce              0.5682360      0.32541914    0.1163090667
#> 232_Ccnd3              0.6473547      0.95648901    0.6803315820
#> 233_Phf8               0.7138456      0.75201157    0.9271815649
#> 234_H671_1g2680        0.6971286      1.00000000    1.0000000000
#> 235_Ell                0.9996224      0.70094781    0.4034772289
#> 236_U2surp             0.8540623      0.99844582    0.2941869072
#> 237_Rps10              0.8327078      0.39726695    0.3974824450
#> 238_Ctdspl2            0.7391776      0.50159122    0.1492654263
#> 239_Top2b              1.0000000      1.00000000    0.4028573087
#> 240_Msantd3            0.9551112      0.66173205    0.5329495286
#> 241_Fam76b             1.0000000      1.00000000    0.6498365687
#> 242_Ppp4r3a            0.4884807      0.77660746    0.4014644565
#> 243_Gpatch4            1.0000000      1.00000000    1.0000000000
#> 244_Nudc               0.9996591      0.66227289    0.7292065106
#> 245_Nol7               1.0000000      1.00000000    0.4987772572
#> 246_Plekho2            1.0000000      1.00000000    1.0000000000
#> 247_Prpf4b             0.9862244      0.52311323    0.1975058450
#> 248_Mta2               1.0000000      1.00000000    1.0000000000
#> 249_U2surp             0.7976077      0.89453090    0.8262137970
#> 250_Ubxn1              0.7595219      0.76934373    0.5103974918
#> 251_Rlim               0.8540623      0.89309712    0.8474656873
#> 252_Atat1              0.4335248      0.63814186    0.4987772572
#> 253_Ubxn1              0.7797809      0.97469270    0.8140711759
#> 254_H671_1g2680        0.8241813      0.97469270    0.7032969196
#> 255_eIF2aK2            0.4937630      0.57340016    0.2758043564
#> 256_Skiv2l             0.5975031      0.57340016    0.3480771364
#> 257_Rpl28              0.8447753      0.32541914    0.1896719181
#> 258_LOC100759640       0.8818738      0.39726695    0.1163090667
#> 259_Gatad2b            0.5731500      0.70094781    0.4028573087
#> 260_NA                 0.9996591      0.97469270    0.5302084872
#> 261_Gprasp1            0.8711111      0.88113882    0.3937197943
#> 262_Luzp1              0.4884807      0.90654204    0.3498488866
#> 263_Slc1a5             0.9996591      0.89453090    0.6416169647
#> 264_LOC113834282       0.8540623      0.39726695    0.3757441469
#> 265_Srsf6              0.8179315      0.77660746    0.9562143923
#> 266_Cdc42ep1           0.7976077      0.59650136    0.2992400844
#> 267_Net1               1.0000000      1.00000000    1.0000000000
#> 268_Caskin2            0.9996591      0.84518166    0.5103974918
#> 269_LOC100759640       0.5731500      0.39726695    0.4590470303
#> 270_Mideas             0.9173402      0.99844582    0.6803315820
#> 271_Luzp1              0.5731500      1.00000000    1.0000000000
#> 272_Emd                0.9996591      0.97469270    0.9291830695
#> 273_Plpp6              0.4937630      0.57340016    0.6180316715
#> 274_LOC100759640       0.7528198      0.89453090    0.2720197870
#> 275_Rps7               0.8478867      0.96609855    0.9178825492
#> 276_Fkbp1a             0.5842771      0.99844582    0.4638862180
#> 277_Gatad2b            0.9409124      0.57340016    0.1706411676
#> 278_Znf385a            0.9937178      0.73042151    0.2728600300
#> 279_Arhgef6            0.9173402      0.58893964    0.3458807846
#> 280_Slirp              0.9996591      0.57340016    0.3183586939
#> 281_Skiv2l             0.7528008      0.62588442    0.3390104519
#> 282_H671_21690         0.9937178      0.91266136    0.9273415365
#> 283_Kat8               1.0000000      1.00000000    0.9174607580
#> 284_Nkap               0.9173402      0.89148887    0.4363021752
#> 285_Gsk3b              0.7987719      0.94004182    0.8085069881
#> 286_Ints1              0.9307018      0.74691171    0.4007519987
#> 287_Gas2l1             0.6968896      0.39726695    1.0000000000
#> 288_LOC100759640       0.9996591      0.78620962    0.3095333237
#> 289_Top2b              0.9996591      0.90047417    0.6353934056
#> 290_Kif20b             0.9996591      0.77660746    0.8232359030
#> 291_Phf8               0.9996591      1.00000000    1.0000000000
#> 292_Snip1              0.5731500      0.94208942    0.2738360002
#> 293_Gsk3b              0.9996591      0.88700774    0.7448404619
#> 294_Caskin2            0.9996591      0.88700774    0.7448404619
#> 295_C3H11orf58         0.4884807      1.00000000    1.0000000000
#> 296_Lrch4              1.0000000      1.00000000    1.0000000000
#> 297_LOC113834282       0.9937178      0.85102646    0.8085069881
#> 298_LOC100750407       0.9967931      0.87647129    0.5329495286
#> 299_LOC113833392       0.9996591      0.57340016    0.2738360002
#> 300_LOC113833882       0.8540623      0.89453090    0.2728600300
#> 301_Ldlrap1            0.9173402      0.91476163    0.9798872443
#> 302_Wee1               0.9996224      0.70094781    0.3437425283
#> 303_Caap1              1.0000000      1.00000000    1.0000000000
#> 304_Eif4ebp2           0.9546909      0.98666147    0.6946987247
#> 305_Ripk2              0.9206541      0.99844582    0.6897028754
#> 306_Srp72              1.0000000      1.00000000    0.1691476201
#> 307_Taok2              0.9634942      1.00000000    1.0000000000
#> 308_Nr2f6              0.9206541      0.65495110    0.4987772572
#> 309_Arhgef40           0.9428566      0.66022793    0.3770190551
#> 310_Gys1               0.7987719      0.84244807    0.9733407167
#> 311_Dlg1               0.8241813      0.85826672    0.9909934070
#> 312_Vapb               0.6904936      0.63814186    0.4262173942
#> 313_LOC100757535       1.0000000      1.00000000    1.0000000000
#> 314_Mkrn2              0.9926304      0.85430905    0.8053523228
#> 315_Eif3b              0.5601500      0.39726695    0.4084783215
#> 316_Isyna1             0.9996591      0.54140969    0.1896719181
#> 317_Prpf4b             0.7461617      0.57340016    0.3624569089
#> 318_LOC113833882       0.4884807      0.96609855    0.8351075753
#> 319_Lrch4              0.4884807      0.39726695    0.8053523228
#> 320_Dbn1               0.5975031      0.89228318    0.9470511571
#> 321_Abcf1              0.9362893      0.57340016    0.2116134333
#> 322_Ints1              0.4884807      0.39726695    0.4518029413
#> 323_C3H11orf58         0.8987058      0.58893964    0.3390788295
#> 324_Psma5              0.5337685      0.50159122    1.0000000000
#> 325_Fundc1             0.5975031      0.42185492    0.6724947438
#> 326_Papolg             0.8192612      0.99844582    0.9562143923
#> 327_Mideas             0.9996591      0.98873241    0.6990690966
#> 328_Ubxn1              0.6748108      0.74404228    0.2988269514
#> 329_Synm               0.9996591      0.85102646    0.6417500446
#> 330_Arhgef6            1.0000000      1.00000000    1.0000000000
#> 331_Ptpn14             1.0000000      1.00000000    1.0000000000
#> 332_Pgrmc1             1.0000000      1.00000000    1.0000000000
#> 333_Myh9               0.9996224      0.96826768    0.8085069881
#> 334_Etv3               1.0000000      0.74796187    0.4358371316
#> 335_Ip6k1              1.0000000      1.00000000    0.4919700587
#> 336_Luzp1              0.9996591      0.93838860    0.9771317229
#> 337_Ptpn14             0.8241813      0.76519496    0.8474656873
#> 338_Caskin2            0.8353686      0.77660746    0.9628343032
#> 339_Chaf1b             1.0000000      1.00000000    1.0000000000
#> 340_Ubxn1              1.0000000      1.00000000    1.0000000000
#> 341_Ube2c              0.6485597      0.89453090    0.5329495286
#> 342_Gins2              0.4884807      0.41581591    0.4945708006
#> 343_Nlgn2              0.5601500      0.74796187    0.6182392042
#> 344_Nf2                0.9996591      0.84244807    0.7100903489
#> 345_Pip4p2             0.6971286      0.65688382    0.7566755605
#> 346_Emd                0.4884807      0.57340016    0.9562143923
#> 347_Top2b              0.9361973      0.84244807    0.7774012636
#> 348_Trim35             0.9173402      0.99844582    0.9798872443
#> 349_NA                 0.8818738      0.97469270    0.8140711759
#> 350_NA                 0.9996591      1.00000000    1.0000000000
#> 351_Mideas             0.9996591      0.61840397    0.3093690624
#> 352_Gas2l1             0.5731500      0.90047417    0.4240325835
#> 353_Ampd2              0.8540623      0.77660746    0.9177365750
#> 354_Calu               1.0000000      0.98168118    0.8053523228
#> 355_Fam76b             1.0000000      0.99844582    0.8364333175
#> 356_Dlg1               0.9996591      0.91476163    0.9562143923
#> 357_Srsf6              0.9551112      0.70094781    0.3458807846
#> 358_Chaf1b             0.9551112      0.70094781    0.3458807846
#> 359_Dbn1               0.9865255      0.97469270    0.8474656873
#> 360_Tcf25              0.8617286      0.70125132    0.3154653432
#> 361_Psip1              0.9996591      0.89453090    0.4459274268
#> 362_Cnpy3              0.7987719      0.73046398    0.1896719181
#> 363_LOC100759640       0.9409124      0.57340016    0.3587160453
#> 364_Zyx                0.8818738      0.90581890    1.0000000000
#> 365_Lrch4              0.5975031      0.63315255    0.7781328172
#> 366_Bola1              1.0000000      1.00000000    1.0000000000
#> 367_Znf385a            0.9996591      0.64588748    0.4945708006
#> 368_Kif20b             0.8443625      0.57340016    0.2895665722
#> 369_Ell                0.8617286      1.00000000    1.0000000000
#> 370_Ell                0.9307018      0.93310639    1.0000000000
#> 371_Srsf6              0.9996591      0.95236696    0.4116522570
#> 372_Pwp1               0.9996591      0.61840397    0.4034772289
#> 373_Def6               0.9556619      0.70094781    0.5512163402
#> 374_Cbx8               0.9996591      1.00000000    1.0000000000
#> 375_Ddx51              0.9362984      0.94208942    0.5274679042
#> 376_Psip1              0.9937178      0.49257001    0.5689471290
#> 377_Arhgef40           0.9937178      0.39726695    0.1175948545
#> 378_Raly               0.5556805      0.63814186    0.3587160453
#> 379_NA                 0.1475220      0.15274857    0.2090589517
#> 380_Lrrfip2            0.4884807      0.32541914    0.1682782731
#> 381_Gnl3               0.9996591      0.85430905    0.9562143923
#> 382_Caskin2            0.5731500      0.39726695    0.6897028754
#> 383_Rragc              0.7619059      0.63814186    0.4945708006
#> 384_Caskin2            0.1475220      0.15911620    0.1650552773
#> 385_Bcar1              0.9996591      0.77220390    0.7589135521
#> 386_Homer3             0.9996591      0.84244807    0.8664983034
#> 387_Luzp1              0.9996591      0.95561081    0.9464121243
#> 388_N4bp1              1.0000000      1.00000000    1.0000000000
#> 389_Ppp4r3a            0.9996591      0.97692950    0.7613779935
#> 390_H671_1g2680        1.0000000      1.00000000    1.0000000000
#> 391_Gnl3               0.7710625      0.62588442    0.1956873478
#> 392_Top2b              0.9996591      0.99844582    0.9562143923
#> 393_Oser1              1.0000000      0.85102646    0.3907241794
#> 394_Snrk               0.9996591      0.99844582    0.4987772572
#> 395_Kat8               1.0000000      1.00000000    1.0000000000
#> 396_Raver1             1.0000000      1.00000000    1.0000000000
#> 397_Pdcd11             0.8447753      0.90047417    0.8268537217
#> 398_Rps20              0.9996591      0.58893964    0.2738360002
#> 399_Bsg                1.0000000      1.00000000    1.0000000000
#> 400_Raly               0.9996591      0.89453090    0.2796179498
#> 401_Pdcd2              0.8447753      0.91560726    0.1706411676
#> 402_Caskin2            0.8540623      0.91476163    0.5302084872
#> 403_LOC100773571       0.9362893      0.76934373    0.2441304415
#> 404_Papolg             0.7528198      0.89309712    0.7046018595
#> 405_LOC100757535       0.8540623      0.83239038    0.9909934070
#> 406_Caap1              0.9996591      0.57340016    0.1163090667
#> 407_Psip1              0.8540623      0.76701587    0.9271815649
#> 408_Dbn1               0.9996591      0.42243195    0.1650552773
#> 409_Mta2               0.7671830      0.84244807    0.5902417364
#> 410_Abcf1              0.9937178      0.89148887    1.0000000000
#> 411_LOC100754108       0.9206541      0.98168118    0.9185882056
#> 412_Slirp              0.9996591      0.89228318    0.4111144282
#> 413_Nelfa              0.9173402      0.98168118    0.9426420690
#> 414_Aggf1              0.9996591      0.62588442    0.4945708006
#> 415_Bap1               0.9173402      0.98261381    0.7448404619
#> 416_Luc7l3             0.6213394      0.57340016    0.3669729209
#> 417_Rrp1               0.5731500      0.57955504    0.8873517778
#> 418_Wrnip1             0.9926304      1.00000000    1.0000000000
#> 419_NA                 1.0000000      1.00000000    0.9973033256
#> 420_Abcf1              0.9996591      0.95871390    0.6987798252
#> 421_Cluap1             0.3045423      0.89309712    0.8474656873
#> 422_Hnrnpc             0.9937178      0.84503288    0.2738360002
#> 423_Ptpn1              0.9504850      0.81115498    0.1706411676
#> 424_Myh9               0.4335248      0.50159122    0.5103974918
#> 425_U2surp             0.7359970      0.63814186    0.5103974918
#> 426_NA                 0.9996591      0.96609855    0.9291830695
#> 427_Arhgef40           1.0000000      1.00000000    0.7448404619
#> 428_Chaf1b             0.9206541      0.80342508    0.3989353266
#> 429_Prpf4b             0.9996591      0.73710357    0.6838225235
#> 430_Epb41l2            1.0000000      1.00000000    1.0000000000
#> 431_Eif3b              0.9862244      0.76511779    0.2441304415
#> 432_Isyna1             0.9996591      0.70125132    0.1542353843
#> 433_U2surp             0.7528198      0.57340016    0.4187510969
#> 434_LOC100765020       0.9996591      0.84244807    0.3676390944
#> 435_Arhgef6            0.9996224      0.84503288    0.3444986722
#> 436_Ptpn1              0.9173402      0.58038856    0.9733407167
#> 437_Prpf4b             0.7619059      0.75201157    0.7347336557
#> 438_Rpl35a             0.5731500      0.65688382    0.8085069881
#> 439_Prpf4b             1.0000000      1.00000000    1.0000000000
#> 440_Zyx                1.0000000      1.00000000    0.7046018595
#> 441_Dbn1               0.4884807      0.76672626    0.4945708006
#> 442_Chaf1b             1.0000000      0.47674732    0.5902417364
#> 443_LOC113834282       0.7976077      0.60337996    0.4927937173
#> 444_Gpsm2              0.9428566      0.90374179    0.7384669751
#> 445_LOC100757535       0.5731500      0.84503288    0.6146985018
#> 446_Cfap410            0.7619059      0.93310639    0.7977575454
#> 447_Epb41l2            0.5975031      0.57340016    0.2125739710
#> 448_Ncbp1              0.9996591      0.72336026    0.3444986722
#> 449_Pacsin1            0.9307018      0.65440610    0.9623469453
#> 450_Cstf2              0.9173402      0.90964926    0.8094385315
#> 451_LOC100769437       0.9996591      0.89453090    0.4391933166
#> 452_eIF2aK2            0.7911358      0.75726680    0.7589135521
#> 453_Kiaa1191           0.9996591      0.76519496    0.7977575454
#> 454_Mepce              0.9361973      0.93838860    0.8094385315
#> 455_Cbx8               0.4884807      0.57340016    0.8474656873
#> 456_Eed                0.7619059      0.61924693    0.1583039647
#> 457_Cdc42ep1           0.3045423      0.76519496    0.3474399848
#> 458_Lrrfip2            0.7619059      0.83239038    0.5103974918
#> 459_Pacsin1            0.8327078      0.83239038    0.6368831503
#> 460_Gpatch4            0.8818738      0.61840397    0.4111144282
#> 461_Plin4              0.8241813      0.74871175    0.4889794954
#> 462_NA                 0.7710625      0.44694320    0.9798872443
#> 463_Snip1              0.9307018      0.84436118    0.2738360002
#> 464_Cyld               0.8540623      0.73042151    0.6654892405
#> 465_Plin4              0.9996591      0.74871175    0.3557894948
#> 466_Twf1               1.0000000      1.00000000    1.0000000000
#> 467_LOC113834282       1.0000000      1.00000000    1.0000000000
#> 468_Snip1              0.9996591      0.88113882    0.9253735310
#> 469_Ppp4r3a            1.0000000      1.00000000    1.0000000000
#> 470_Psip1              0.8702312      0.96609855    0.5423506258
#> 471_Dnajc5             0.8478867      0.78620962    0.4945708006
#> 472_Phf8               0.9996591      0.41581591    0.1163090667
#> 473_Bola1              0.9085372      0.93310639    0.9973033256
#> 474_Cdc42ep1           0.4884807      0.32541914    0.3390104519
#> 475_Eif4ebp2           0.6160334      0.91266136    0.7083822981
#> 476_Prpf38b            0.7976077      0.63814186    0.6187384142
#> 477_Klhl26             0.9307018      0.89453090    0.3908588040
#> 478_Hsph1              1.0000000      1.00000000    1.0000000000
#> 479_Snip1              0.5654106      0.95648901    0.3474399848
#> 480_Caskin2            0.9996591      0.69844769    0.3458807846
#> 481_Plpp6              0.8915063      0.39726695    0.2738360002
#> 482_NA                 0.8478867      0.39726695    1.0000000000
#> 483_Mlh1               0.4884807      0.62588442    0.8053523228
#> 484_Gys1               0.4884807      0.76519496    0.9973033256
#> 485_Tfg                0.6473547      0.99844582    0.4034772289
#> 486_Arhgef6            1.0000000      1.00000000    0.3277644801
#> 487_Mphosph10          1.0000000      1.00000000    0.3277644801
#> 488_Hoxc10             0.6267469      0.32541914    0.1272404753
#> 489_LOC100759640       0.9771509      0.96826768    0.9914761578
#> 490_Arhgef40           1.0000000      1.00000000    1.0000000000
#> 491_Dnajc5             0.8711111      0.47176584    0.4919700587
#> 492_Tbc1d23            0.9996591      1.00000000    1.0000000000
#> 493_Ubxn1              0.4884807      0.81301115    0.2738360002
#> 494_Rab1a              1.0000000      1.00000000    1.0000000000
#> 495_Eif3b              0.8192612      0.89148887    0.7046018595
#> 496_Tceal8             0.9996591      0.99844582    0.7942831879
#> 497_Dlgap4             0.6968896      0.99844582    0.3770190551
#> 498_Smim13             0.9996591      0.89148887    0.2649667117
#> 499_NA                 0.6968896      0.57340016    0.4240325835
#> 500_Lrch4              0.9937178      0.86446114    0.1492654263
#> 501_Bola1              0.7619059      0.95179760    0.2895665722
#> 502_NA                 0.9926304      0.44694320    0.1343299522
#> 503_Ptpn14             0.9996591      0.96707922    0.7448404619
#> 504_LOC100759640       0.4884807      0.84138511    0.3779801995
#> 505_Rps10              0.5601500      0.61840397    0.8611757851
#> 506_Top2b              0.8156605      0.87450350    0.4159943815
#> 507_Ssr3               1.0000000      0.99844582    0.5329495286
#> 508_Homer3             0.8327078      0.98261381    0.6803315820
#> 509_Phf8               0.7976077      0.74796187    0.3937197943
#> 510_LOC100767716       0.7976077      0.78620962    0.2225177202
#> 511_Xpa                0.5731500      0.88113882    0.5696238511
#> 512_H671_21690         0.9996591      0.94208942    0.4945708006
#> 513_LOC100769471       0.8128341      0.77660746    0.3154653432
#> 514_Gas2l1             0.4884807      0.57340016    0.6353934056
#> 515_Luzp1              0.3045423      0.50709661    0.2090589517
#> 516_Gpbp1              0.5731500      0.76519496    0.9623341869
#> 517_Gatad2b            1.0000000      1.00000000    1.0000000000
#> 518_Gys1               0.8632347      0.95016268    0.5025865483
#> 519_Top2b              0.9996591      0.97469270    0.8262137970
#> 520_LOC100757535       0.5971522      1.00000000    1.0000000000
#> 521_Lpcat4             0.5975031      0.75726680    0.8967754084
#> 522_Arhgef6            0.4251762      0.87647129    0.8611757851
#> 523_Cavin3             0.9996591      0.97411817    0.9642694022
#> 524_Gpatch4            0.9996591      0.71801079    0.7514243941
#> 525_Prpf38b            0.4884807      0.63814186    0.8591498257
#> 526_Timm8a             0.4116137      0.49183911    0.2738360002
#> 527_Cavin3             0.8241813      0.88113882    0.5103974918
#> 528_Mkrn2              0.8039534      0.75726680    0.3458807846
#> 529_Oser1              0.5601500      0.57340016    0.3587160453
#> 530_Gsk3b              0.4116137      0.49183911    0.2738360002
#> 531_Eef1b2             0.6933908      1.00000000    1.0000000000
#> 532_Ampd2              0.8241813      0.99652195    0.7361793136
#> 533_Lrrfip2            0.9085372      0.74796187    0.6324435208
#> 534_Ring1              0.4884807      0.89453090    0.6654892405
#> 535_Rlim               1.0000000      1.00000000    1.0000000000
#> 536_LOC100759640       0.9409124      0.99844582    0.5512163402
#> 537_LOC100759640       0.7619059      0.84764375    0.9630822526
#> 538_Atp5pf             0.9146177      0.96609855    0.6314147294
#> 539_Max                0.6499875      0.76519496    0.8094385315
#> 540_Bap1               0.9337346      0.93190716    0.9234152431
#> 541_Nsfl1c             0.4884807      0.57955504    0.1163090667
#> 542_Prpf4b             1.0000000      1.00000000    0.9354251433
#> 543_LOC100757535       0.7911358      0.84244807    0.9273415365
#> 544_Mtmr10             0.9937178      0.84764375    0.4028573087
#> 545_Hoxc10             1.0000000      1.00000000    1.0000000000
#> 546_Trim35             1.0000000      1.00000000    1.0000000000
#> 547_Eif4ebp2           0.7987719      0.44808277    0.2941869072
#> 548_Dlgap4             0.4884807      0.57340016    1.0000000000
#> 549_Gys1               0.5975031      0.54140969    0.4084783215
#> 550_Sgtb               0.9996591      0.89309712    0.5302084872
#> 551_Eri2               0.7391776      0.89309712    0.9460896228
#> 552_Ccnd3              0.9996591      1.00000000    1.0000000000
#> 553_Smim13             0.6485597      0.89453090    0.3390104519
#> 554_Snrk               0.6738408      0.90442509    0.2738360002
#> 555_Caskin2            0.4884807      0.39726695    0.6266024511
#> 556_Pdcd11             1.0000000      1.00000000    1.0000000000
#> 557_Pgam5              0.9996591      0.96609855    0.8473517850
#> 558_Mphosph10          0.9085372      0.88246087    0.8474656873
#> 559_Mideas             0.9996591      0.99844582    0.7865438845
#> 560_Top2b              0.9996591      0.84503288    0.9291830695
#> 561_LOC100763014       0.7057658      0.47674732    0.3169935031
#> 562_Snip1              0.9173402      0.84244807    0.5778503376
#> 563_Ubxn1              1.0000000      0.78620962    0.5457472072
#> 564_LOC100750407       1.0000000      1.00000000    1.0000000000
#> 565_Morf4l2            1.0000000      1.00000000    1.0000000000
#> 566_Ctdspl2            0.9996591      0.89453090    0.9291830695
#> 567_Cwf19l1            0.5193881      0.50159122    0.7764469293
#> 568_Eef1b2             0.8751329      1.00000000    1.0000000000
#> 569_C1H12orf45         0.8751329      1.00000000    1.0000000000
#> 570_Znf367             1.0000000      1.00000000    1.0000000000
#> 571_Ankrd34a           0.7886993      0.57340016    0.5319516213
#> 572_Mllt11             0.7391776      0.97692950    0.9701805996
#> 573_LOC100774792       0.9996591      0.98300959    0.8179717078
#> 574_NA                 0.5975031      0.99652195    0.3610206615
#> 575_Cbx8               0.7461617      0.97469270    0.6182392042
#> 576_Bckdk              0.7391776      0.74703751    0.3365248466
#> 577_Snip1              1.0000000      1.00000000    1.0000000000
#> 578_Nsfl1c             1.0000000      1.00000000    1.0000000000
#> 579_Gas2l1             1.0000000      0.61641970    0.6928203220
#> 580_Nudc               0.7528198      0.63315255    0.7604572285
#> 581_Epb41l2            0.9307018      1.00000000    1.0000000000
#> 582_Mtmr6              0.6726106      0.69493031    0.8094385315
#> 583_Znf668             0.7619059      0.32541914    0.2116134333
#> 584_Hsph1              1.0000000      1.00000000    0.5329495286
#> 585_LOC113834282       0.8447753      0.97469270    0.9271815649
#> 586_Ctdspl2            1.0000000      1.00000000    1.0000000000
#> 587_Foxf1              0.8711111      0.69097390    0.3183586939
#> 588_Luzp1              0.8539081      0.98226136    0.7083822981
#> 589_Xpa                0.9937178      0.93770530    0.9701805996
#> 590_Psip1              0.8241813      0.76934373    0.3974824450
#> 591_Rbm7               0.4884807      0.50159122    0.8473517850
#> 592_Mtrex              0.8705424      0.44694320    0.9174607580
#> 593_Arhgef40           1.0000000      1.00000000    1.0000000000
#> 594_Plekho2            0.8818738      0.88113882    0.9623341869
#> 595_Bckdk              0.5975031      0.78620962    0.6629729761
#> 596_Dut                0.4884807      0.57340016    0.9104099835
#> 597_Abcf1              0.4884807      0.57340016    0.9104099835
#> 598_Txnl1              0.9996591      0.94478306    0.8085069881
#> 599_Nudc               0.8818738      0.94517210    0.8311625987
#> 600_Sh3gl1             0.7994190      0.39726695    0.1709825636
#> 601_Gatad2b            0.9996591      0.89309712    0.9973033256
#> 602_Homer3             1.0000000      1.00000000    0.6182392042
#> 603_Septin6            0.9996591      0.93838860    0.8271678824
#> 604_Smim13             0.6971286      0.34703515    0.3175586501
#> 605_Arhgef40           0.9996591      0.99480848    0.5457472072
#> 606_Rpl32              0.9996591      1.00000000    1.0000000000
#> 607_Tomm34             0.6968896      0.46668044    0.5329495286
#> 608_Mlh1               0.7619059      0.68895932    0.2555563503
#> 609_Tbcc               1.0000000      1.00000000    1.0000000000
#> 610_Eif3d              1.0000000      1.00000000    0.7589135521
#> 611_Snrk               1.0000000      0.96609855    0.6803315820
#> 612_Bckdk              1.0000000      1.00000000    0.9456073392
#> 613_Wdr3               1.0000000      1.00000000    1.0000000000
#> 614_LOC100757535       0.9634942      0.81115498    0.8794786048
#> 615_Dlg1               0.4116137      0.64588748    0.4987772572
#> 616_LOC100767716       0.9996591      0.93998435    0.4945708006
#> 617_Hnrnpc             0.4335248      0.53918279    0.9456073392
#> 618_Mphosph10          0.5731500      0.58902875    0.9273415365
#> 619_Eif3b              0.4884807      0.57340016    0.6724947438
#> 620_Emd                0.5975031      0.84244807    0.4987772572
#> 621_Txlng              1.0000000      1.00000000    1.0000000000
#> 622_Prpf4b             0.7391776      0.65495110    0.5103974918
#> 623_Rlim               0.9862244      0.81115498    0.6241079841
#> 624_Eef1b2             0.9996591      0.84503288    0.8094385315
#> 625_Def6               0.4335248      0.97469270    0.3344749025
#> 626_LOC100765020       0.9996591      0.95016268    0.3989353266
#> 627_U2surp             0.9937178      0.94004182    0.4919700587
#> 628_Elf2               0.9926304      0.93144523    0.4987772572
#> 629_Slc1a5             0.7528198      0.97469270    0.5696238511
#> 630_NA                 0.8540623      0.65495110    0.9252842391
#> 631_Tfg                0.7391776      0.97469270    0.9529197807
#> 632_Top2b              0.5842771      0.89228318    0.5319516213
#> 633_Pip4p2             0.9362984      0.90047417    0.9672759026
#> 634_Cdc42ep1           1.0000000      1.00000000    1.0000000000
#> 635_Hsph1              0.9937178      0.89453090    1.0000000000
#> 636_Twf1               0.9085372      1.00000000    1.0000000000
#> 637_Nbn                0.5975031      0.84244807    0.8377870644
#> 638_Psmd4              0.9206541      0.39726695    0.1492654263
#> 639_Bap1               0.9996591      0.61840397    0.3147010633
#> 640_Mepce              0.5842771      0.99844582    0.9863023699
#> 641_Mideas             1.0000000      1.00000000    1.0000000000
#> 642_LOC100759640       0.5731500      0.39726695    0.6352867542
#> 643_Epb41l2            0.6738408      0.93309878    0.7589135521
#> 644_Sav1               0.5731500      0.44694320    0.4014644565
#> 645_Prpf4b             0.5975031      0.57340016    0.4834221588
#> 646_Gnas               0.7619059      0.78620962    0.4823171066
#> 647_Mllt1              0.7585314      0.97469270    0.4987772572
#> 648_Poldip3            0.4335248      0.90047417    0.5655994197
#> 649_Aldoa              0.7391776      0.57340016    0.4830442104
#> 650_Rbbp8              0.9206541      0.89453090    0.9914761578
#> 651_LOC113834282       0.5731500      0.44694320    0.9071163669
#> 652_Gys1               0.9996591      0.93450401    0.9456073392
#> 653_Hnrnpc             0.8577026      0.51929656    0.3458807846
#> 654_Vps35              0.6445819      0.39726695    0.1870155117
#> 655_Miga2              0.9738333      0.81115498    0.8474656873
#> 656_Epb41l2            0.9996591      0.96609855    0.9973033256
#> 657_Tob2               1.0000000      1.00000000    1.0000000000
#> 658_Lamtor1            0.9996591      0.89188709    0.9984999026
#> 659_LOC100759640       0.6473547      0.65495110    0.6606817456
#> 660_Epb41l2            1.0000000      1.00000000    1.0000000000
#> 661_Rlim               1.0000000      1.00000000    0.9174607580
#> 662_Gys1               0.4884807      0.69361567    0.2758043564
#> 663_LOC100750437       0.9085372      0.84244807    0.3175586501
#> 664_NA                 0.7528198      0.50159122    0.1163090667
#> 665_Nbn                0.9996591      0.63814186    0.6956285583
#> 666_Tyw3               0.9085372      0.96685495    0.6796661112
#> 667_Gas2l1             0.5975031      0.92644080    0.3458807846
#> 668_Fus                0.9173402      1.00000000    1.0000000000
#> 669_Prpf38b            0.9937178      0.89453090    0.4945708006
#> 670_Calu               1.0000000      1.00000000    1.0000000000
#> 671_Rras2              0.9361973      0.91299684    0.3390104519
#> 672_Prpf4b             0.4884807      0.41581591    0.8967754084
#> 673_Nelfa              0.5731500      0.88113882    0.7036922531
#> 674_LOC100754077       0.5536475      0.89228318    0.6992386870
#> 675_Rbm28              0.8711111      0.99844582    0.8004961183
#> 676_Nsfl1c             0.6317540      0.86446114    0.8726332264
#> 677_Rnf126             0.3686822      0.81576416    0.9906710175
#> 678_Eme1               0.6471239      0.75726680    0.9174607580
#> 679_Nbn                0.4335248      0.94208942    0.1163090667
#> 680_Eif4ebp2           0.8540623      0.81647594    0.2738360002
#> 681_Wee1               0.5536475      0.42243195    0.1492654263
#> 682_Prpf38b            0.3045423      0.34703515    0.8094385315
#> 683_Luzp1              0.5975031      0.94478306    0.6164629145
#> 684_Gas2l1             0.7976077      0.39726695    0.1170839101
#> 685_Pdcd11             0.5975031      0.91266136    0.8053523228
#> 686_Chaf1b             0.3045423      0.22149083    0.6803315820
#> 687_Pycr1              0.4116137      0.39726695    0.6647202777
#> 688_Phf8               0.8789937      0.49995158    0.5105655186
#> 689_Raver1             0.9206541      0.97469270    0.8754200386
#> 690_Dbn1               0.7057658      0.72336026    0.8094385315
#> 691_Dut                0.8327078      0.41581591    0.1163090667
#> 692_Prpf4b             0.9996591      0.97469270    0.6790956907
#> 693_Prpf4b             0.9996591      0.89148887    0.6187384142
#> 694_Efs                0.4884807      0.88113882    0.6199397072
#> 695_NA                 0.5731500      0.97469270    0.3390104519
#> 696_Ppp2r5b            0.6485597      0.99844582    0.7448404619
#> 697_Caskin2            0.4335248      0.58893964    0.2007870458
#> 698_Arhgef40           1.0000000      1.00000000    1.0000000000
#> 699_Zyx                0.3533344      0.57340016    0.1163090667
#> 700_Mphosph10          0.9996591      0.58893964    0.1395399686
#> 701_LOC113833392       1.0000000      0.84503288    0.8474656873
#> 702_Cdc42ep1           0.9206541      0.74345593    0.4014644565
#> 703_Snrpa1             1.0000000      1.00000000    1.0000000000
#> 704_Ncbp1              0.7976077      0.86641680    0.9973033256
#> 705_Gas2l1             0.6904936      1.00000000    1.0000000000
#> 706_Gas2l1             0.9379294      0.83239038    0.8677435804
#> 707_Bap1               0.7239466      0.83239038    0.6724077964
#> 708_LOC100759640       0.9937178      0.69571122    0.7229251756
#> 709_Cherp              1.0000000      0.97469270    0.9909934070
#> 710_Nbn                0.0104809      0.16977096    0.5183821396
#> 711_LOC100759640       0.9996591      1.00000000    1.0000000000
#> 712_NA                 0.8540623      0.83239038    0.4014644565
#> 713_Eif3b              1.0000000      1.00000000    1.0000000000
#> 714_Miga2              0.4884807      0.84764375    0.6458262979
#> 715_Prpf4b             0.7765806      0.76519496    0.1542353843
#> 716_Dbn1               0.8443625      0.58893964    0.1272404753
#> 717_Ppp2r5b            0.7765806      0.89453090    0.3138611909
#> 718_Exosc9             1.0000000      1.00000000    0.1492654263
#> 719_Eif3b              0.9173402      0.96609855    0.4184815759
#> 720_Ripk2              0.9996591      0.97469270    0.6801679133
#> 721_Dlg1               0.6473547      0.95912280    0.8394579780
#> 722_N4bp1              0.7528198      0.57340016    0.5509669344
#> 723_Nudc               0.9996591      0.84503288    0.8744753671
#> 724_Znf367             1.0000000      0.83022221    1.0000000000
#> 725_Ring1              0.9307018      0.83239038    0.8187916031
#> 726_Snrpa1             0.7976077      0.92922311    0.9174607580
#> 727_U2surp             0.5975031      0.75726680    0.8620480770
#> 728_LOC100764225       0.7976077      0.78620962    0.9460896228
#> 729_Cdc42ep1           0.9996591      0.76519496    0.9446298214
#> 730_Znf385a            0.9150685      0.58893964    0.2441304415
#> 731_Ints1              0.9150685      0.58893964    0.2441304415
#> 732_LOC113833392       0.9150685      0.58893964    0.2441304415
#> 733_Lrch4              1.0000000      0.70094781    0.8515168795
#> 734_Ctdspl2            0.9996591      0.83239038    0.2384727084
#> 735_Prpf4b             0.8540623      0.83239038    0.5457472072
#> 736_Luzp1              0.7129418      0.84244807    0.9797810656
#> 737_Eif3b              0.1475220      0.77220390    0.3458807846
#> 738_Ptpn14             0.7528198      0.85102646    0.9701805996
#> 739_Rrp1               0.9937178      0.97808534    0.5959542038
#> 740_Lrrfip2            0.3045423      0.39726695    0.2738360002
#> 741_Nsfl1c             0.8913195      0.73710357    0.4463656788
#> 742_Ddx51              0.9996591      0.96609855    0.8474656873
#> 743_Prpf38b            0.5629666      0.52179637    0.2738360002
#> 744_Eef1b2             0.9307018      0.97469270    0.7962693514
#> 745_Znf385a            1.0000000      0.70094781    0.4314319337
#> 746_Map9               0.7239466      0.84674366    0.3669729209
#> 747_Rflnb              0.7239466      0.84674366    0.3669729209
#> 748_NA                 0.4884807      0.98226136    0.5103974918
#> 749_C1H12orf45         0.4884807      0.98226136    0.5103974918
#> 750_U2surp             0.7528198      0.91266136    0.3807330765
#> 751_Caskin2            0.8443625      0.97469270    0.5103974918
#> 752_Eri1               0.9173402      0.90047417    0.4358371316
#> 753_Gsk3b              0.9206541      0.98261381    0.6647202777
#> 754_LOC100766946       0.8789016      0.91266136    0.9240238739
#> 755_Cnpy3              1.0000000      1.00000000    0.2729399625
#> 756_Hnrnpc             0.7391776      0.83239038    1.0000000000
#> 757_Ptpn14             0.7391776      0.83239038    1.0000000000
#> 758_Slc7a11            0.6968896      0.90756350    0.5619464051
#> 759_Hnrnpc             0.7461617      0.95912280    0.5693079184
#> 760_Cdc37l1            0.4116137      0.32541914    0.3587160453
#> 761_LOC100768405       0.4884807      0.63814186    0.8311625987
#> 762_Rragc              0.4884807      0.63814186    0.8311625987
#> 763_LOC113834282       0.9996591      0.69493031    0.4028573087
#> 764_Fus                1.0000000      1.00000000    0.6629729761
#> 765_Ubxn1              0.4884807      0.84244807    0.9906710175
#> 766_Mmut               0.4884807      0.95236696    0.3911676182
#> 767_Pdcd11             0.4470790      0.70094781    0.6897028754
#> 768_LOC100757535       0.8353686      0.10784466    0.1406775311
#> 769_Eif3b              0.8039534      0.83239038    0.2895665722
#> 770_Rnf113a            0.9996591      0.89309712    0.3770190551
#> 771_Sytl4              0.9996591      0.58893964    0.3183586939
#> 772_Tlnrd1             0.8608189      0.96609855    0.9456073392
#> 773_H671_1g1131        0.9173402      0.74796187    0.5457472072
#> 774_Neurl1             0.7619059      0.61840397    0.3449452914
#> 775_Zyx                1.0000000      1.00000000    1.0000000000
#> 776_Ctdspl2            0.7619059      0.76519496    0.6272688737
#> 777_Chaf1b             0.9996591      0.57340016    0.1163090667
#> 778_Rragc              0.9937178      0.96609855    0.6724947438
#> 779_Srfbp1             0.4335248      0.51929656    0.1896719181
#> 780_Gys1               0.7585314      0.81647594    0.9292673283
#> 781_Usp15              0.6473547      0.99844582    0.5521486360
#> 782_Arhgef40           0.6485597      0.58893964    0.8120905602
#> 783_Gigyf1             1.0000000      1.00000000    1.0000000000
#> 784_Minar1             0.8540623      0.80301683    0.3458807846
#> 785_Dus2               0.9996591      0.57340016    0.4927937173
#> 786_Gatad2b            0.9173402      0.41853713    0.1662906723
#> 787_Eif5               1.0000000      0.27486659    0.1163090667
#> 788_Epb41l2            1.0000000      1.00000000    1.0000000000
#> 789_Arl6ip4            0.9173402      0.94517210    0.7759084747
#> 790_Plin4              0.9361973      0.91767997    0.9591634617
#> 791_Elf2               0.9996591      0.90654204    0.9090160656
#> 792_Plin4              1.0000000      0.57340016    0.3183586939
#> 793_Snip1              0.6933908      0.83239038    0.5423506258
#> 794_Txlng              0.5731500      0.84244807    0.7530007652
#> 795_LOC100769437       0.9307018      0.34703515    0.1163090667
#> 796_Caskin2            0.9996591      0.67796309    0.3183586939
#> 797_NA                 0.7987719      0.97462156    0.3557894948
#> 798_Synm               0.9428566      0.97469270    0.6803315820
#> 799_Synm               0.9206541      0.98168118    0.9690539503
#> 800_Ube2c              0.8447753      1.00000000    1.0000000000
#> 801_Sgtb               0.7528198      0.57340016    0.6275469607
#> 802_Prpf4b             0.5601500      0.84764375    0.5457472072
#> 803_Epb41l2            0.8241813      0.49257001    0.8453895336
#> 804_Mllt1              0.5975031      0.68895932    0.5778503376
#> 805_LOC100759640       1.0000000      1.00000000    1.0000000000
#> 806_Epb41l2            0.8478867      0.78620962    0.3814379204
#> 807_Znf280b            0.9996591      0.90654204    0.8473517850
#> 808_Kiaa1143           1.0000000      1.00000000    1.0000000000
#> 809_Gas2l1             0.5975031      0.32541914    0.8474656873
#> 810_Srp72              0.4884807      0.54140969    0.6460689701
#> 811_Tomm22             0.5975031      0.92922311    0.3498488866
#> 812_Psip1              0.4884807      0.97469270    0.3548546154
#> 813_Arhgef37           0.8478867      0.89453090    0.2796179498
#> 814_Bckdk              0.9615458      0.77660746    0.4927937173
#> 815_Strip1             0.8340945      0.74703751    0.1975058450
#> 816_Usp15              0.9173402      0.57955504    0.2738360002
#> 817_Ssr3               0.9996591      0.63814186    0.3243591041
#> 818_Strip1             0.7619059      0.70094781    0.3154653432
#> 819_Eif3b              0.4884807      0.98261381    0.1761106619
#> 820_U2surp             0.5682360      0.84518166    0.2102464538
#> 821_Bend3              0.9485913      0.57340016    0.2441304415
#> 822_Rps10              1.0000000      1.00000000    1.0000000000
#> 823_Rpl23a             0.8340945      0.34703515    0.1395399686
#> 824_Nbn                0.7992518      0.57340016    0.3879479941
#> 825_Rpap3              1.0000000      0.88113882    0.9460896228
#> 826_LOC100759640       0.9937178      0.89979241    0.2928150105
#> 827_Ric8a              0.9240383      0.80301683    0.2738360002
#> 828_Hsph1              0.9937178      0.66216815    0.4188878049
#> 829_LOC100759640       0.9409124      0.84764375    0.9273415365
#> 830_LOC100757535       1.0000000      1.00000000    1.0000000000
#> 831_Gigyf1             0.6507661      0.57955504    0.1272404753
#> 832_Dbn1               1.0000000      1.00000000    1.0000000000
#> 833_Snrk               1.0000000      1.00000000    1.0000000000
#> 834_Prpf38b            0.7224935      0.59567667    0.3458807846
#> 835_LOC100766868       1.0000000      1.00000000    1.0000000000
#> 836_LOC100766868       0.9206541      0.96609855    0.1732450585
#> 837_Wbp11              0.7671830      0.95236696    0.5619464051
#> 838_Rusc2              0.7976077      0.91629651    0.5103974918
#> 839_Eif3b              0.8711111      0.80301683    0.7764469293
#> 840_Ptpn14             0.9996591      0.89148887    0.6928203220
#> 841_Rlim               0.9996591      0.65495110    0.6654892405
#> 842_Ints1              0.6473547      0.84244807    0.6359380798
#> 843_Chaf1b             0.9206541      0.97469270    0.2496049280
#> 844_Dlg1               0.7239466      1.00000000    1.0000000000
#> 845_Lamtor1            0.7619059      0.76519496    0.2738360002
#> 846_Tab1               0.7671830      0.72336026    0.3558907560
#> 847_Dbn1               0.9362893      0.44694320    0.1896719181
#> 848_Psip1              0.9362893      0.44694320    0.1896719181
#> 849_Dbn1               0.9362893      0.84244807    0.4945708006
#> 850_Pabpc1             0.9524993      0.90047417    0.6987798252
#> 851_Hnrnpc             0.9361973      0.61641970    0.6595569957
#> 852_Emd                0.8049314      0.76519496    0.5319516213
#> 853_LOC100764225       0.8540623      0.99844582    0.8085069881
#> 854_Nup50              0.9442833      0.97469270    0.6619954261
#> 855_Ctcf               0.9173402      0.89883798    0.8444428354
#> 856_Raly               0.9307018      0.99844582    0.5103974918
#> 857_Bard1              0.9996591      0.84138511    0.6595569957
#> 858_Ptpn14             0.8540623      0.91548160    0.6262551451
#> 859_LOC100757535       0.5536475      0.62588442    0.1691476201
#> 860_Psmd2              1.0000000      1.00000000    1.0000000000
#> 861_Junb               0.9937178      0.81938063    0.3444986722
#> 862_C1qbp              0.8241813      0.73046398    0.5696238511
#> 863_Lrch4              0.8327078      0.89453090    0.6353934056
#> 864_CUNH14orf93        1.0000000      1.00000000    1.0000000000
#> 865_U2surp             0.9996224      0.63315255    0.3879479941
#> 866_Raly               1.0000000      1.00000000    1.0000000000
#> 867_LOC100774417       0.5842771      0.57340016    0.6803315820
#> 868_Srp72              0.9085372      1.00000000    1.0000000000
#> 869_LOC100764225       0.9409124      0.91266136    1.0000000000
#> 870_Morf4l2            0.7976077      0.58893964    0.3669755656
#> 871_CUNH9orf40         0.5975031      0.41581591    0.2441304415
#> 872_Gas2l1             0.9996591      0.70194592    0.8179717078
#> 873_Atp5pf             0.9996591      0.73152943    0.7126478180
#> 874_Lrrfip2            0.7976077      0.59567667    0.7962693514
#> 875_Prpf4b             1.0000000      0.84244807    0.2738360002
#> 876_Top2b              0.9996591      0.75726680    0.1542353843
#> 877_Mepce              0.4335248      0.75726680    0.8085069881
#> 878_Ptpn14             0.5731500      0.89453090    0.8474656873
#> 879_Dnajc25            0.6473547      0.90654204    0.1492654263
#> 880_Cbx8               0.9361973      0.57340016    0.2738360002
#> 881_Synm               0.6982773      0.96609855    0.5103974918
#> 882_Def6               1.0000000      1.00000000    1.0000000000
#> 883_Gys1               0.9996591      0.92123506    0.9773291621
#> 884_Luzp1              0.9996591      0.84764375    0.3437425283
#> 885_Synm               0.8241813      0.65688382    0.4655348841
#> 886_Snip1              0.7239466      0.39726695    0.2738360002
#> 887_Top2b              0.5601500      0.75726680    0.9273415365
#> 888_NA                 0.8566323      0.84518166    1.0000000000
#> 889_Trim35             1.0000000      1.00000000    1.0000000000
#> 890_Znf385a            0.9996591      0.73046398    0.4209451907
#> 891_Chaf1b             0.6971286      0.84244807    0.6801679133
#> 892_Abcf1              1.0000000      1.00000000    1.0000000000
#> 893_Pdcd11             0.5601500      0.68895932    0.8744753671
#> 894_Dlg1               0.9361973      0.68895932    0.4021681980
#> 895_Dbn1               0.9996224      0.90047417    0.8041289506
#> 896_LOC100752363       0.9173402      0.89453090    0.9562143923
#> 897_Ppp4r3a            0.9996591      0.89453090    0.3955081305
#> 898_Gas2l1             0.7239466      0.57340016    0.4945708006
#> 899_Mtmr10             1.0000000      1.00000000    1.0000000000
#> 900_Cyld               0.7528198      0.97469270    0.9642694022
#> 901_NA                 0.9937178      0.63814186    0.8138575422
#> 902_Rnf113a            0.9996591      0.99844582    0.4945708006
#> 903_Nelfa              1.0000000      0.99844582    0.7068411866
#> 904_Zkscan1            0.9996591      0.99844582    0.7448404619
#> 905_Chaf1b             0.9206541      0.77660746    0.2441304415
#> 906_Eif3b              0.7585314      0.96609855    0.2796179498
#> 907_Top2b              0.8312584      0.77660746    0.4823171066
#> 908_Chaf1b             1.0000000      1.00000000    0.3458807846
#> 909_Epb41l2            1.0000000      1.00000000    1.0000000000
#> 910_C3H11orf58         1.0000000      1.00000000    0.3154653432
#> 911_Top2b              0.8540623      0.89228318    0.4987772572
#> 912_Wee1               0.5975031      0.74590123    0.9413836018
#> 913_Raly               0.9173402      0.57340016    0.4184815759
#> 914_H671_1g2680        0.9996591      0.89148887    0.6926345811
#> 915_Eef1b2             0.9085372      0.72336026    0.5103974918
#> 916_Gas2l1             0.7619059      0.47674732    0.1163090667
#> 917_Epb41l2            0.8540623      0.40322706    0.1395399686
#> 918_Rpl23a             0.8540623      0.84503288    0.3342221489
#> 919_Chmp2b             1.0000000      0.34703515    0.3669729209
#> 920_Lrrfip2            0.7585314      0.99844582    0.8611757851
#> 921_Aldoa              0.9634942      0.70868902    0.2090589517
#> 922_Cby1               0.4884807      0.49995158    0.5805572519
#> 923_LOC100759640       0.9426037      0.57340016    0.4187969375
#> 924_Rbm28              0.3045423      0.32541914    0.8053523228
#> 925_Skiv2l             1.0000000      1.00000000    0.1761106619
#> 926_Ints1              0.9996591      0.89797917    0.8474656873
#> 927_Ehd1               0.8617286      0.39726695    0.4889794954
#> 928_Nr2f6              0.8093895      0.74503582    0.8053523228
#> 929_Top2b              0.9206541      0.83239038    0.6803315820
#> 930_Lrrfip2            0.9937178      0.86446114    0.6320453800
#> 931_Pip4p2             0.9996591      0.76519496    0.3458807846
#> 932_Srp72              0.7789744      0.95648901    0.7448404619
#> 933_Mtmr9              0.9996591      0.57340016    0.1552486432
#> 934_Gigyf1             0.9996591      0.80301683    0.7347336557
#> 935_Rbm7               0.9173402      0.62588442    0.2738360002
#> 936_LOC100773565       0.8156605      0.77660746    0.4919700587
#> 937_Trim35             0.8327078      0.59876431    0.3676390944
#> 938_Cbx8               0.9996591      0.92123506    0.4830442104
#> 939_Rplp0              0.9918258      1.00000000    1.0000000000
#> 940_Aldoa              0.9996591      0.97469270    0.9460896228
#> 941_NA                 0.6213394      0.50159122    0.3277644801
#> 942_Zyx                1.0000000      0.94208942    0.7566755605
#> 943_Psip1              0.5731500      0.57340016    0.5532033310
#> 944_Slc7a11            0.5682360      0.76519496    0.8179717078
#> 945_Miga2              0.5805890      0.39726695    0.1552486432
#> 946_Arhgef6            0.7886993      0.96609855    0.2988269514
#> 947_Dlgap4             0.9937178      0.73152943    0.8873517778
#> 948_Ampd2              1.0000000      1.00000000    1.0000000000
#> 949_Luzp1              0.5445966      0.80301683    0.9174607580
#> 950_Camlg              1.0000000      1.00000000    1.0000000000
#> 951_Pfkfb3             0.8705424      0.85102646    0.3458807846
#> 952_NA                 1.0000000      1.00000000    1.0000000000
#> 953_Raly               1.0000000      1.00000000    1.0000000000
#> 954_Kiaa1143           0.5731500      0.41581591    0.1395399686
#> 955_Bcar1              0.9937178      0.62588442    0.4014644565
#> 956_Gatad2b            0.4884807      0.57955504    0.5319516213
#> 957_Eif4ebp2           0.3045423      0.32541914    0.4358371316
#> 958_Fam76b             0.4884807      0.57340016    0.3458807846
#> 959_Camlg              0.8539081      0.67937672    0.3770190551
#> 960_LOC100754077       0.9996591      0.63814186    0.1343299522
#> 961_NA                 0.9996591      0.29247726    0.0237123566
#> 962_Epb41l2            0.8754165      0.99480848    0.4945708006
#> 963_Ankrd34a           1.0000000      1.00000000    1.0000000000
#> 964_Zc3h15             0.9996591      0.78620962    0.2610754681
#> 965_Def6               0.4335248      0.92434829    0.4034772289
#> 966_Srsf6              1.0000000      1.00000000    0.5502559641
#> 967_H671_4g11480       0.5975031      0.74796187    0.7613779935
#> 968_Top2b              1.0000000      1.00000000    0.9464121243
#> 969_LOC100769471       0.9996591      0.83239038    0.4889794954
#> 970_Raver1             0.4884807      0.67222253    0.4987772572
#> 971_Etv3               1.0000000      0.96707635    1.0000000000
#> 972_Psd                0.9996591      0.89797917    0.8042667469
#> 973_Usp15              1.0000000      1.00000000    0.7211170537
#> 974_Nol7               0.9159876      0.74938527    0.8085069881
#> 975_Stk38              0.9996591      0.58893964    0.1691476201
#> 976_Smim13             0.4937630      0.76934373    0.3812768780
#> 977_Etv3               0.4884807      0.10784466    0.1374908882
#> 978_Synm               0.4884807      0.39726695    0.6449466883
#> 979_Pwp1               0.9996591      0.75726680    0.4111144282
#> 980_Fus                0.9150685      0.89309712    0.9520677088
#> 981_Junb               0.8739042      0.90581890    0.7589135521
#> 982_Phf8               0.7482648      0.70094781    0.9164066845
#> 983_Nelfa              1.0000000      1.00000000    1.0000000000
#> 984_Prpf4b             1.0000000      1.00000000    1.0000000000
#> 985_Abraxas1           0.4470790      0.84244807    0.6840169429
#> 986_Prpf4b             0.9996591      1.00000000    1.0000000000
#> 987_Raver1             1.0000000      1.00000000    1.0000000000
#> 988_Caap1              0.9159876      0.95648901    0.9016849838
#> 989_Rpap3              0.9937178      0.94475915    0.9174607580
#> 990_Hsph1              0.9937178      0.74796187    0.4028573087
#> 991_LOC100750437       1.0000000      1.00000000    1.0000000000
#> 992_Mepce              0.8156605      0.98300959    0.6582681907
#> 993_Efs                0.9937178      0.68895932    0.8591498257
#> 994_Epb41l2            0.9996591      0.90756350    0.9628343032
#> 995_Abcf1              0.4884807      0.50159122    0.4823171066
#> 996_NA                 0.9085372      0.96609855    0.6724947438
#> 997_Eif4ebp2           1.0000000      1.00000000    1.0000000000
#> 998_Pfkfb3             0.6295160      0.91266136    0.8120905602
#> 999_Hnrnpc             0.7224935      0.44694320    0.3587160453
#> 1000_Psmd2             0.9307018      0.50159122    0.4967507613
#>                  120_vs_neighbors
#> 1_Top2b                0.99332811
#> 2_NA                   0.65273244
#> 3_Snip1                1.00000000
#> 4_Tomm34               0.86859820
#> 5_Pus3                 0.96317292
#> 6_Ints1                1.00000000
#> 7_Mlh1                 0.95822477
#> 8_LOC100750437         0.78837220
#> 9_Pabpc1               0.81550950
#> 10_Top2b               0.87684327
#> 11_Gorasp1             0.57624609
#> 12_Ints1               1.00000000
#> 13_Syvn1               0.98272300
#> 14_Znf280b             0.71378983
#> 15_Mrnip               0.87125683
#> 16_Rragc               0.87684327
#> 17_Gorasp1             0.74897796
#> 18_Tomm34              1.00000000
#> 19_LOC100757430        0.86650073
#> 20_Ubxn1               0.89691748
#> 21_H671_1g1131         1.00000000
#> 22_Luzp1               0.84433926
#> 23_Efs                 0.96309797
#> 24_Mta2                0.59299794
#> 25_Nedd1               0.53998241
#> 26_Gigyf1              0.95401692
#> 27_Myh9                0.96317292
#> 28_Caskin2             0.61804675
#> 29_Papolg              0.53998241
#> 30_Tfg                 1.00000000
#> 31_Rpl34               0.89691748
#> 32_Mideas              0.53647630
#> 33_Gys1                0.53998241
#> 34_Arhgef6             0.53998241
#> 35_Ctdspl2             0.99040703
#> 36_Ptpn14              0.76552820
#> 37_Raly                1.00000000
#> 38_Znhit3              0.82285588
#> 39_LOC113833392        1.00000000
#> 40_Luc7l3              0.99040703
#> 41_Rplp0               1.00000000
#> 42_Gys1                0.65792365
#> 43_Rpl22l1             0.61804675
#> 44_Eif3b               0.80741720
#> 45_Med26               0.71891081
#> 46_Mepce               0.71378983
#> 47_Pdcd11              0.97937571
#> 48_Twf1                1.00000000
#> 49_LOC100759640        1.00000000
#> 50_Wrnip1              0.76463145
#> 51_Poldip3             0.89691748
#> 52_Ampd2               0.64034789
#> 53_Mea1                0.81456911
#> 54_Dbn1                0.89691748
#> 55_Snip1               0.61127133
#> 56_Srsf6               0.84680256
#> 57_LOC113834282        0.11563729
#> 58_Map9                0.69586387
#> 59_Cdc42ep1            0.53998241
#> 60_Poldip3             0.96317292
#> 61_LOC100764225        0.54836325
#> 62_Epb41l2             1.00000000
#> 63_H671_4g11480        0.76463145
#> 64_Nbn                 1.00000000
#> 65_U2surp              1.00000000
#> 66_Gigyf1              0.59107576
#> 67_NA                  0.57624609
#> 68_Luc7l3              0.57624609
#> 69_LOC100752363        0.86650073
#> 70_Ampd2               0.53647630
#> 71_LOC100759640        0.53647630
#> 72_Stam                0.89691748
#> 73_Nsfl1c              0.71891081
#> 74_Pfkfb3              0.96317292
#> 75_Rad23a              0.61804675
#> 76_Elf2                0.97937571
#> 77_Crem                0.87125683
#> 78_Rragc               0.82801076
#> 79_Lrrfip2             0.89691748
#> 80_Zyx                 0.46858598
#> 81_Lrrfip2             1.00000000
#> 82_Gatad2b             0.75084189
#> 83_Bcar1               0.82531804
#> 84_Ehd1                0.61127133
#> 85_LOC113834282        0.71891081
#> 86_Tmem230             0.76552820
#> 87_Ncbp1               0.89691748
#> 88_Mllt1               1.00000000
#> 89_Stk17b              0.53647630
#> 90_Dlgap4              0.76552820
#> 91_Papolg              0.67058556
#> 92_Cyld                0.61127133
#> 93_Gigyf1              0.99790083
#> 94_Lrrfip2             0.95365528
#> 95_Lrrfip2             0.76463145
#> 96_Rlim                0.89977838
#> 97_Eif3b               0.95365528
#> 98_Mphosph10           0.65273244
#> 99_Gatad2b             0.96317292
#> 100_Srsf6              0.58004623
#> 101_Zyx                0.80007140
#> 102_Mphosph10          0.59107576
#> 103_Psip1              0.89691748
#> 104_Fbl                1.00000000
#> 105_H671_1g2680        0.98272300
#> 106_Sgtb               0.99040703
#> 107_Gnl3               0.96317292
#> 108_Eif3b              0.69265742
#> 109_Serpinb1           0.61127133
#> 110_N4bp1              0.79370312
#> 111_Snip1              0.82801076
#> 112_Psip1              0.86686270
#> 113_Mlh1               0.87684327
#> 114_Bsg                0.80519954
#> 115_Tnpo1              1.00000000
#> 116_H671_1g2680        0.53998241
#> 117_Cbx8               0.71378983
#> 118_Mideas             0.53998241
#> 119_Mideas             0.96317292
#> 120_Dcun1d3            0.96317292
#> 121_Dlg1               0.96317292
#> 122_Rad23a             0.89691748
#> 123_Srsf6              0.87125683
#> 124_Stx7               0.84433926
#> 125_Pdcd11             0.92281298
#> 126_Kiaa1958           0.55994835
#> 127_Pwp1               0.67058556
#> 128_Txlng              0.76552820
#> 129_Junb               0.86650073
#> 130_LOC100759640       0.89691748
#> 131_Dbn1               0.61127133
#> 132_Top2b              0.59107576
#> 133_Rusc2              0.80007140
#> 134_NA                 1.00000000
#> 135_LOC113837251       0.97937571
#> 136_Fam76b             0.98272300
#> 137_Ptpn14             0.76463145
#> 138_Chmp4b             0.98272300
#> 139_Prpf4b             0.61127133
#> 140_Eif3b              1.00000000
#> 141_Nsfl1c             0.87684327
#> 142_Pdlim7             0.95731924
#> 143_Rnf113a            0.80822413
#> 144_Epb41l2            0.71891081
#> 145_Hnrnpc             0.75084189
#> 146_LOC113834282       0.89691748
#> 147_Plekho2            0.86319749
#> 148_Med26              0.98272300
#> 149_Arhgef40           0.95013062
#> 150_NA                 0.53647630
#> 151_Phf8               0.78837220
#> 152_Minar1             1.00000000
#> 153_H671_21690         0.87435458
#> 154_Arhgef40           0.70198358
#> 155_Chaf1b             0.84519728
#> 156_Prpf4b             0.89691748
#> 157_Znf367             0.65273244
#> 158_Luzp1              0.87435458
#> 159_LOC113833882       0.81550950
#> 160_Hnrnpc             0.61127133
#> 161_Mepce              0.76552820
#> 162_Ubxn1              0.96317292
#> 163_Mllt1              0.71891081
#> 164_Chaf1b             1.00000000
#> 165_Raly               0.89691748
#> 166_Gas2l1             0.74161192
#> 167_Dlg1               0.46858598
#> 168_Hoxc10             1.00000000
#> 169_Gigyf1             1.00000000
#> 170_Luzp1              0.89691748
#> 171_Srp72              0.81284520
#> 172_LOC100771461       1.00000000
#> 173_Chaf1b             0.65273244
#> 174_C3H11orf58         0.73477355
#> 175_Pdcd11             0.99040703
#> 176_Psip1              0.81456911
#> 177_Prpf4b             0.53998241
#> 178_Rnf113a            0.61127133
#> 179_Irf3               0.87684327
#> 180_Smim13             0.89691748
#> 181_Gnl3               0.89691748
#> 182_Psma5              0.55097037
#> 183_Ptpn14             0.97937571
#> 184_Prpf4b             0.78837220
#> 185_Top2b              0.91498011
#> 186_Prpf38b            0.89691748
#> 187_Epb41l2            1.00000000
#> 188_Eif3b              0.95731924
#> 189_Hnrnpc             0.67058556
#> 190_LOC100758278       0.96425465
#> 191_Prpf4b             0.76463145
#> 192_Caskin2            0.31331155
#> 193_LOC100752363       0.96317292
#> 194_Septin6            0.61127133
#> 195_Max                0.59299794
#> 196_Mid1ip1            0.88019189
#> 197_NA                 0.76552820
#> 198_Hsph1              0.67058556
#> 199_Nol7               0.87684327
#> 200_Raly               0.96317292
#> 201_Smim13             0.85184664
#> 202_LOC100757535       0.87684327
#> 203_Net1               0.95831218
#> 204_LOC100754077       0.99821395
#> 205_Snip1              0.76552820
#> 206_Hnrnpc             0.01545919
#> 207_Ldlrap1            0.89691748
#> 208_Luzp1              0.61127133
#> 209_Rpl26              0.87684327
#> 210_Epb41l2            0.78427438
#> 211_Znf367             0.54346611
#> 212_Dlgap4             0.61127133
#> 213_Plekho2            1.00000000
#> 214_Zpr1               0.72214605
#> 215_Dlgap4             1.00000000
#> 216_Def6               0.76552820
#> 217_Eif4ebp2           0.61127133
#> 218_Eef1b2             0.89691748
#> 219_Rad23a             0.80007140
#> 220_Morf4l2            0.61127133
#> 221_Arhgef40           0.82069736
#> 222_NA                 0.87684327
#> 223_LOC100773565       0.81550950
#> 224_Dus2               0.95897789
#> 225_Pip4p2             0.92872331
#> 226_Top2b              0.87684327
#> 227_Znf280b            0.76463145
#> 228_Pdcd11             0.71378983
#> 229_Bckdk              0.76552820
#> 230_Arhgef40           0.53647630
#> 231_Mepce              0.53998241
#> 232_Ccnd3              0.89691748
#> 233_Phf8               0.87684327
#> 234_H671_1g2680        1.00000000
#> 235_Ell                0.96317292
#> 236_U2surp             0.53998241
#> 237_Rps10              0.86686270
#> 238_Ctdspl2            0.46858598
#> 239_Top2b              0.61947174
#> 240_Msantd3            0.87684327
#> 241_Fam76b             0.96317292
#> 242_Ppp4r3a            0.64518630
#> 243_Gpatch4            1.00000000
#> 244_Nudc               0.99040703
#> 245_Nol7               0.96425465
#> 246_Plekho2            0.86650073
#> 247_Prpf4b             0.53647630
#> 248_Mta2               1.00000000
#> 249_U2surp             0.79131251
#> 250_Ubxn1              0.71891081
#> 251_Rlim               0.97937571
#> 252_Atat1              0.76552820
#> 253_Ubxn1              0.96317292
#> 254_H671_1g2680        0.76552820
#> 255_eIF2aK2            0.68492351
#> 256_Skiv2l             0.71378983
#> 257_Rpl28              0.82260728
#> 258_LOC100759640       0.53998241
#> 259_Gatad2b            0.80822413
#> 260_NA                 0.89691748
#> 261_Gprasp1            0.82285588
#> 262_Luzp1              0.88372742
#> 263_Slc1a5             0.81550950
#> 264_LOC113834282       0.71378983
#> 265_Srsf6              0.98272300
#> 266_Cdc42ep1           0.53998241
#> 267_Net1               1.00000000
#> 268_Caskin2            0.55406551
#> 269_LOC100759640       0.82243218
#> 270_Mideas             0.91879652
#> 271_Luzp1              1.00000000
#> 272_Emd                0.91734510
#> 273_Plpp6              0.81907261
#> 274_LOC100759640       0.74897796
#> 275_Rps7               0.74897796
#> 276_Fkbp1a             0.97937571
#> 277_Gatad2b            0.53998241
#> 278_Znf385a            0.65792365
#> 279_Arhgef6            0.84433926
#> 280_Slirp              0.85184664
#> 281_Skiv2l             0.71378983
#> 282_H671_21690         0.89691748
#> 283_Kat8               0.95822477
#> 284_Nkap               0.89691748
#> 285_Gsk3b              0.95731924
#> 286_Ints1              0.74897796
#> 287_Gas2l1             1.00000000
#> 288_LOC100759640       0.81456911
#> 289_Top2b              0.65273244
#> 290_Kif20b             0.87684327
#> 291_Phf8               1.00000000
#> 292_Snip1              0.61127133
#> 293_Gsk3b              0.96317292
#> 294_Caskin2            0.96317292
#> 295_C3H11orf58         1.00000000
#> 296_Lrch4              0.85184664
#> 297_LOC113834282       0.80822413
#> 298_LOC100750407       0.87125683
#> 299_LOC113833392       0.89691748
#> 300_LOC113833882       0.76463145
#> 301_Ldlrap1            0.87684327
#> 302_Wee1               0.74897796
#> 303_Caap1              1.00000000
#> 304_Eif4ebp2           0.81284520
#> 305_Ripk2              0.86563801
#> 306_Srp72              0.59299794
#> 307_Taok2              1.00000000
#> 308_Nr2f6              0.82285588
#> 309_Arhgef40           0.93783336
#> 310_Gys1               0.76552820
#> 311_Dlg1               0.86650073
#> 312_Vapb               0.78483347
#> 313_LOC100757535       1.00000000
#> 314_Mkrn2              0.99040703
#> 315_Eif3b              0.81456911
#> 316_Isyna1             0.61804675
#> 317_Prpf4b             0.81912718
#> 318_LOC113833882       0.98272300
#> 319_Lrch4              0.96425465
#> 320_Dbn1               0.93345068
#> 321_Abcf1              0.71891081
#> 322_Ints1              0.76463145
#> 323_C3H11orf58         0.78837220
#> 324_Psma5              1.00000000
#> 325_Fundc1             0.98272300
#> 326_Papolg             0.78837220
#> 327_Mideas             0.61198963
#> 328_Ubxn1              0.89691748
#> 329_Synm               1.00000000
#> 330_Arhgef6            1.00000000
#> 331_Ptpn14             1.00000000
#> 332_Pgrmc1             1.00000000
#> 333_Myh9               0.89691748
#> 334_Etv3               0.81550950
#> 335_Ip6k1              0.59299794
#> 336_Luzp1              0.96425465
#> 337_Ptpn14             0.78837220
#> 338_Caskin2            0.96317292
#> 339_Chaf1b             0.61127133
#> 340_Ubxn1              1.00000000
#> 341_Ube2c              0.57624609
#> 342_Gins2              0.98272300
#> 343_Nlgn2              0.71425490
#> 344_Nf2                0.76463145
#> 345_Pip4p2             0.90249140
#> 346_Emd                0.89691748
#> 347_Top2b              0.85184664
#> 348_Trim35             0.99790083
#> 349_NA                 0.81550950
#> 350_NA                 1.00000000
#> 351_Mideas             0.61804675
#> 352_Gas2l1             0.96317292
#> 353_Ampd2              0.99040703
#> 354_Calu               0.85779284
#> 355_Fam76b             0.86650073
#> 356_Dlg1               1.00000000
#> 357_Srsf6              0.89691748
#> 358_Chaf1b             0.89691748
#> 359_Dbn1               0.96302284
#> 360_Tcf25              0.61804675
#> 361_Psip1              0.73733493
#> 362_Cnpy3              0.80822413
#> 363_LOC100759640       0.83332028
#> 364_Zyx                1.00000000
#> 365_Lrch4              0.86650073
#> 366_Bola1              1.00000000
#> 367_Znf385a            0.89691748
#> 368_Kif20b             0.87435458
#> 369_Ell                1.00000000
#> 370_Ell                1.00000000
#> 371_Srsf6              0.61127133
#> 372_Pwp1               0.87684327
#> 373_Def6               0.98272300
#> 374_Cbx8               1.00000000
#> 375_Ddx51              0.75084189
#> 376_Psip1              0.89691748
#> 377_Arhgef40           0.96317292
#> 378_Raly               0.53647630
#> 379_NA                 0.99040703
#> 380_Lrrfip2            0.70208246
#> 381_Gnl3               0.98272300
#> 382_Caskin2            0.81456911
#> 383_Rragc              0.89691748
#> 384_Caskin2            0.61127133
#> 385_Bcar1              0.98272300
#> 386_Homer3             0.95731924
#> 387_Luzp1              0.96317292
#> 388_N4bp1              1.00000000
#> 389_Ppp4r3a            0.91526116
#> 390_H671_1g2680        1.00000000
#> 391_Gnl3               0.61127133
#> 392_Top2b              0.71425490
#> 393_Oser1              0.98272300
#> 394_Snrk               0.57624609
#> 395_Kat8               1.00000000
#> 396_Raver1             1.00000000
#> 397_Pdcd11             0.95365528
#> 398_Rps20              0.76552820
#> 399_Bsg                1.00000000
#> 400_Raly               0.53998241
#> 401_Pdcd2              0.65273244
#> 402_Caskin2            0.61127133
#> 403_LOC100773571       0.53647630
#> 404_Papolg             0.89691748
#> 405_LOC100757535       0.84139140
#> 406_Caap1              0.53998241
#> 407_Psip1              0.73477355
#> 408_Dbn1               0.81456911
#> 409_Mta2               0.96317292
#> 410_Abcf1              1.00000000
#> 411_LOC100754108       0.97937571
#> 412_Slirp              0.61127133
#> 413_Nelfa              0.98891432
#> 414_Aggf1              0.99040703
#> 415_Bap1               0.86859820
#> 416_Luc7l3             0.82069736
#> 417_Rrp1               0.88019189
#> 418_Wrnip1             1.00000000
#> 419_NA                 0.87684327
#> 420_Abcf1              0.77944157
#> 421_Cluap1             0.89691748
#> 422_Hnrnpc             0.53647630
#> 423_Ptpn1              0.46858598
#> 424_Myh9               0.89691748
#> 425_U2surp             0.76552820
#> 426_NA                 0.98272300
#> 427_Arhgef40           0.86859820
#> 428_Chaf1b             0.94766133
#> 429_Prpf4b             1.00000000
#> 430_Epb41l2            1.00000000
#> 431_Eif3b              0.76552820
#> 432_Isyna1             0.55097037
#> 433_U2surp             0.99040703
#> 434_LOC100765020       0.53998241
#> 435_Arhgef6            0.53647630
#> 436_Ptpn1              0.70202424
#> 437_Prpf4b             0.99040703
#> 438_Rpl35a             1.00000000
#> 439_Prpf4b             0.99404010
#> 440_Zyx                0.89691748
#> 441_Dbn1               0.87684327
#> 442_Chaf1b             0.74897796
#> 443_LOC113834282       0.76463145
#> 444_Gpsm2              0.81550950
#> 445_LOC100757535       0.82254097
#> 446_Cfap410            0.98891432
#> 447_Epb41l2            0.53998241
#> 448_Ncbp1              0.76463145
#> 449_Pacsin1            0.53647630
#> 450_Cstf2              0.87684327
#> 451_LOC100769437       0.84745243
#> 452_eIF2aK2            0.76552820
#> 453_Kiaa1191           0.99040703
#> 454_Mepce              0.61127133
#> 455_Cbx8               0.97937571
#> 456_Eed                0.53998241
#> 457_Cdc42ep1           0.89691748
#> 458_Lrrfip2            0.57624609
#> 459_Pacsin1            0.70453687
#> 460_Gpatch4            0.78427438
#> 461_Plin4              0.89691748
#> 462_NA                 0.95731924
#> 463_Snip1              0.43271982
#> 464_Cyld               0.61127133
#> 465_Plin4              0.98272300
#> 466_Twf1               1.00000000
#> 467_LOC113834282       1.00000000
#> 468_Snip1              0.87684327
#> 469_Ppp4r3a            1.00000000
#> 470_Psip1              0.71378983
#> 471_Dnajc5             0.89691748
#> 472_Phf8               0.46858598
#> 473_Bola1              0.87684327
#> 474_Cdc42ep1           0.86650073
#> 475_Eif4ebp2           0.80605441
#> 476_Prpf38b            0.89691748
#> 477_Klhl26             0.67058556
#> 478_Hsph1              0.96317292
#> 479_Snip1              0.80007140
#> 480_Caskin2            0.80519954
#> 481_Plpp6              0.98272300
#> 482_NA                 1.00000000
#> 483_Mlh1               0.96317292
#> 484_Gys1               0.99790083
#> 485_Tfg                0.92779889
#> 486_Arhgef6            0.65273244
#> 487_Mphosph10          0.65273244
#> 488_Hoxc10             0.83634429
#> 489_LOC100759640       0.67058556
#> 490_Arhgef40           0.99040703
#> 491_Dnajc5             0.98272300
#> 492_Tbc1d23            1.00000000
#> 493_Ubxn1              0.59422876
#> 494_Rab1a              0.89691748
#> 495_Eif3b              0.87684327
#> 496_Tceal8             0.80350549
#> 497_Dlgap4             0.78837220
#> 498_Smim13             0.53647630
#> 499_NA                 1.00000000
#> 500_Lrch4              0.46858598
#> 501_Bola1              0.59107576
#> 502_NA                 0.68324247
#> 503_Ptpn14             0.97937571
#> 504_LOC100759640       0.61127133
#> 505_Rps10              0.99040703
#> 506_Top2b              0.99040703
#> 507_Ssr3               0.61127133
#> 508_Homer3             0.61127133
#> 509_Phf8               0.87684327
#> 510_LOC100767716       0.86617797
#> 511_Xpa                0.53998241
#> 512_H671_21690         0.76463145
#> 513_LOC100769471       0.74897796
#> 514_Gas2l1             0.57624609
#> 515_Luzp1              0.53998241
#> 516_Gpbp1              0.57624609
#> 517_Gatad2b            1.00000000
#> 518_Gys1               0.89691748
#> 519_Top2b              0.89391002
#> 520_LOC100757535       1.00000000
#> 521_Lpcat4             0.98891432
#> 522_Arhgef6            0.65273244
#> 523_Cavin3             0.61127133
#> 524_Gpatch4            0.99333603
#> 525_Prpf38b            0.76463145
#> 526_Timm8a             0.53998241
#> 527_Cavin3             0.82069736
#> 528_Mkrn2              0.76552820
#> 529_Oser1              0.87684327
#> 530_Gsk3b              0.53998241
#> 531_Eef1b2             1.00000000
#> 532_Ampd2              0.97937571
#> 533_Lrrfip2            0.86650073
#> 534_Ring1              0.79131251
#> 535_Rlim               1.00000000
#> 536_LOC100759640       0.76552820
#> 537_LOC100759640       0.74897796
#> 538_Atp5pf             0.94766133
#> 539_Max                0.76552820
#> 540_Bap1               0.96309797
#> 541_Nsfl1c             0.53998241
#> 542_Prpf4b             0.89691748
#> 543_LOC100757535       0.87684327
#> 544_Mtmr10             0.66402033
#> 545_Hoxc10             1.00000000
#> 546_Trim35             0.61127133
#> 547_Eif4ebp2           0.92779889
#> 548_Dlgap4             1.00000000
#> 549_Gys1               0.86686270
#> 550_Sgtb               0.97937571
#> 551_Eri2               0.95365528
#> 552_Ccnd3              1.00000000
#> 553_Smim13             0.64380397
#> 554_Snrk               0.82069736
#> 555_Caskin2            0.70453687
#> 556_Pdcd11             1.00000000
#> 557_Pgam5              0.89691748
#> 558_Mphosph10          0.86650073
#> 559_Mideas             0.84433926
#> 560_Top2b              0.76463145
#> 561_LOC100763014       0.61117857
#> 562_Snip1              0.53998241
#> 563_Ubxn1              0.86859820
#> 564_LOC100750407       1.00000000
#> 565_Morf4l2            1.00000000
#> 566_Ctdspl2            0.89691748
#> 567_Cwf19l1            0.89691748
#> 568_Eef1b2             1.00000000
#> 569_C1H12orf45         1.00000000
#> 570_Znf367             0.53647630
#> 571_Ankrd34a           0.69265742
#> 572_Mllt11             0.99616460
#> 573_LOC100774792       0.88253951
#> 574_NA                 0.64270521
#> 575_Cbx8               0.86319749
#> 576_Bckdk              0.99040703
#> 577_Snip1              1.00000000
#> 578_Nsfl1c             1.00000000
#> 579_Gas2l1             0.65273244
#> 580_Nudc               0.46858598
#> 581_Epb41l2            1.00000000
#> 582_Mtmr6              0.74897796
#> 583_Znf668             0.95731924
#> 584_Hsph1              0.70198358
#> 585_LOC113834282       0.92669903
#> 586_Ctdspl2            1.00000000
#> 587_Foxf1              0.61127133
#> 588_Luzp1              0.82243218
#> 589_Xpa                0.88019189
#> 590_Psip1              0.57624609
#> 591_Rbm7               0.82148606
#> 592_Mtrex              0.91740623
#> 593_Arhgef40           1.00000000
#> 594_Plekho2            0.95731924
#> 595_Bckdk              0.76463145
#> 596_Dut                0.87684327
#> 597_Abcf1              0.87684327
#> 598_Txnl1              0.99500130
#> 599_Nudc               0.77944157
#> 600_Sh3gl1             0.99040703
#> 601_Gatad2b            0.82254097
#> 602_Homer3             0.89691748
#> 603_Septin6            0.89691748
#> 604_Smim13             0.85184664
#> 605_Arhgef40           0.81456911
#> 606_Rpl32              1.00000000
#> 607_Tomm34             0.96317292
#> 608_Mlh1               0.46858598
#> 609_Tbcc               1.00000000
#> 610_Eif3d              0.89691748
#> 611_Snrk               0.65273244
#> 612_Bckdk              0.97937571
#> 613_Wdr3               1.00000000
#> 614_LOC100757535       0.81153692
#> 615_Dlg1               0.99040703
#> 616_LOC100767716       0.70453687
#> 617_Hnrnpc             0.98272300
#> 618_Mphosph10          0.95365528
#> 619_Eif3b              0.95831218
#> 620_Emd                0.94326412
#> 621_Txlng              1.00000000
#> 622_Prpf4b             0.76552820
#> 623_Rlim               0.71891081
#> 624_Eef1b2             0.81907261
#> 625_Def6               0.71378983
#> 626_LOC100765020       0.61127133
#> 627_U2surp             0.81153692
#> 628_Elf2               1.00000000
#> 629_Slc1a5             0.99040703
#> 630_NA                 0.96317292
#> 631_Tfg                0.81153692
#> 632_Top2b              0.81550950
#> 633_Pip4p2             0.82203853
#> 634_Cdc42ep1           0.66219472
#> 635_Hsph1              1.00000000
#> 636_Twf1               1.00000000
#> 637_Nbn                0.95365528
#> 638_Psmd4              0.98272300
#> 639_Bap1               0.89691748
#> 640_Mepce              0.73741749
#> 641_Mideas             1.00000000
#> 642_LOC100759640       0.81550950
#> 643_Epb41l2            0.89691748
#> 644_Sav1               0.89691748
#> 645_Prpf4b             0.59107576
#> 646_Gnas               0.59107576
#> 647_Mllt1              0.67058556
#> 648_Poldip3            0.97937571
#> 649_Aldoa              0.61127133
#> 650_Rbbp8              0.98272300
#> 651_LOC113834282       0.76463145
#> 652_Gys1               0.86650073
#> 653_Hnrnpc             0.95365528
#> 654_Vps35              0.80822413
#> 655_Miga2              0.83144533
#> 656_Epb41l2            0.87684327
#> 657_Tob2               1.00000000
#> 658_Lamtor1            0.89691748
#> 659_LOC100759640       1.00000000
#> 660_Epb41l2            1.00000000
#> 661_Rlim               0.70202424
#> 662_Gys1               0.53647630
#> 663_LOC100750437       0.69265742
#> 664_NA                 0.46858598
#> 665_Nbn                0.89691748
#> 666_Tyw3               0.71891081
#> 667_Gas2l1             0.76463145
#> 668_Fus                1.00000000
#> 669_Prpf38b            0.71891081
#> 670_Calu               1.00000000
#> 671_Rras2              0.53647630
#> 672_Prpf4b             0.80519954
#> 673_Nelfa              0.98272300
#> 674_LOC100754077       1.00000000
#> 675_Rbm28              0.97937571
#> 676_Nsfl1c             0.87684327
#> 677_Rnf126             0.70134699
#> 678_Eme1               0.98272300
#> 679_Nbn                0.11563729
#> 680_Eif4ebp2           0.53647630
#> 681_Wee1               0.80741720
#> 682_Prpf38b            0.89691748
#> 683_Luzp1              0.95365528
#> 684_Gas2l1             1.00000000
#> 685_Pdcd11             0.98272300
#> 686_Chaf1b             0.89691748
#> 687_Pycr1              0.85184664
#> 688_Phf8               0.61804675
#> 689_Raver1             0.98272300
#> 690_Dbn1               0.86650073
#> 691_Dut                0.65273244
#> 692_Prpf4b             0.78401202
#> 693_Prpf4b             0.86650073
#> 694_Efs                0.96317292
#> 695_NA                 0.82069736
#> 696_Ppp2r5b            0.98272300
#> 697_Caskin2            0.46858598
#> 698_Arhgef40           0.65273244
#> 699_Zyx                0.89691748
#> 700_Mphosph10          0.53998241
#> 701_LOC113833392       0.89793287
#> 702_Cdc42ep1           0.90249140
#> 703_Snrpa1             0.74897796
#> 704_Ncbp1              0.95822477
#> 705_Gas2l1             1.00000000
#> 706_Gas2l1             0.93504308
#> 707_Bap1               0.89691748
#> 708_LOC100759640       0.96012819
#> 709_Cherp              0.89691748
#> 710_Nbn                0.53647630
#> 711_LOC100759640       1.00000000
#> 712_NA                 0.53647630
#> 713_Eif3b              1.00000000
#> 714_Miga2              0.86650073
#> 715_Prpf4b             0.61127133
#> 716_Dbn1               0.53998241
#> 717_Ppp2r5b            0.66402033
#> 718_Exosc9             0.75851223
#> 719_Eif3b              0.89691748
#> 720_Ripk2              0.81153692
#> 721_Dlg1               0.91740623
#> 722_N4bp1              0.97848956
#> 723_Nudc               0.61127133
#> 724_Znf367             1.00000000
#> 725_Ring1              1.00000000
#> 726_Snrpa1             0.65273244
#> 727_U2surp             0.80822413
#> 728_LOC100764225       0.59107576
#> 729_Cdc42ep1           0.59107576
#> 730_Znf385a            0.61127133
#> 731_Ints1              0.61127133
#> 732_LOC113833392       0.61127133
#> 733_Lrch4              0.83634429
#> 734_Ctdspl2            0.53647630
#> 735_Prpf4b             0.89691748
#> 736_Luzp1              0.84433926
#> 737_Eif3b              0.98272300
#> 738_Ptpn14             0.76552820
#> 739_Rrp1               0.89691748
#> 740_Lrrfip2            0.46858598
#> 741_Nsfl1c             0.89691748
#> 742_Ddx51              0.82801076
#> 743_Prpf38b            0.55406551
#> 744_Eef1b2             0.89691748
#> 745_Znf385a            0.95365528
#> 746_Map9               0.89691748
#> 747_Rflnb              0.89691748
#> 748_NA                 0.96425465
#> 749_C1H12orf45         0.96425465
#> 750_U2surp             0.70134699
#> 751_Caskin2            0.61127133
#> 752_Eri1               0.80519954
#> 753_Gsk3b              0.81456911
#> 754_LOC100766946       0.98272300
#> 755_Cnpy3              0.99123043
#> 756_Hnrnpc             1.00000000
#> 757_Ptpn14             1.00000000
#> 758_Slc7a11            0.96425465
#> 759_Hnrnpc             0.94766133
#> 760_Cdc37l1            0.76463145
#> 761_LOC100768405       0.82069736
#> 762_Rragc              0.82069736
#> 763_LOC113834282       0.71200897
#> 764_Fus                0.87684327
#> 765_Ubxn1              0.94365155
#> 766_Mmut               0.99821395
#> 767_Pdcd11             0.99040703
#> 768_LOC100757535       0.75851223
#> 769_Eif3b              0.53998241
#> 770_Rnf113a            0.96317292
#> 771_Sytl4              0.74825817
#> 772_Tlnrd1             0.90249140
#> 773_H671_1g1131        0.80007140
#> 774_Neurl1             0.86563801
#> 775_Zyx                1.00000000
#> 776_Ctdspl2            0.64380397
#> 777_Chaf1b             0.53647630
#> 778_Rragc              0.65421849
#> 779_Srfbp1             0.46858598
#> 780_Gys1               0.92300444
#> 781_Usp15              0.98272300
#> 782_Arhgef40           0.87684327
#> 783_Gigyf1             0.95365528
#> 784_Minar1             0.84433926
#> 785_Dus2               0.76552820
#> 786_Gatad2b            0.97937571
#> 787_Eif5               0.89691748
#> 788_Epb41l2            1.00000000
#> 789_Arl6ip4            0.87684327
#> 790_Plin4              0.82069736
#> 791_Elf2               0.85184664
#> 792_Plin4              0.54346611
#> 793_Snip1              0.81456911
#> 794_Txlng              0.62865083
#> 795_LOC100769437       0.53647630
#> 796_Caskin2            0.65273244
#> 797_NA                 0.65273244
#> 798_Synm               0.81550950
#> 799_Synm               1.00000000
#> 800_Ube2c              1.00000000
#> 801_Sgtb               0.92872331
#> 802_Prpf4b             0.61127133
#> 803_Epb41l2            0.71378983
#> 804_Mllt1              0.43271982
#> 805_LOC100759640       1.00000000
#> 806_Epb41l2            0.95731924
#> 807_Znf280b            0.87684327
#> 808_Kiaa1143           1.00000000
#> 809_Gas2l1             0.71891081
#> 810_Srp72              0.87684327
#> 811_Tomm22             0.87684327
#> 812_Psip1              0.92267626
#> 813_Arhgef37           0.53647630
#> 814_Bckdk              0.97937571
#> 815_Strip1             0.71378983
#> 816_Usp15              0.78218056
#> 817_Ssr3               0.82973183
#> 818_Strip1             0.89691748
#> 819_Eif3b              0.61804675
#> 820_U2surp             0.81456911
#> 821_Bend3              0.89691748
#> 822_Rps10              1.00000000
#> 823_Rpl23a             0.78605465
#> 824_Nbn                0.70208246
#> 825_Rpap3              0.96317292
#> 826_LOC100759640       0.53998241
#> 827_Ric8a              0.58004623
#> 828_Hsph1              0.84433926
#> 829_LOC100759640       0.96317292
#> 830_LOC100757535       0.91418289
#> 831_Gigyf1             0.71378983
#> 832_Dbn1               1.00000000
#> 833_Snrk               1.00000000
#> 834_Prpf38b            0.99040703
#> 835_LOC100766868       1.00000000
#> 836_LOC100766868       0.53998241
#> 837_Wbp11              0.77944157
#> 838_Rusc2              0.61127133
#> 839_Eif3b              0.98969498
#> 840_Ptpn14             0.96425465
#> 841_Rlim               1.00000000
#> 842_Ints1              0.71378983
#> 843_Chaf1b             0.46858598
#> 844_Dlg1               1.00000000
#> 845_Lamtor1            0.97937571
#> 846_Tab1               0.99040703
#> 847_Dbn1               0.53998241
#> 848_Psip1              0.53998241
#> 849_Dbn1               0.99790083
#> 850_Pabpc1             0.99040703
#> 851_Hnrnpc             0.89691748
#> 852_Emd                0.76463145
#> 853_LOC100764225       0.61127133
#> 854_Nup50              0.76463145
#> 855_Ctcf               0.71891081
#> 856_Raly               0.96317292
#> 857_Bard1              0.87684327
#> 858_Ptpn14             0.87125683
#> 859_LOC100757535       0.43271982
#> 860_Psmd2              1.00000000
#> 861_Junb               0.53998241
#> 862_C1qbp              0.97937571
#> 863_Lrch4              0.98891432
#> 864_CUNH14orf93        0.82069736
#> 865_U2surp             0.66402033
#> 866_Raly               0.98272300
#> 867_LOC100774417       0.80605441
#> 868_Srp72              1.00000000
#> 869_LOC100764225       1.00000000
#> 870_Morf4l2            0.68492351
#> 871_CUNH9orf40         0.71378983
#> 872_Gas2l1             0.57624609
#> 873_Atp5pf             0.87684327
#> 874_Lrrfip2            0.71378983
#> 875_Prpf4b             0.62405509
#> 876_Top2b              0.53998241
#> 877_Mepce              0.84433926
#> 878_Ptpn14             0.91734510
#> 879_Dnajc25            0.46858598
#> 880_Cbx8               0.95731924
#> 881_Synm               0.82203853
#> 882_Def6               1.00000000
#> 883_Gys1               0.96317292
#> 884_Luzp1              0.76463145
#> 885_Synm               0.87184148
#> 886_Snip1              0.65273244
#> 887_Top2b              0.80741720
#> 888_NA                 1.00000000
#> 889_Trim35             1.00000000
#> 890_Znf385a            0.96317292
#> 891_Chaf1b             0.71891081
#> 892_Abcf1              1.00000000
#> 893_Pdcd11             0.95365528
#> 894_Dlg1               0.79131251
#> 895_Dbn1               0.69476960
#> 896_LOC100752363       0.92669903
#> 897_Ppp4r3a            0.59107576
#> 898_Gas2l1             0.76552820
#> 899_Mtmr10             1.00000000
#> 900_Cyld               0.61127133
#> 901_NA                 0.76552820
#> 902_Rnf113a            0.53647630
#> 903_Nelfa              0.76463145
#> 904_Zkscan1            0.81550950
#> 905_Chaf1b             0.65273244
#> 906_Eif3b              1.00000000
#> 907_Top2b              0.65273244
#> 908_Chaf1b             0.98272300
#> 909_Epb41l2            0.99616460
#> 910_C3H11orf58         0.76463145
#> 911_Top2b              0.99040703
#> 912_Wee1               0.99040703
#> 913_Raly               0.98891432
#> 914_H671_1g2680        0.91418289
#> 915_Eef1b2             0.65273244
#> 916_Gas2l1             0.46858598
#> 917_Epb41l2            0.61127133
#> 918_Rpl23a             0.97937571
#> 919_Chmp2b             0.98272300
#> 920_Lrrfip2            0.46858598
#> 921_Aldoa              0.61127133
#> 922_Cby1               0.76552820
#> 923_LOC100759640       0.96317292
#> 924_Rbm28              0.99616460
#> 925_Skiv2l             0.76552820
#> 926_Ints1              0.95731924
#> 927_Ehd1               1.00000000
#> 928_Nr2f6              0.84433926
#> 929_Top2b              0.53998241
#> 930_Lrrfip2            0.96317292
#> 931_Pip4p2             0.92300444
#> 932_Srp72              0.76463145
#> 933_Mtmr9              0.65273244
#> 934_Gigyf1             0.98380907
#> 935_Rbm7               0.87684327
#> 936_LOC100773565       0.74897796
#> 937_Trim35             0.86650073
#> 938_Cbx8               0.81456911
#> 939_Rplp0              1.00000000
#> 940_Aldoa              0.74897796
#> 941_NA                 0.94365155
#> 942_Zyx                0.76463145
#> 943_Psip1              0.95897789
#> 944_Slc7a11            0.96317292
#> 945_Miga2              0.61127133
#> 946_Arhgef6            0.53998241
#> 947_Dlgap4             0.95365528
#> 948_Ampd2              1.00000000
#> 949_Luzp1              0.91740623
#> 950_Camlg              0.81456911
#> 951_Pfkfb3             0.73477355
#> 952_NA                 1.00000000
#> 953_Raly               1.00000000
#> 954_Kiaa1143           0.43271982
#> 955_Bcar1              0.76552820
#> 956_Gatad2b            0.74897796
#> 957_Eif4ebp2           0.98272300
#> 958_Fam76b             0.57736706
#> 959_Camlg              0.82069736
#> 960_LOC100754077       0.31331155
#> 961_NA                 0.43271982
#> 962_Epb41l2            0.99500130
#> 963_Ankrd34a           1.00000000
#> 964_Zc3h15             0.53647630
#> 965_Def6               0.76552820
#> 966_Srsf6              0.69265742
#> 967_H671_4g11480       0.76463145
#> 968_Top2b              0.88019189
#> 969_LOC100769471       0.87684327
#> 970_Raver1             0.81550950
#> 971_Etv3               1.00000000
#> 972_Psd                0.69476960
#> 973_Usp15              0.99500130
#> 974_Nol7               0.89691748
#> 975_Stk38              0.89691748
#> 976_Smim13             0.98272300
#> 977_Etv3               0.80741720
#> 978_Synm               0.96317292
#> 979_Pwp1               0.76463145
#> 980_Fus                0.87684327
#> 981_Junb               0.89691748
#> 982_Phf8               0.66699560
#> 983_Nelfa              0.98272300
#> 984_Prpf4b             0.87435458
#> 985_Abraxas1           0.87684327
#> 986_Prpf4b             1.00000000
#> 987_Raver1             1.00000000
#> 988_Caap1              0.94365155
#> 989_Rpap3              0.89691748
#> 990_Hsph1              0.69265742
#> 991_LOC100750437       1.00000000
#> 992_Mepce              0.75851223
#> 993_Efs                0.95731924
#> 994_Epb41l2            0.92593443
#> 995_Abcf1              0.80822413
#> 996_NA                 0.86319749
#> 997_Eif4ebp2           1.00000000
#> 998_Pfkfb3             0.98272300
#> 999_Hnrnpc             0.89691748
#> 1000_Psmd2             0.99500130
#> 
#> $Exponential$pvc_pattern_summary
#>   -60 15 60 90 120 240
#> p   0  0  0  1   1   0
#> v   0  1  0  2   0   0
#> b   0  0  1  0   0   0
#> t   0  0  0  0   0   0
#> 
#> $Exponential$plots
#> list()
#> 
#> 
#> $Stationary
#> $Stationary$pvc_adj_pvals
#>                  15_vs_neighbors 60_vs_neighbors 90_vs_neighbors
#> 1_Top2b               0.47623282     0.010069485       0.2528074
#> 2_NA                  1.00000000     1.000000000       1.0000000
#> 3_Snip1               1.00000000     1.000000000       0.9914096
#> 4_Tomm34              0.97399727     0.730884695       0.8594514
#> 5_Pus3                0.93408439     0.612667440       1.0000000
#> 6_Ints1               0.90515009     0.860650480       0.9790592
#> 7_Mlh1                0.93408439     0.949911982       0.7816496
#> 8_LOC100750437        1.00000000     1.000000000       0.2528074
#> 9_Pabpc1              0.54502481     1.000000000       0.7774182
#> 10_Top2b              0.99508555     0.848200886       0.7605719
#> 11_Gorasp1            0.53591151     0.212803322       0.8452398
#> 12_Ints1              1.00000000     1.000000000       1.0000000
#> 13_Syvn1              0.60118421     0.457565692       0.8538114
#> 14_Znf280b            0.63289895     0.536829699       0.7656499
#> 15_Mrnip              0.79319446     0.616363428       0.9998830
#> 16_Rragc              0.90515009     0.839109126       0.7012143
#> 17_Gorasp1            0.71724978     0.119293142       0.7305634
#> 18_Tomm34             0.83839993     0.481902036       0.6563396
#> 19_LOC100757430       0.62388094     0.127441307       0.6923074
#> 20_Ubxn1              0.88292911     0.368330199       0.8049503
#> 21_H671_1g1131        0.54502481     0.400510484       0.2456534
#> 22_Luzp1              0.58745270     0.502705746       0.9298273
#> 23_Efs                0.54502481     0.048783400       0.3971096
#> 24_Mta2               0.94427730     0.523907271       0.7685783
#> 25_Nedd1              0.71607132     0.127428892       0.3648974
#> 26_Gigyf1             0.54325567     0.130931101       0.7416436
#> 27_Myh9               0.90515009     0.093865695       0.2456534
#> 28_Caskin2            0.64073199     0.028561245       0.2528074
#> 29_Papolg             0.54325567     0.310570702       0.3648974
#> 30_Tfg                0.65406828     0.834107317       0.6170691
#> 31_Rpl34              0.52376249     0.113354095       0.5515318
#> 32_Mideas             0.90515009     0.769933446       0.6464593
#> 33_Gys1               0.29488064     0.004123866       0.1865649
#> 34_Arhgef6            0.85192356     0.170774257       0.2792817
#> 35_Ctdspl2            1.00000000     1.000000000       1.0000000
#> 36_Ptpn14             1.00000000     1.000000000       0.6451309
#> 37_Raly               1.00000000     0.996871421       0.9965477
#> 38_Znhit3             0.63851869     0.082558457       0.3226632
#> 39_LOC113833392       0.38373964     0.024970798       1.0000000
#> 40_Luc7l3             0.61241133     0.118235739       0.4873530
#> 41_Rplp0              1.00000000     1.000000000       1.0000000
#> 42_Gys1               0.54502481     0.507654023       0.9790592
#> 43_Rpl22l1            0.49554493     0.299407772       0.7888888
#> 44_Eif3b              0.90515009     0.818415544       0.7073961
#> 45_Med26              0.69624464     0.040469441       0.5214596
#> 46_Mepce              0.97194962     0.874141663       0.7728969
#> 47_Pdcd11             0.63993993     0.307212373       0.8613130
#> 48_Twf1               0.53917994     0.217102161       0.2528074
#> 49_LOC100759640       0.73891031     1.000000000       1.0000000
#> 50_Wrnip1             0.74053880     0.329225607       0.7605719
#> 51_Poldip3            0.52376249     0.047911897       0.4629072
#> 52_Ampd2              0.88137715     0.554525769       0.7038252
#> 53_Mea1               0.79677029     0.949911982       0.9905692
#> 54_Dbn1               0.66796513     0.469833614       0.8939102
#> 55_Snip1              0.97194962     0.981399784       0.9368268
#> 56_Srsf6              0.52063115     0.011705607       0.3975070
#> 57_LOC113834282       1.00000000     1.000000000       1.0000000
#> 58_Map9               0.63385599     0.085099598       0.5004048
#> 59_Cdc42ep1           0.66005075     0.761503234       0.9790592
#> 60_Poldip3            1.00000000     1.000000000       1.0000000
#> 61_LOC100764225       1.00000000     1.000000000       1.0000000
#> 62_Epb41l2            1.00000000     1.000000000       1.0000000
#> 63_H671_4g11480       0.76257173     0.958954018       0.9757581
#> 64_Nbn                1.00000000     1.000000000       1.0000000
#> 65_U2surp             1.00000000     1.000000000       1.0000000
#> 66_Gigyf1             0.54325567     0.073928051       0.3779698
#> 67_NA                 0.92364395     0.554712047       0.7121590
#> 68_Luc7l3             1.00000000     1.000000000       1.0000000
#> 69_LOC100752363       0.72613118     0.101906724       0.4771061
#> 70_Ampd2              1.00000000     1.000000000       1.0000000
#> 71_LOC100759640       1.00000000     1.000000000       1.0000000
#> 72_Stam               0.92291089     0.355297643       0.3614870
#> 73_Nsfl1c             0.54502481     0.678847456       0.4980875
#> 74_Pfkfb3             0.93408439     0.208936896       0.2788533
#> 75_Rad23a             0.92291089     0.041040655       0.2528074
#> 76_Elf2               0.62288573     0.009592916       0.2456534
#> 77_Crem               0.90515009     0.872113131       0.9316856
#> 78_Rragc              0.57155299     0.277582354       0.4412340
#> 79_Lrrfip2            0.93408439     0.934269330       0.9998830
#> 80_Zyx                0.90515009     0.383464190       0.6505823
#> 81_Lrrfip2            0.46779026     0.085215420       0.4586578
#> 82_Gatad2b            0.89679055     0.115365078       0.5098348
#> 83_Bcar1              0.92291089     0.709226993       0.9546590
#> 84_Ehd1               0.54325567     0.025047244       0.3226632
#> 85_LOC113834282       0.53917994     0.048783400       0.4046589
#> 86_Tmem230            0.97878789     0.218109162       0.7605719
#> 87_Ncbp1              0.62945080     0.044152701       0.3966787
#> 88_Mllt1              1.00000000     1.000000000       1.0000000
#> 89_Stk17b             0.73942586     0.099296146       0.3756990
#> 90_Dlgap4             0.97878789     0.218109162       0.7605719
#> 91_Papolg             0.84802386     0.093254560       0.2456534
#> 92_Cyld               0.80885508     0.218109162       0.7605719
#> 93_Gigyf1             0.64073199     0.854188318       0.9164056
#> 94_Lrrfip2            1.00000000     1.000000000       1.0000000
#> 95_Lrrfip2            0.97194962     0.100570871       0.3779698
#> 96_Rlim               0.70303031     0.990436120       0.7969119
#> 97_Eif3b              0.71607132     0.980519065       0.7857770
#> 98_Mphosph10          0.64552746     0.375556171       1.0000000
#> 99_Gatad2b            0.39686202     0.161989284       0.7686332
#> 100_Srsf6             0.60118421     0.062796911       0.5474512
#> 101_Zyx               0.88292911     0.386375070       0.7774182
#> 102_Mphosph10         1.00000000     1.000000000       1.0000000
#> 103_Psip1             0.81804469     0.834107317       0.9542486
#> 104_Fbl               0.99508555     0.835850865       0.9316856
#> 105_H671_1g2680       0.82649872     0.818669157       0.4959790
#> 106_Sgtb              0.62288573     0.881643420       0.9905692
#> 107_Gnl3              0.93408439     0.218109162       0.6170691
#> 108_Eif3b             0.63385599     0.492326464       0.9090431
#> 109_Serpinb1          1.00000000     1.000000000       1.0000000
#> 110_N4bp1             1.00000000     1.000000000       0.7991360
#> 111_Snip1             1.00000000     1.000000000       1.0000000
#> 112_Psip1             1.00000000     1.000000000       1.0000000
#> 113_Mlh1              0.90515009     0.473245549       0.2792817
#> 114_Bsg               0.97194962     0.646473799       0.9606716
#> 115_Tnpo1             0.72880676     0.809159384       0.8049503
#> 116_H671_1g2680       1.00000000     1.000000000       1.0000000
#> 117_Cbx8              1.00000000     1.000000000       1.0000000
#> 118_Mideas            1.00000000     1.000000000       1.0000000
#> 119_Mideas            0.97194962     0.118412423       0.3678174
#> 120_Dcun1d3           0.97194962     0.118412423       0.3678174
#> 121_Dlg1              1.00000000     1.000000000       1.0000000
#> 122_Rad23a            0.54325567     0.024200859       0.5683353
#> 123_Srsf6             0.72613118     0.231705463       0.9790592
#> 124_Stx7              1.00000000     1.000000000       1.0000000
#> 125_Pdcd11            0.88292911     0.726453304       0.6138471
#> 126_Kiaa1958          0.90515009     0.799298832       0.9790592
#> 127_Pwp1              0.89639358     0.822406373       0.8487388
#> 128_Txlng             0.71607132     0.790034699       0.9998830
#> 129_Junb              1.00000000     1.000000000       1.0000000
#> 130_LOC100759640      1.00000000     1.000000000       1.0000000
#> 131_Dbn1              0.45454337     0.048653042       0.7305634
#> 132_Top2b             0.62288573     0.461331485       0.8184153
#> 133_Rusc2             0.93408439     0.399784993       0.6138471
#> 134_NA                0.94606430     1.000000000       1.0000000
#> 135_LOC113837251      0.88563466     0.588310352       0.9857687
#> 136_Fam76b            1.00000000     1.000000000       1.0000000
#> 137_Ptpn14            0.75827918     0.098525135       0.3779698
#> 138_Chmp4b            0.54325567     0.028456138       0.8258043
#> 139_Prpf4b            1.00000000     1.000000000       1.0000000
#> 140_Eif3b             1.00000000     1.000000000       1.0000000
#> 141_Nsfl1c            1.00000000     1.000000000       1.0000000
#> 142_Pdlim7            1.00000000     1.000000000       1.0000000
#> 143_Rnf113a           1.00000000     1.000000000       1.0000000
#> 144_Epb41l2           0.54094413     0.222577072       0.9905692
#> 145_Hnrnpc            1.00000000     1.000000000       1.0000000
#> 146_LOC113834282      0.54325567     0.044152701       0.6170691
#> 147_Plekho2           0.88292911     0.322882788       0.9724011
#> 148_Med26             0.93408439     0.936900048       0.6675377
#> 149_Arhgef40          0.90515009     0.796500636       0.7321951
#> 150_NA                1.00000000     1.000000000       0.8566426
#> 151_Phf8              0.46779026     0.010131878       0.2211880
#> 152_Minar1            0.79319446     1.000000000       1.0000000
#> 153_H671_21690        0.88292911     0.553896530       0.5683353
#> 154_Arhgef40          0.51979551     0.018650257       0.4959790
#> 155_Chaf1b            0.99983034     0.589597408       0.5683353
#> 156_Prpf4b            0.62479340     0.429382386       0.8786434
#> 157_Znf367            0.64828195     0.838455031       0.9740321
#> 158_Luzp1             0.70908163     0.804669079       0.6138471
#> 159_LOC113833882      0.57155299     0.234312222       0.7605719
#> 160_Hnrnpc            1.00000000     1.000000000       1.0000000
#> 161_Mepce             0.53591151     0.095752128       0.3678174
#> 162_Ubxn1             0.93408439     0.947257152       0.6260825
#> 163_Mllt1             0.05820201     0.024970798       0.6749909
#> 164_Chaf1b            0.54094413     0.268296359       0.9905692
#> 165_Raly              0.67549601     0.183136137       0.6669812
#> 166_Gas2l1            0.66502052     0.368796409       0.7605719
#> 167_Dlg1              0.88292911     0.804749173       0.8303413
#> 168_Hoxc10            1.00000000     1.000000000       0.6995470
#> 169_Gigyf1            0.67190462     1.000000000       1.0000000
#> 170_Luzp1             0.71607132     0.495954955       0.7660779
#> 171_Srp72             0.64073199     0.152487189       0.7658717
#> 172_LOC100771461      0.54325567     0.342023429       1.0000000
#> 173_Chaf1b            0.71794780     0.568194019       0.9316856
#> 174_C3H11orf58        0.78137419     0.799298832       1.0000000
#> 175_Pdcd11            1.00000000     1.000000000       1.0000000
#> 176_Psip1             0.62288573     0.822406373       0.2528074
#> 177_Prpf4b            0.97194962     0.231298837       0.3614870
#> 178_Rnf113a           0.57155299     0.044152701       0.4980875
#> 179_Irf3              0.94427730     0.752570747       0.5683353
#> 180_Smim13            0.88292911     0.834107317       0.5369183
#> 181_Gnl3              0.94012802     0.589597408       0.4441837
#> 182_Psma5             0.93408439     0.728733616       0.6223072
#> 183_Ptpn14            0.57155299     0.044152701       0.1973920
#> 184_Prpf4b            1.00000000     0.132571176       0.3917706
#> 185_Top2b             1.00000000     1.000000000       1.0000000
#> 186_Prpf38b           0.65406828     0.072705811       0.2456534
#> 187_Epb41l2           0.93408439     0.864692664       0.6995470
#> 188_Eif3b             1.00000000     1.000000000       1.0000000
#> 189_Hnrnpc            0.85192356     0.962032310       0.7423490
#> 190_LOC100758278      0.89679055     0.260257950       0.3913066
#> 191_Prpf4b            0.92291089     0.818669157       0.5734605
#> 192_Caskin2           0.58496599     0.477961599       0.4010473
#> 193_LOC100752363      0.54325567     0.700740104       0.8184153
#> 194_Septin6           0.93408439     0.818415544       0.6532590
#> 195_Max               0.62388094     0.329225607       0.9688836
#> 196_Mid1ip1           0.61241133     0.093652409       0.5515318
#> 197_NA                0.94792345     0.399784993       0.8391381
#> 198_Hsph1             1.00000000     1.000000000       1.0000000
#> 199_Nol7              0.32526276     0.015678224       0.2528074
#> 200_Raly              0.52376249     0.405200601       0.9149144
#> 201_Smim13            0.70908163     0.944038423       0.5382447
#> 202_LOC100757535      0.65406828     0.093254560       0.2535009
#> 203_Net1              0.89679055     0.072788586       0.2456534
#> 204_LOC100754077      1.00000000     1.000000000       1.0000000
#> 205_Snip1             0.46779026     0.249506635       0.7309207
#> 206_Hnrnpc            0.62288573     0.312784861       0.9905692
#> 207_Ldlrap1           0.85192356     0.868830731       0.8613130
#> 208_Luzp1             0.93408439     0.118412423       0.2792817
#> 209_Rpl26             0.93408439     0.847488627       0.3779698
#> 210_Epb41l2           0.89639358     0.262881922       0.6532590
#> 211_Znf367            0.81448384     0.945572261       0.9585265
#> 212_Dlgap4            1.00000000     1.000000000       1.0000000
#> 213_Plekho2           0.76257173     1.000000000       1.0000000
#> 214_Zpr1              0.92291089     0.804069909       0.8109574
#> 215_Dlgap4            1.00000000     1.000000000       1.0000000
#> 216_Def6              0.94427730     0.754718745       0.9740321
#> 217_Eif4ebp2          0.62479340     0.412411042       0.7605719
#> 218_Eef1b2            0.57311992     0.690767762       0.8594514
#> 219_Rad23a            0.76257173     0.616238304       0.6995470
#> 220_Morf4l2           0.72613118     0.818415544       0.9790592
#> 221_Arhgef40          0.72613118     0.447603231       0.8711279
#> 222_NA                0.78223484     0.585388668       0.4144947
#> 223_LOC100773565      0.77895142     0.799298832       0.7305634
#> 224_Dus2              0.78315628     0.088940513       0.6563396
#> 225_Pip4p2            0.45454337     0.141221624       0.7564088
#> 226_Top2b             0.61035300     0.934269330       0.9894765
#> 227_Znf280b           0.72251712     0.306129267       0.9790592
#> 228_Pdcd11            0.89867219     0.368796409       0.7495548
#> 229_Bckdk             0.79319446     0.937388320       1.0000000
#> 230_Arhgef40          1.00000000     1.000000000       1.0000000
#> 231_Mepce             1.00000000     1.000000000       1.0000000
#> 232_Ccnd3             0.98509912     0.548421073       0.6995470
#> 233_Phf8              1.00000000     1.000000000       1.0000000
#> 234_H671_1g2680       0.99983034     0.813074192       0.9790592
#> 235_Ell               0.93408439     0.991723987       0.8594514
#> 236_U2surp            0.52376249     0.422107012       0.8487388
#> 237_Rps10             0.76257173     0.809465395       0.9905692
#> 238_Ctdspl2           0.54502481     0.144987365       0.9316856
#> 239_Top2b             1.00000000     1.000000000       1.0000000
#> 240_Msantd3           0.62945080     0.938325214       0.3614870
#> 241_Fam76b            0.97399727     0.754718745       0.9790592
#> 242_Ppp4r3a           1.00000000     1.000000000       1.0000000
#> 243_Gpatch4           1.00000000     1.000000000       0.7660779
#> 244_Nudc              0.97194962     0.804749173       0.7245342
#> 245_Nol7              0.82667600     0.460414349       0.5528622
#> 246_Plekho2           0.89639358     0.049525675       0.3614870
#> 247_Prpf4b            0.71794780     0.116295046       0.7073961
#> 248_Mta2              0.64073199     0.285687994       0.8487388
#> 249_U2surp            0.62388094     0.082341425       0.9740321
#> 250_Ubxn1             1.00000000     1.000000000       1.0000000
#> 251_Rlim              0.97194962     0.337140363       0.2817729
#> 252_Atat1             0.54094413     0.237143782       0.7321951
#> 253_Ubxn1             0.98602325     0.798374776       0.8757057
#> 254_H671_1g2680       1.00000000     0.619146987       0.5098348
#> 255_eIF2aK2           0.93975237     0.100881143       0.2792817
#> 256_Skiv2l            0.63851869     0.110015956       1.0000000
#> 257_Rpl28             0.54325567     0.645812581       1.0000000
#> 258_LOC100759640      1.00000000     1.000000000       1.0000000
#> 259_Gatad2b           0.97878789     0.834331354       0.8613130
#> 260_NA                0.62945080     0.383464190       0.9998830
#> 261_Gprasp1           0.93975237     0.589927693       0.8883499
#> 262_Luzp1             1.00000000     1.000000000       1.0000000
#> 263_Slc1a5            1.00000000     1.000000000       1.0000000
#> 264_LOC113834282      0.86639744     0.834107317       0.6110255
#> 265_Srsf6             0.65406828     0.914920513       0.7605719
#> 266_Cdc42ep1          0.72880676     0.378206060       0.6923074
#> 267_Net1              1.00000000     1.000000000       0.5048485
#> 268_Caskin2           0.24835879     0.082341425       0.7964405
#> 269_LOC100759640      0.88292911     0.046953528       0.2762895
#> 270_Mideas            0.70908163     0.859444878       0.4959790
#> 271_Luzp1             0.72613118     0.763870844       0.7746328
#> 272_Emd               1.00000000     1.000000000       1.0000000
#> 273_Plpp6             0.53917994     0.044182198       0.7685783
#> 274_LOC100759640      1.00000000     1.000000000       1.0000000
#> 275_Rps7              0.71607132     0.158835844       0.5528622
#> 276_Fkbp1a            0.54094413     0.554525769       0.6138471
#> 277_Gatad2b           0.97878789     0.307629261       0.3913066
#> 278_Znf385a           0.45454337     0.056503461       0.6494172
#> 279_Arhgef6           0.93408439     0.427929351       0.9905692
#> 280_Slirp             0.94276860     0.935348601       0.8580075
#> 281_Skiv2l            0.75621416     0.985476771       0.8315399
#> 282_H671_21690        0.89639358     0.368796409       0.5098348
#> 283_Kat8              0.54325567     0.078777292       0.6923074
#> 284_Nkap              0.82649872     0.352472407       0.8013220
#> 285_Gsk3b             0.88292911     0.114460768       0.2841615
#> 286_Ints1             0.54325567     0.137187744       0.6563396
#> 287_Gas2l1            0.97935078     1.000000000       1.0000000
#> 288_LOC100759640      1.00000000     1.000000000       1.0000000
#> 289_Top2b             0.72880676     0.024200859       0.3022223
#> 290_Kif20b            0.56867904     0.256751492       0.5853013
#> 291_Phf8              1.00000000     1.000000000       1.0000000
#> 292_Snip1             1.00000000     0.164093919       0.2528074
#> 293_Gsk3b             1.00000000     1.000000000       1.0000000
#> 294_Caskin2           1.00000000     1.000000000       1.0000000
#> 295_C3H11orf58        0.53917994     0.251777985       0.4258505
#> 296_Lrch4             0.71724978     0.249506635       0.7788113
#> 297_LOC113834282      0.54502481     0.218109162       0.6140788
#> 298_LOC100750407      1.00000000     1.000000000       1.0000000
#> 299_LOC113833392      1.00000000     1.000000000       1.0000000
#> 300_LOC113833882      0.89410612     0.256944312       0.4771061
#> 301_Ldlrap1           0.54325567     0.145603610       0.3175328
#> 302_Wee1              0.88292911     0.761386308       0.9790592
#> 303_Caap1             1.00000000     1.000000000       1.0000000
#> 304_Eif4ebp2          0.99508555     0.937520812       0.9790592
#> 305_Ripk2             0.97301371     0.447195064       0.3975070
#> 306_Srp72             0.76257173     0.388957672       0.5683353
#> 307_Taok2             0.85518808     0.941668344       0.6604538
#> 308_Nr2f6             0.45454337     0.010428941       0.2528074
#> 309_Arhgef40          0.27742871     0.079178318       0.5098348
#> 310_Gys1              1.00000000     1.000000000       1.0000000
#> 311_Dlg1              1.00000000     1.000000000       1.0000000
#> 312_Vapb              0.79319446     0.839336901       0.5683353
#> 313_LOC100757535      0.53917994     0.010428941       0.2456534
#> 314_Mkrn2             0.90515009     0.913964935       0.7434850
#> 315_Eif3b             0.87609268     0.225209194       0.5515318
#> 316_Isyna1            0.94543672     1.000000000       1.0000000
#> 317_Prpf4b            0.71724978     0.044182198       0.3175328
#> 318_LOC113833882      0.24835879     0.015678224       0.7073961
#> 319_Lrch4             0.56592593     0.181894212       0.8585766
#> 320_Dbn1              0.69475214     0.140667781       0.5528622
#> 321_Abcf1             0.52376249     0.044152701       0.3175328
#> 322_Ints1             0.32526276     0.011705607       0.2528074
#> 323_C3H11orf58        1.00000000     1.000000000       1.0000000
#> 324_Psma5             0.71524930     0.302100917       0.5677689
#> 325_Fundc1            0.98602325     0.778244792       0.7605719
#> 326_Papolg            0.83839993     1.000000000       1.0000000
#> 327_Mideas            0.93408439     0.368194560       0.8038231
#> 328_Ubxn1             0.71633877     0.546856546       0.8210004
#> 329_Synm              0.85192356     0.774335787       0.9790592
#> 330_Arhgef6           1.00000000     1.000000000       1.0000000
#> 331_Ptpn14            1.00000000     0.343423518       0.7073961
#> 332_Pgrmc1            1.00000000     0.343423518       0.7073961
#> 333_Myh9              0.88292911     0.915527348       0.6563396
#> 334_Etv3              1.00000000     1.000000000       1.0000000
#> 335_Ip6k1             1.00000000     1.000000000       1.0000000
#> 336_Luzp1             1.00000000     1.000000000       1.0000000
#> 337_Ptpn14            1.00000000     1.000000000       1.0000000
#> 338_Caskin2           1.00000000     1.000000000       1.0000000
#> 339_Chaf1b            1.00000000     1.000000000       0.5956977
#> 340_Ubxn1             1.00000000     1.000000000       1.0000000
#> 341_Ube2c             1.00000000     1.000000000       1.0000000
#> 342_Gins2             0.70908163     0.902538619       0.9032349
#> 343_Nlgn2             0.62238996     0.676683041       0.9724011
#> 344_Nf2               0.54325567     0.105583422       0.5683353
#> 345_Pip4p2            0.54502481     0.285227044       0.7073961
#> 346_Emd               1.00000000     1.000000000       1.0000000
#> 347_Top2b             0.81864075     0.393617907       0.7660779
#> 348_Trim35            0.93408439     0.285818903       0.2528074
#> 349_NA                0.66796513     0.635527827       0.3226632
#> 350_NA                0.75827918     0.137491563       1.0000000
#> 351_Mideas            1.00000000     1.000000000       1.0000000
#> 352_Gas2l1            1.00000000     0.124564801       0.2841615
#> 353_Ampd2             0.88292911     0.453321923       0.9332729
#> 354_Calu              1.00000000     1.000000000       1.0000000
#> 355_Fam76b            1.00000000     1.000000000       1.0000000
#> 356_Dlg1              0.65300837     0.079545642       0.4629072
#> 357_Srsf6             0.94792345     0.120554949       0.4555868
#> 358_Chaf1b            0.94792345     0.120554949       0.4555868
#> 359_Dbn1              0.66796513     0.337140363       0.3913066
#> 360_Tcf25             0.93975237     0.834107317       0.8452398
#> 361_Psip1             0.12918526     0.009592916       0.4629072
#> 362_Cnpy3             0.60253962     0.224001282       0.7882100
#> 363_LOC100759640      0.61241133     0.453321923       0.3975070
#> 364_Zyx               1.00000000     1.000000000       0.6771960
#> 365_Lrch4             0.62479340     0.912339550       0.6995470
#> 366_Bola1             1.00000000     1.000000000       1.0000000
#> 367_Znf385a           1.00000000     0.640415997       0.6260825
#> 368_Kif20b            0.98096368     0.457525130       0.6260825
#> 369_Ell               0.54502481     0.462802123       0.6433664
#> 370_Ell               0.77895142     0.316343659       0.7611982
#> 371_Srsf6             0.57186421     0.915557817       0.6715757
#> 372_Pwp1              0.88292911     0.033535898       0.2456534
#> 373_Def6              0.97878789     0.090761275       0.3226632
#> 374_Cbx8              0.66930492     0.752570747       0.8480630
#> 375_Ddx51             0.79319446     0.140667781       0.3447812
#> 376_Psip1             0.93408439     0.295813395       0.7415917
#> 377_Arhgef40          0.63454585     0.713568546       0.7858065
#> 378_Raly              0.89639358     0.818415544       0.9790592
#> 379_NA                0.54325567     0.226989403       0.6767307
#> 380_Lrrfip2           0.93408439     0.833342900       0.5034136
#> 381_Gnl3              0.58496599     0.133601735       0.5474512
#> 382_Caskin2           0.71633877     0.130931101       0.5515318
#> 383_Rragc             0.98100329     0.852754600       0.7305634
#> 384_Caskin2           0.97194962     0.950493002       0.9905692
#> 385_Bcar1             0.54094413     1.000000000       1.0000000
#> 386_Homer3            0.65640867     0.009592916       0.2456534
#> 387_Luzp1             0.78913563     0.048653042       0.6715757
#> 388_N4bp1             0.92291089     0.460414349       1.0000000
#> 389_Ppp4r3a           0.71607132     0.265196832       0.6171612
#> 390_H671_1g2680       0.94427730     0.848200886       0.5853013
#> 391_Gnl3              0.72459850     0.089170899       0.6464664
#> 392_Top2b             0.52376249     0.115348038       0.8613130
#> 393_Oser1             1.00000000     1.000000000       1.0000000
#> 394_Snrk              0.97301371     0.949911982       0.9657031
#> 395_Kat8              0.97763126     0.085358109       0.2456534
#> 396_Raver1            0.65300837     0.438379999       0.8613130
#> 397_Pdcd11            0.52063115     0.020193592       0.2841615
#> 398_Rps20             0.65406828     0.800283545       0.6494172
#> 399_Bsg               0.60253962     0.484162075       0.9332729
#> 400_Raly              0.76257173     0.340197749       0.7660779
#> 401_Pdcd2             0.97194962     0.876503074       0.5683353
#> 402_Caskin2           0.73942586     0.494439642       0.1973920
#> 403_LOC100773571      0.93408439     0.174869256       0.6995470
#> 404_Papolg            1.00000000     1.000000000       1.0000000
#> 405_LOC100757535      0.97878789     0.098525135       0.2528074
#> 406_Caap1             0.89688873     0.788507330       0.6178613
#> 407_Psip1             0.65406828     0.796500636       0.9790592
#> 408_Dbn1              0.71607132     0.267648475       0.8993151
#> 409_Mta2              0.46779026     0.086647449       0.8184153
#> 410_Abcf1             0.46779026     1.000000000       1.0000000
#> 411_LOC100754108      0.84802386     0.378690117       0.5515318
#> 412_Slirp             0.16916527     0.009592916       0.7714479
#> 413_Nelfa             0.85691664     0.120554949       0.3614870
#> 414_Aggf1             0.86772434     0.834107317       0.5567138
#> 415_Bap1              0.69625731     0.015227573       0.2456534
#> 416_Luc7l3            0.54325567     0.062795263       0.4017848
#> 417_Rrp1              0.45454337     0.011705607       0.3631106
#> 418_Wrnip1            1.00000000     1.000000000       0.7073961
#> 419_NA                0.62479340     0.103781228       0.3226632
#> 420_Abcf1             0.62288573     0.272113625       0.7888888
#> 421_Cluap1            0.76670578     0.240172943       0.6464593
#> 422_Hnrnpc            0.87608828     0.942618066       0.5643064
#> 423_Ptpn1             0.97194962     0.806223317       0.9790592
#> 424_Myh9              0.52376249     0.035324687       0.1973920
#> 425_U2surp            0.88009809     0.589597408       0.7251536
#> 426_NA                0.63564757     0.067817361       0.3320750
#> 427_Arhgef40          1.00000000     1.000000000       1.0000000
#> 428_Chaf1b            0.88292911     0.312784861       0.5378301
#> 429_Prpf4b            0.54094413     0.326203498       0.6972062
#> 430_Epb41l2           0.53917994     0.072843266       0.6138471
#> 431_Eif3b             1.00000000     1.000000000       1.0000000
#> 432_Isyna1            0.53917994     0.650300774       0.3226632
#> 433_U2surp            0.92291089     0.123733874       0.2514204
#> 434_LOC100765020      0.63385599     0.798306192       0.9032349
#> 435_Arhgef6           0.52376249     0.919930357       0.8566426
#> 436_Ptpn1             0.93279831     0.231705463       0.8786434
#> 437_Prpf4b            0.93490147     0.368796409       0.8639562
#> 438_Rpl35a            0.64073199     0.272113625       0.7920882
#> 439_Prpf4b            1.00000000     1.000000000       1.0000000
#> 440_Zyx               0.93408439     0.957142372       0.8487388
#> 441_Dbn1              0.86772434     0.279979469       0.2456534
#> 442_Chaf1b            1.00000000     1.000000000       0.7685783
#> 443_LOC113834282      0.89431253     0.158835844       0.5853013
#> 444_Gpsm2             0.45454337     0.011705607       0.1351724
#> 445_LOC100757535      1.00000000     1.000000000       1.0000000
#> 446_Cfap410           0.51985473     0.051078092       0.5853013
#> 447_Epb41l2           0.66796513     0.343423518       0.7605719
#> 448_Ncbp1             1.00000000     1.000000000       1.0000000
#> 449_Pacsin1           1.00000000     1.000000000       1.0000000
#> 450_Cstf2             0.61888944     0.093865695       0.7251536
#> 451_LOC100769437      0.76257173     0.981399784       0.9790592
#> 452_eIF2aK2           0.92291089     0.872113131       0.9894765
#> 453_Kiaa1191          0.97194962     0.799298832       0.9998830
#> 454_Mepce             0.54325567     0.678847456       0.8594514
#> 455_Cbx8              0.34788533     0.011705607       0.4629072
#> 456_Eed               0.62288573     0.494439642       0.9790592
#> 457_Cdc42ep1          0.88292911     0.227624172       0.5515318
#> 458_Lrrfip2           1.00000000     1.000000000       1.0000000
#> 459_Pacsin1           0.61241133     0.536829699       0.9790592
#> 460_Gpatch4           0.73967958     0.167212385       0.8303413
#> 461_Plin4             1.00000000     1.000000000       1.0000000
#> 462_NA                0.66687948     0.044182198       0.6715757
#> 463_Snip1             0.94543672     0.329225607       0.9546590
#> 464_Cyld              0.97194962     0.312011245       0.8184153
#> 465_Plin4             0.99508555     0.796160277       0.7895709
#> 466_Twf1              0.86772434     1.000000000       1.0000000
#> 467_LOC113834282      0.86772434     1.000000000       1.0000000
#> 468_Snip1             0.72880676     0.521352863       0.7580342
#> 469_Ppp4r3a           0.90515009     0.833342900       0.9790592
#> 470_Psip1             0.97399727     0.455298716       0.6747728
#> 471_Dnajc5            1.00000000     1.000000000       1.0000000
#> 472_Phf8              0.93408439     0.738183804       0.9998830
#> 473_Bola1             0.53591151     0.269847805       0.9047093
#> 474_Cdc42ep1          1.00000000     0.839266628       0.8690306
#> 475_Eif4ebp2          0.63385599     0.282421060       0.9201656
#> 476_Prpf38b           0.99508555     0.330329679       0.4711825
#> 477_Klhl26            0.83989471     0.065943776       0.2456534
#> 478_Hsph1             1.00000000     1.000000000       1.0000000
#> 479_Snip1             0.97194962     0.938325214       0.9164056
#> 480_Caskin2           0.54325567     0.065565791       0.8480630
#> 481_Plpp6             0.57311992     0.249506635       0.7857770
#> 482_NA                0.61241133     0.061202699       0.2456534
#> 483_Mlh1              0.79604118     0.181894212       0.7038252
#> 484_Gys1              0.57155299     0.009592916       0.2456534
#> 485_Tfg               0.16916527     0.009592916       0.5477890
#> 486_Arhgef6           1.00000000     1.000000000       1.0000000
#> 487_Mphosph10         1.00000000     1.000000000       1.0000000
#> 488_Hoxc10            0.92291089     0.322882788       0.8476647
#> 489_LOC100759640      1.00000000     1.000000000       1.0000000
#> 490_Arhgef40          0.93975237     0.355297643       0.6260825
#> 491_Dnajc5            1.00000000     1.000000000       1.0000000
#> 492_Tbc1d23           1.00000000     1.000000000       1.0000000
#> 493_Ubxn1             0.99251976     0.864692664       0.9047093
#> 494_Rab1a             1.00000000     1.000000000       1.0000000
#> 495_Eif3b             1.00000000     1.000000000       1.0000000
#> 496_Tceal8            0.57155299     0.093254560       0.5515318
#> 497_Dlgap4            0.76257173     0.457565692       0.9316856
#> 498_Smim13            0.92364395     0.853078223       0.7605719
#> 499_NA                1.00000000     1.000000000       1.0000000
#> 500_Lrch4             1.00000000     1.000000000       1.0000000
#> 501_Bola1             1.00000000     1.000000000       1.0000000
#> 502_NA                0.87609268     0.452259090       0.5382447
#> 503_Ptpn14            0.84177776     0.492332140       0.7434850
#> 504_LOC100759640      0.78328468     0.606386718       0.7855578
#> 505_Rps10             0.98509912     0.713568546       0.9509986
#> 506_Top2b             1.00000000     1.000000000       1.0000000
#> 507_Ssr3              0.86774666     0.554812277       0.8181409
#> 508_Homer3            1.00000000     1.000000000       1.0000000
#> 509_Phf8              0.58745270     0.492326464       0.6592642
#> 510_LOC100767716      0.58496599     0.046804463       0.5055045
#> 511_Xpa               0.86863314     0.962032310       0.8711279
#> 512_H671_21690        0.54502481     0.108920594       0.3975070
#> 513_LOC100769471      0.90515009     0.487247401       0.8303413
#> 514_Gas2l1            0.54325567     0.359604268       0.8487388
#> 515_Luzp1             0.71607132     0.043582798       0.6260825
#> 516_Gpbp1             1.00000000     0.285818903       0.9905692
#> 517_Gatad2b           1.00000000     1.000000000       1.0000000
#> 518_Gys1              0.70880247     0.598849644       0.5435710
#> 519_Top2b             0.72880676     0.251648545       0.4164747
#> 520_LOC100757535      1.00000000     0.329225607       0.5535905
#> 521_Lpcat4            0.88292911     0.118235739       0.5341971
#> 522_Arhgef6           1.00000000     1.000000000       1.0000000
#> 523_Cavin3            0.93975237     0.378690117       0.7605719
#> 524_Gpatch4           0.54502481     0.915557817       0.4960784
#> 525_Prpf38b           0.60502460     0.178156017       0.5683353
#> 526_Timm8a            0.64073199     0.254215304       0.4959790
#> 527_Cavin3            0.52376249     0.036591030       0.3922119
#> 528_Mkrn2             0.64828195     0.161989284       0.5098348
#> 529_Oser1             0.78999981     0.211749285       0.6171612
#> 530_Gsk3b             0.64073199     0.254215304       0.4959790
#> 531_Eef1b2            1.00000000     1.000000000       1.0000000
#> 532_Ampd2             0.54094413     0.251648545       0.6532590
#> 533_Lrrfip2           0.47623282     0.598849644       0.8402942
#> 534_Ring1             0.59557209     0.018987223       0.2528074
#> 535_Rlim              0.47623282     0.027421699       0.2528074
#> 536_LOC100759640      0.89679055     0.101906724       0.5546358
#> 537_LOC100759640      0.99508555     0.441780823       0.4428999
#> 538_Atp5pf            0.69624464     0.390006870       0.7321951
#> 539_Max               0.79319446     0.165335918       0.9164056
#> 540_Bap1              0.78328468     1.000000000       1.0000000
#> 541_Nsfl1c            0.54820576     0.536829699       0.4959790
#> 542_Prpf4b            1.00000000     0.936428449       0.8184153
#> 543_LOC100757535      0.93408439     0.712149728       0.6068931
#> 544_Mtmr10            0.64828195     0.354124539       0.4010473
#> 545_Hoxc10            0.47504094     0.095752128       0.7658717
#> 546_Trim35            1.00000000     1.000000000       1.0000000
#> 547_Eif4ebp2          0.77929811     0.145603610       0.5683353
#> 548_Dlgap4            0.72880676     0.507654023       1.0000000
#> 549_Gys1              0.54502481     0.318408177       0.8594514
#> 550_Sgtb              0.71607132     0.223869493       0.9332729
#> 551_Eri2              1.00000000     1.000000000       0.8639562
#> 552_Ccnd3             0.93408439     1.000000000       1.0000000
#> 553_Smim13            0.61241133     0.825474804       0.9047093
#> 554_Snrk              0.54325567     0.195989817       0.4555868
#> 555_Caskin2           1.00000000     1.000000000       1.0000000
#> 556_Pdcd11            0.47504094     0.095752128       0.7658717
#> 557_Pgam5             0.88292911     0.319757459       0.6995470
#> 558_Mphosph10         1.00000000     1.000000000       1.0000000
#> 559_Mideas            1.00000000     1.000000000       1.0000000
#> 560_Top2b             0.60118421     0.190393736       0.9790592
#> 561_LOC100763014      0.92291089     0.181894212       0.3397790
#> 562_Snip1             1.00000000     1.000000000       1.0000000
#> 563_Ubxn1             0.88292911     0.755595789       0.9740321
#> 564_LOC100750407      0.60580980     0.015227573       0.3141525
#> 565_Morf4l2           1.00000000     1.000000000       1.0000000
#> 566_Ctdspl2           1.00000000     1.000000000       1.0000000
#> 567_Cwf19l1           1.00000000     1.000000000       1.0000000
#> 568_Eef1b2            1.00000000     1.000000000       1.0000000
#> 569_C1H12orf45        1.00000000     1.000000000       1.0000000
#> 570_Znf367            1.00000000     0.589597408       0.7154080
#> 571_Ankrd34a          0.54325567     0.118235739       0.7038252
#> 572_Mllt11            0.72613118     0.921066131       0.7858065
#> 573_LOC100774792      0.67424968     0.399639137       0.7309207
#> 574_NA                0.94427730     0.108920594       0.2515247
#> 575_Cbx8              0.70908163     0.598849644       0.6138471
#> 576_Bckdk             0.89639358     0.452259090       0.9790592
#> 577_Snip1             0.46033394     1.000000000       1.0000000
#> 578_Nsfl1c            0.97194962     0.839336901       0.8356342
#> 579_Gas2l1            1.00000000     1.000000000       1.0000000
#> 580_Nudc              0.90515009     0.936428449       0.9790592
#> 581_Epb41l2           0.54094413     0.058615300       0.2456534
#> 582_Mtmr6             0.70731481     1.000000000       1.0000000
#> 583_Znf668            0.53591151     0.047572498       0.6256015
#> 584_Hsph1             1.00000000     1.000000000       0.9905692
#> 585_LOC113834282      0.71607132     0.770566581       0.8594514
#> 586_Ctdspl2           0.54325567     0.072843266       0.7605719
#> 587_Foxf1             0.87608828     0.351779941       0.8480630
#> 588_Luzp1             0.99508555     0.272113625       0.7321951
#> 589_Xpa               0.63993993     0.799298832       1.0000000
#> 590_Psip1             0.51979551     0.030851148       0.5721803
#> 591_Rbm7              0.65406828     0.062795263       0.5341971
#> 592_Mtrex             0.97194962     0.625117006       0.6923074
#> 593_Arhgef40          0.54502481     0.231705463       0.6715757
#> 594_Plekho2           0.89679055     0.803553113       0.5965493
#> 595_Bckdk             0.64073199     0.660718983       0.9870435
#> 596_Dut               1.00000000     1.000000000       1.0000000
#> 597_Abcf1             1.00000000     1.000000000       1.0000000
#> 598_Txnl1             1.00000000     1.000000000       1.0000000
#> 599_Nudc              1.00000000     1.000000000       1.0000000
#> 600_Sh3gl1            0.47623282     0.028743673       0.3614870
#> 601_Gatad2b           0.60253962     0.279979469       0.6171612
#> 602_Homer3            0.66796513     0.438379999       0.8391381
#> 603_Septin6           0.47623282     0.100570871       0.6835401
#> 604_Smim13            0.99176908     0.834107317       0.6464593
#> 605_Arhgef40          0.95686779     0.880579860       0.7416436
#> 606_Rpl32             0.97878789     0.209461038       0.5956977
#> 607_Tomm34            0.71607132     0.359986277       0.8303413
#> 608_Mlh1              0.92364395     0.756572162       0.9790592
#> 609_Tbcc              0.71724978     0.990436120       0.5474512
#> 610_Eif3d             0.54325567     0.809465395       0.5148762
#> 611_Snrk              0.62288573     0.685883250       0.7858065
#> 612_Bckdk             1.00000000     1.000000000       1.0000000
#> 613_Wdr3              0.86427681     0.388957672       0.6464664
#> 614_LOC100757535      1.00000000     1.000000000       1.0000000
#> 615_Dlg1              0.86772434     0.532690647       0.5956977
#> 616_LOC100767716      0.62479340     0.322882788       0.6995470
#> 617_Hnrnpc            0.54502481     0.452259090       0.9998830
#> 618_Mphosph10         0.05820201     0.103875602       0.8594514
#> 619_Eif3b             0.93490147     0.707124470       0.8367759
#> 620_Emd               1.00000000     1.000000000       1.0000000
#> 621_Txlng             0.71607132     0.860650480       0.6138471
#> 622_Prpf4b            1.00000000     1.000000000       1.0000000
#> 623_Rlim              0.79319446     1.000000000       1.0000000
#> 624_Eef1b2            1.00000000     1.000000000       0.9998830
#> 625_Def6              0.52376249     0.364712626       0.9869263
#> 626_LOC100765020      0.64552746     0.623579485       0.9998830
#> 627_U2surp            0.75854480     0.995949366       0.7901560
#> 628_Elf2              0.54094413     0.011705607       0.2456534
#> 629_Slc1a5            0.70908163     0.124564801       0.5098348
#> 630_NA                0.45904970     0.367092213       0.9058986
#> 631_Tfg               1.00000000     1.000000000       1.0000000
#> 632_Top2b             0.93408439     0.254215304       0.3756990
#> 633_Pip4p2            0.99176908     0.585388668       0.4892722
#> 634_Cdc42ep1          1.00000000     1.000000000       1.0000000
#> 635_Hsph1             1.00000000     1.000000000       1.0000000
#> 636_Twf1              1.00000000     1.000000000       1.0000000
#> 637_Nbn               0.71466746     0.317281054       0.6464664
#> 638_Psmd4             0.68196409     0.756572162       0.9795073
#> 639_Bap1              0.89639358     0.188600072       0.6860277
#> 640_Mepce             0.93408439     0.368796409       0.4398391
#> 641_Mideas            0.89047058     0.093583162       0.2456534
#> 642_LOC100759640      0.54325567     0.081632755       0.5214596
#> 643_Epb41l2           1.00000000     1.000000000       1.0000000
#> 644_Sav1              0.75827918     0.926824194       0.7888888
#> 645_Prpf4b            0.84802386     0.219937022       0.7305634
#> 646_Gnas              0.24835879     0.048653042       0.6110255
#> 647_Mllt1             1.00000000     1.000000000       1.0000000
#> 648_Poldip3           0.47623282     0.044152701       0.3975070
#> 649_Aldoa             0.14642828     0.009592916       0.2762895
#> 650_Rbbp8             0.46779026     0.018987223       0.2788533
#> 651_LOC113834282      0.71607132     0.072705811       0.2528074
#> 652_Gys1              0.94427730     0.752570747       0.4959790
#> 653_Hnrnpc            0.54502481     0.099296146       0.7423490
#> 654_Vps35             0.90515009     0.612667440       0.8594514
#> 655_Miga2             0.82187420     0.339218607       0.6093150
#> 656_Epb41l2           0.85192356     0.329225607       0.7321951
#> 657_Tob2              1.00000000     1.000000000       1.0000000
#> 658_Lamtor1           1.00000000     1.000000000       1.0000000
#> 659_LOC100759640      1.00000000     1.000000000       1.0000000
#> 660_Epb41l2           1.00000000     0.023674400       0.3175328
#> 661_Rlim              0.49554493     0.011705607       0.2657991
#> 662_Gys1              0.26292836     1.000000000       1.0000000
#> 663_LOC100750437      0.77895142     0.224001282       0.2762895
#> 664_NA                0.99983034     0.170774257       0.9998830
#> 665_Nbn               0.98602325     0.058324968       0.4555868
#> 666_Tyw3              1.00000000     1.000000000       1.0000000
#> 667_Gas2l1            0.92291089     0.752570747       0.8391381
#> 668_Fus               1.00000000     1.000000000       1.0000000
#> 669_Prpf38b           0.53917994     0.487247401       0.5683353
#> 670_Calu              1.00000000     1.000000000       0.6138471
#> 671_Rras2             0.54502481     0.272113625       0.9130044
#> 672_Prpf4b            0.81804469     0.860650480       0.9661887
#> 673_Nelfa             0.54325567     0.009592916       0.1865649
#> 674_LOC100754077      0.54094413     0.009592916       0.1865649
#> 675_Rbm28             0.98100329     0.874141663       0.6464593
#> 676_Nsfl1c            0.66930492     0.030810786       0.3141525
#> 677_Rnf126            0.89047058     0.161989284       0.6260825
#> 678_Eme1              0.64073199     0.616363428       0.5956977
#> 679_Nbn               0.82932424     0.044182198       0.2456534
#> 680_Eif4ebp2          0.62288573     0.044152701       0.2792817
#> 681_Wee1              1.00000000     1.000000000       1.0000000
#> 682_Prpf38b           1.00000000     0.569001074       0.9407237
#> 683_Luzp1             0.60118421     0.251648545       0.8480630
#> 684_Gas2l1            1.00000000     1.000000000       1.0000000
#> 685_Pdcd11            1.00000000     1.000000000       1.0000000
#> 686_Chaf1b            0.66796513     0.041040655       0.4306813
#> 687_Pycr1             0.47623282     0.033239223       0.4586578
#> 688_Phf8              1.00000000     0.021259362       0.1973920
#> 689_Raver1            1.00000000     1.000000000       1.0000000
#> 690_Dbn1              0.62388094     0.015227573       0.4959790
#> 691_Dut               0.46779026     0.078754217       0.3397691
#> 692_Prpf4b            0.97878789     0.589597408       1.0000000
#> 693_Prpf4b            0.90515009     0.885437511       0.9672425
#> 694_Efs               0.99508555     0.219937022       0.6532590
#> 695_NA                0.82561145     0.113354095       0.3071868
#> 696_Ppp2r5b           1.00000000     1.000000000       1.0000000
#> 697_Caskin2           0.99508555     0.727132902       0.9790592
#> 698_Arhgef40          0.72399299     0.769933446       0.9998830
#> 699_Zyx               1.00000000     1.000000000       1.0000000
#> 700_Mphosph10         0.64828195     0.532733162       0.6532590
#> 701_LOC113833392      0.49486143     0.080714181       0.5382447
#> 702_Cdc42ep1          0.52376249     0.004123866       0.1351724
#> 703_Snrpa1            0.67549601     1.000000000       1.0000000
#> 704_Ncbp1             0.88292911     0.990436120       0.9790592
#> 705_Gas2l1            1.00000000     1.000000000       1.0000000
#> 706_Gas2l1            0.14642828     0.004123866       0.2456534
#> 707_Bap1              0.94276860     0.343423518       0.6464664
#> 708_LOC100759640      0.64828195     0.100570871       0.6068931
#> 709_Cherp             1.00000000     1.000000000       1.0000000
#> 710_Nbn               0.64828195     0.080948689       0.2792817
#> 711_LOC100759640      1.00000000     1.000000000       1.0000000
#> 712_NA                1.00000000     1.000000000       1.0000000
#> 713_Eif3b             1.00000000     1.000000000       1.0000000
#> 714_Miga2             0.87581193     0.616363428       0.6563396
#> 715_Prpf4b            0.52376249     0.616363428       0.8939102
#> 716_Dbn1              0.71724978     0.635079569       0.9090431
#> 717_Ppp2r5b           0.47623282     0.007058261       0.2456534
#> 718_Exosc9            1.00000000     1.000000000       1.0000000
#> 719_Eif3b             0.62388094     0.114460768       0.9790592
#> 720_Ripk2             1.00000000     1.000000000       1.0000000
#> 721_Dlg1              0.65406828     0.158394831       0.6464664
#> 722_N4bp1             0.83654835     0.319858659       0.8487388
#> 723_Nudc              0.94792345     0.100570871       0.4010473
#> 724_Znf367            0.62288573     0.015227573       0.6464664
#> 725_Ring1             0.62388094     0.839109126       0.5214596
#> 726_Snrpa1            0.62945080     0.737458633       0.4714871
#> 727_U2surp            0.45454337     0.080268242       0.4555868
#> 728_LOC100764225      0.29488064     0.018737895       0.5683353
#> 729_Cdc42ep1          1.00000000     1.000000000       1.0000000
#> 730_Znf385a           1.00000000     1.000000000       1.0000000
#> 731_Ints1             1.00000000     1.000000000       1.0000000
#> 732_LOC113833392      1.00000000     1.000000000       1.0000000
#> 733_Lrch4             0.86639744     0.756572162       0.9796476
#> 734_Ctdspl2           0.54094413     0.268296359       0.7027710
#> 735_Prpf4b            0.97399727     0.834331354       0.8480630
#> 736_Luzp1             1.00000000     1.000000000       1.0000000
#> 737_Eif3b             1.00000000     1.000000000       1.0000000
#> 738_Ptpn14            1.00000000     1.000000000       1.0000000
#> 739_Rrp1              0.69625731     0.319757459       0.7305634
#> 740_Lrrfip2           0.66796513     1.000000000       1.0000000
#> 741_Nsfl1c            1.00000000     1.000000000       1.0000000
#> 742_Ddx51             0.97194962     0.291608060       0.7245342
#> 743_Prpf38b           0.71607132     0.589597408       0.8184153
#> 744_Eef1b2            0.97878789     0.885437511       0.5535905
#> 745_Znf385a           1.00000000     1.000000000       1.0000000
#> 746_Map9              0.53917994     0.081632755       0.8129386
#> 747_Rflnb             0.53917994     0.081632755       0.8129386
#> 748_NA                0.65406828     0.588411064       0.3975070
#> 749_C1H12orf45        0.65406828     0.588411064       0.3975070
#> 750_U2surp            0.82667600     0.180007676       0.4629072
#> 751_Caskin2           0.47623282     0.009592916       0.2528074
#> 752_Eri1              0.71633877     0.874141663       0.7044591
#> 753_Gsk3b             0.72880676     0.848200886       0.7660779
#> 754_LOC100766946      0.52376249     0.028456138       0.6530086
#> 755_Cnpy3             1.00000000     1.000000000       0.5214596
#> 756_Hnrnpc            1.00000000     1.000000000       0.7321951
#> 757_Ptpn14            1.00000000     1.000000000       0.7321951
#> 758_Slc7a11           0.63385599     0.018987223       0.2528074
#> 759_Hnrnpc            0.71633877     0.065157196       0.3648974
#> 760_Cdc37l1           0.97194962     0.663905553       0.8303413
#> 761_LOC100768405      0.93408439     0.554812277       0.7305634
#> 762_Rragc             0.93408439     0.554812277       0.7305634
#> 763_LOC113834282      0.46779026     0.030685499       0.3779698
#> 764_Fus               0.62288573     0.799298832       0.6171612
#> 765_Ubxn1             0.45454337     0.011705607       0.3975070
#> 766_Mmut              0.93408439     0.209461038       0.4144947
#> 767_Pdcd11            0.88924209     0.499536019       0.5554079
#> 768_LOC100757535      1.00000000     1.000000000       1.0000000
#> 769_Eif3b             0.84802386     0.721771211       0.9790592
#> 770_Rnf113a           1.00000000     1.000000000       1.0000000
#> 771_Sytl4             1.00000000     1.000000000       1.0000000
#> 772_Tlnrd1            0.54502481     0.417597821       0.7321951
#> 773_H671_1g1131       0.92496501     0.809738231       0.9657031
#> 774_Neurl1            0.88292911     0.838455031       0.9755035
#> 775_Zyx               0.52376249     0.108821927       0.2528074
#> 776_Ctdspl2           0.70908163     0.158835844       0.7201570
#> 777_Chaf1b            0.89867219     0.350582445       0.5148762
#> 778_Rragc             1.00000000     1.000000000       1.0000000
#> 779_Srfbp1            0.65406828     0.815429174       0.9316856
#> 780_Gys1              1.00000000     0.822406373       0.6974470
#> 781_Usp15             0.93408439     0.901620154       0.8452398
#> 782_Arhgef40          0.97935078     0.452259090       0.4873530
#> 783_Gigyf1            1.00000000     0.368796409       0.7787875
#> 784_Minar1            0.57451925     0.124140958       0.2792817
#> 785_Dus2              1.00000000     1.000000000       1.0000000
#> 786_Gatad2b           0.72613118     0.813074192       0.7564088
#> 787_Eif5              0.94427730     0.985476771       0.9740321
#> 788_Epb41l2           0.57448002     0.255538068       0.9790592
#> 789_Arl6ip4           0.87608828     0.721593970       0.9998830
#> 790_Plin4             0.88292911     0.368796409       0.9918682
#> 791_Elf2              0.99275666     0.401225862       0.5535905
#> 792_Plin4             0.97399727     0.342023429       1.0000000
#> 793_Snip1             1.00000000     1.000000000       1.0000000
#> 794_Txlng             1.00000000     1.000000000       1.0000000
#> 795_LOC100769437      1.00000000     0.996871421       0.9788357
#> 796_Caskin2           0.67729723     0.322882788       0.6434996
#> 797_NA                0.24835879     0.025997938       0.3346013
#> 798_Synm              1.00000000     1.000000000       1.0000000
#> 799_Synm              0.97194962     0.295813395       0.8184153
#> 800_Ube2c             1.00000000     0.609625186       0.9905692
#> 801_Sgtb              0.66796513     0.082558457       0.5430705
#> 802_Prpf4b            0.58745270     0.797826368       0.9585265
#> 803_Epb41l2           0.14642828     0.161679006       0.9740321
#> 804_Mllt1             0.84802386     0.834107317       0.7776935
#> 805_LOC100759640      1.00000000     1.000000000       1.0000000
#> 806_Epb41l2           0.63851869     0.633144523       0.8786434
#> 807_Znf280b           1.00000000     0.388957672       0.8452398
#> 808_Kiaa1143          1.00000000     1.000000000       1.0000000
#> 809_Gas2l1            1.00000000     1.000000000       1.0000000
#> 810_Srp72             0.64828195     0.359667514       0.8939102
#> 811_Tomm22            0.93408439     0.422286817       0.4549154
#> 812_Psip1             0.93975237     0.071271034       0.3036211
#> 813_Arhgef37          0.92291089     0.632610498       0.9905692
#> 814_Bckdk             0.93975237     0.804749173       0.9449472
#> 815_Strip1            0.47504094     0.009592916       0.2762895
#> 816_Usp15             1.00000000     1.000000000       1.0000000
#> 817_Ssr3              0.47623282     0.616363428       0.5341971
#> 818_Strip1            0.54094413     0.141221624       0.5515318
#> 819_Eif3b             0.32526276     0.008305263       0.2514204
#> 820_U2surp            0.93408439     0.076614751       0.2528074
#> 821_Bend3             0.64073199     0.018987223       0.3320750
#> 822_Rps10             0.82649872     0.634193211       0.5956977
#> 823_Rpl23a            0.71633877     0.700740104       0.7400550
#> 824_Nbn               0.97194962     0.848200886       0.5515318
#> 825_Rpap3             0.72613118     0.316008364       0.6563396
#> 826_LOC100759640      1.00000000     1.000000000       1.0000000
#> 827_Ric8a             0.94792345     0.690767762       0.7073961
#> 828_Hsph1             0.67729723     0.551322191       0.8487388
#> 829_LOC100759640      0.66796513     0.137187744       0.3175328
#> 830_LOC100757535      0.89639358     0.838455031       0.7656499
#> 831_Gigyf1            0.94427730     0.195989817       0.9790592
#> 832_Dbn1              0.71724978     0.105583422       0.6260825
#> 833_Snrk              0.88292911     0.119427976       0.6715757
#> 834_Prpf38b           0.54094413     0.246083797       0.9728359
#> 835_LOC100766868      0.05820201     0.004123866       0.6171612
#> 836_LOC100766868      0.87530023     0.492052392       0.8993151
#> 837_Wbp11             0.93408439     0.752570747       0.7888888
#> 838_Rusc2             0.90515009     0.044182198       0.1316251
#> 839_Eif3b             0.62288573     0.378206060       0.4909478
#> 840_Ptpn14            0.72880676     0.492052392       0.9606716
#> 841_Rlim              0.83839993     0.941385992       0.9332729
#> 842_Ints1             0.64828195     1.000000000       1.0000000
#> 843_Chaf1b            0.54502481     0.598849644       0.9790592
#> 844_Dlg1              1.00000000     1.000000000       1.0000000
#> 845_Lamtor1           0.56985618     0.028561245       0.2528074
#> 846_Tab1              0.65194228     0.044152701       0.3226632
#> 847_Dbn1              1.00000000     1.000000000       1.0000000
#> 848_Psip1             1.00000000     1.000000000       1.0000000
#> 849_Dbn1              0.66796513     0.541728497       0.9332729
#> 850_Pabpc1            0.76257173     0.339218607       0.7714479
#> 851_Hnrnpc            0.91576521     0.158835844       0.2528074
#> 852_Emd               0.54325567     0.025997938       0.1865649
#> 853_LOC100764225      0.54325567     0.083810939       0.4412340
#> 854_Nup50             0.90515009     0.322882788       0.6013451
#> 855_Ctcf              0.71724978     0.124065223       0.6138471
#> 856_Raly              0.94792345     0.246900570       0.3966787
#> 857_Bard1             0.45454337     0.022880728       0.4586578
#> 858_Ptpn14            0.86863314     0.460414349       0.8258043
#> 859_LOC100757535      1.00000000     0.279979469       0.8049503
#> 860_Psmd2             0.88292911     0.640060350       1.0000000
#> 861_Junb              0.97399727     0.532733162       0.5515318
#> 862_C1qbp             0.99983034     0.804069909       0.9790592
#> 863_Lrch4             0.97878789     0.045386395       0.2456534
#> 864_CUNH14orf93       0.61241133     0.591600668       1.0000000
#> 865_U2surp            0.76257173     0.788096638       0.5535905
#> 866_Raly              1.00000000     1.000000000       1.0000000
#> 867_LOC100774417      0.93408439     0.187920561       0.3226632
#> 868_Srp72             0.98602325     0.860650480       0.8402942
#> 869_LOC100764225      1.00000000     1.000000000       1.0000000
#> 870_Morf4l2           0.88137715     0.809465395       0.6485884
#> 871_CUNH9orf40        0.66005075     0.079545642       0.4808161
#> 872_Gas2l1            0.16916527     0.004123866       0.1316251
#> 873_Atp5pf            0.72251712     0.099723851       0.2841615
#> 874_Lrrfip2           0.64073199     0.116295046       0.5683353
#> 875_Prpf4b            1.00000000     1.000000000       1.0000000
#> 876_Top2b             1.00000000     1.000000000       1.0000000
#> 877_Mepce             0.93975237     0.958954018       0.8711279
#> 878_Ptpn14            0.72880676     0.145603610       0.9332729
#> 879_Dnajc25           0.97194962     0.998833004       0.7743909
#> 880_Cbx8              0.62288573     0.373980374       0.9998830
#> 881_Synm              0.65406828     0.220944799       0.6171612
#> 882_Def6              0.53793240     0.241416749       0.6260825
#> 883_Gys1              0.89639358     0.462476170       0.8303413
#> 884_Luzp1             0.99508555     0.532733162       0.7012143
#> 885_Synm              0.79677029     0.566071310       0.9740321
#> 886_Snip1             1.00000000     1.000000000       1.0000000
#> 887_Top2b             0.47623282     0.005603759       0.1351724
#> 888_NA                0.54094413     0.113354095       0.8594514
#> 889_Trim35            0.72613118     0.981399784       0.8391381
#> 890_Znf385a           0.54325567     0.093254560       0.7337186
#> 891_Chaf1b            0.93408439     0.756572162       0.5683353
#> 892_Abcf1             1.00000000     0.269441152       0.6260825
#> 893_Pdcd11            1.00000000     1.000000000       1.0000000
#> 894_Dlg1              0.70908163     0.634718999       0.7321951
#> 895_Dbn1              1.00000000     1.000000000       1.0000000
#> 896_LOC100752363      0.79677029     0.113354095       0.5341971
#> 897_Ppp4r3a           1.00000000     1.000000000       1.0000000
#> 898_Gas2l1            0.97194962     0.227624172       0.4555868
#> 899_Mtmr10            1.00000000     1.000000000       1.0000000
#> 900_Cyld              0.54502481     0.087519180       0.4808161
#> 901_NA                0.54325567     0.048653042       0.4157832
#> 902_Rnf113a           1.00000000     1.000000000       1.0000000
#> 903_Nelfa             1.00000000     1.000000000       1.0000000
#> 904_Zkscan1           1.00000000     1.000000000       1.0000000
#> 905_Chaf1b            0.94013089     0.859444878       0.7685402
#> 906_Eif3b             0.64073199     0.848200886       0.6592642
#> 907_Top2b             0.52376249     0.252674563       0.9998830
#> 908_Chaf1b            1.00000000     1.000000000       0.7415917
#> 909_Epb41l2           1.00000000     1.000000000       1.0000000
#> 910_C3H11orf58        1.00000000     0.722042209       0.8402942
#> 911_Top2b             0.64828195     0.017767918       0.1865649
#> 912_Wee1              0.99176908     0.399784993       0.4258505
#> 913_Raly              0.64828195     0.756572162       0.7479895
#> 914_H671_1g2680       0.87530023     0.074667937       0.1973920
#> 915_Eef1b2            1.00000000     1.000000000       1.0000000
#> 916_Gas2l1            0.54325567     1.000000000       1.0000000
#> 917_Epb41l2           0.71607132     0.460414349       0.8639562
#> 918_Rpl23a            0.71607132     0.011705607       0.1316251
#> 919_Chmp2b            0.86772434     0.336096117       0.3779698
#> 920_Lrrfip2           0.95524306     0.553896530       0.7564088
#> 921_Aldoa             0.63454585     0.786737186       0.6260825
#> 922_Cby1              0.97935078     0.462802123       0.6485884
#> 923_LOC100759640      0.88292911     0.187885817       0.2528074
#> 924_Rbm28             0.47623282     0.337140363       0.7073961
#> 925_Skiv2l            0.62288573     0.460926236       0.9914096
#> 926_Ints1             0.97194962     0.243179869       0.3779698
#> 927_Ehd1              0.62388094     0.913964935       0.9905692
#> 928_Nr2f6             0.54325567     0.011493459       0.1772224
#> 929_Top2b             0.86062850     0.240009924       0.4711825
#> 930_Lrrfip2           0.70908163     0.099723851       0.3975070
#> 931_Pip4p2            0.72880676     0.250783169       0.9332729
#> 932_Srp72             0.66687948     0.495954955       0.8049503
#> 933_Mtmr9             0.53917994     0.067817361       0.2515247
#> 934_Gigyf1            0.97878789     0.796500636       0.7079214
#> 935_Rbm7              0.86772434     0.205373417       0.7305634
#> 936_LOC100773565      0.52376249     0.080268242       0.4469415
#> 937_Trim35            0.72613118     0.835850865       0.7166115
#> 938_Cbx8              0.79319446     0.583764133       0.3320750
#> 939_Rplp0             0.52376249     0.219937022       0.5602351
#> 940_Aldoa             0.05820201     0.144987365       0.9905692
#> 941_NA                0.71724978     0.011493459       0.1865649
#> 942_Zyx               0.54502481     0.024200859       0.1316251
#> 943_Psip1             0.70026237     0.025997938       0.2528074
#> 944_Slc7a11           0.52376249     0.011705607       0.1973920
#> 945_Miga2             0.86772434     0.023674400       0.1316251
#> 946_Arhgef6           1.00000000     1.000000000       1.0000000
#> 947_Dlgap4            0.65384212     0.152487189       0.6995470
#> 948_Ampd2             1.00000000     1.000000000       1.0000000
#> 949_Luzp1             0.60118421     0.756572162       0.8095756
#> 950_Camlg             0.92291089     0.368330199       0.6923074
#> 951_Pfkfb3            1.00000000     1.000000000       1.0000000
#> 952_NA                0.75854480     0.296230394       0.8129386
#> 953_Raly              0.93975237     0.603016843       0.8939102
#> 954_Kiaa1143          0.78315628     0.429382386       0.7079214
#> 955_Bcar1             0.46779026     0.011705607       0.3779698
#> 956_Gatad2b           1.00000000     1.000000000       1.0000000
#> 957_Eif4ebp2          0.90069073     0.676683041       0.6563396
#> 958_Fam76b            1.00000000     1.000000000       1.0000000
#> 959_Camlg             0.24684158     0.018982808       0.5956977
#> 960_LOC100754077      1.00000000     1.000000000       0.9998830
#> 961_NA                0.86427681     0.769933446       0.5853013
#> 962_Epb41l2           0.52376249     0.605252603       0.8711279
#> 963_Ankrd34a          0.62479340     0.072705811       0.2792817
#> 964_Zc3h15            0.93975237     0.553896530       0.7660779
#> 965_Def6              0.93975237     0.247766927       0.7816251
#> 966_Srsf6             0.59313749     0.553896530       1.0000000
#> 967_H671_4g11480      0.89639358     0.249506635       0.5341971
#> 968_Top2b             0.97194962     0.839266628       0.8997262
#> 969_LOC100769471      0.46779026     0.268195707       0.9661887
#> 970_Raver1            0.62288573     0.589597408       0.8049503
#> 971_Etv3              0.93408439     0.020193592       0.2528074
#> 972_Psd               0.89679055     0.521352863       0.4950697
#> 973_Usp15             1.00000000     1.000000000       0.8480630
#> 974_Nol7              1.00000000     1.000000000       0.9740321
#> 975_Stk38             0.14868675     0.025997938       0.5757619
#> 976_Smim13            1.00000000     1.000000000       1.0000000
#> 977_Etv3              0.29488064     0.009592916       0.2792817
#> 978_Synm              0.76257173     0.158835844       0.5214596
#> 979_Pwp1              0.98684650     0.650300774       0.8613130
#> 980_Fus               0.89639358     0.834331354       0.8129386
#> 981_Junb              0.92291089     0.839266628       0.7605719
#> 982_Phf8              0.89410612     0.913964935       1.0000000
#> 983_Nelfa             1.00000000     1.000000000       1.0000000
#> 984_Prpf4b            1.00000000     1.000000000       1.0000000
#> 985_Abraxas1          0.47623282     0.009592916       0.2515247
#> 986_Prpf4b            0.84802386     0.352472407       0.8258043
#> 987_Raver1            0.81948788     0.589597408       0.8402942
#> 988_Caap1             0.72815581     0.937520812       0.8184153
#> 989_Rpap3             1.00000000     1.000000000       1.0000000
#> 990_Hsph1             0.89639358     0.813074192       0.3175328
#> 991_LOC100750437      0.97194962     0.399784993       1.0000000
#> 992_Mepce             0.58980148     0.233136373       0.9790592
#> 993_Efs               0.89639358     0.452259090       0.6995470
#> 994_Epb41l2           0.13952113     0.009592916       0.2528074
#> 995_Abcf1             0.99508555     0.438379999       0.9790592
#> 996_NA                0.73942586     0.269847805       0.9790592
#> 997_Eif4ebp2          0.79604118     0.590957303       0.9790592
#> 998_Pfkfb3            0.71607132     0.319858659       0.6236418
#> 999_Hnrnpc            0.78328468     0.913964935       0.9816002
#> 1000_Psmd2            0.85570913     0.752570747       0.8997262
#>                  120_vs_neighbors
#> 1_Top2b                 0.9974120
#> 2_NA                    1.0000000
#> 3_Snip1                 0.9974120
#> 4_Tomm34                0.9974120
#> 5_Pus3                  1.0000000
#> 6_Ints1                 0.9974120
#> 7_Mlh1                  0.9974120
#> 8_LOC100750437          0.9974120
#> 9_Pabpc1                0.9974120
#> 10_Top2b                0.9974120
#> 11_Gorasp1              0.9974120
#> 12_Ints1                1.0000000
#> 13_Syvn1                0.9974120
#> 14_Znf280b              0.9974120
#> 15_Mrnip                0.9974120
#> 16_Rragc                0.9974120
#> 17_Gorasp1              0.9974120
#> 18_Tomm34               0.9974120
#> 19_LOC100757430         0.9974120
#> 20_Ubxn1                0.9974120
#> 21_H671_1g1131          0.9974120
#> 22_Luzp1                0.9974120
#> 23_Efs                  0.9974120
#> 24_Mta2                 0.9974120
#> 25_Nedd1                0.9974120
#> 26_Gigyf1               0.9974120
#> 27_Myh9                 0.9974120
#> 28_Caskin2              0.9974120
#> 29_Papolg               0.9974120
#> 30_Tfg                  0.9974120
#> 31_Rpl34                0.9974120
#> 32_Mideas               0.9974120
#> 33_Gys1                 0.9974120
#> 34_Arhgef6              0.9974120
#> 35_Ctdspl2              1.0000000
#> 36_Ptpn14               0.9974120
#> 37_Raly                 0.9974120
#> 38_Znhit3               0.9974120
#> 39_LOC113833392         1.0000000
#> 40_Luc7l3               0.9974120
#> 41_Rplp0                1.0000000
#> 42_Gys1                 0.9974120
#> 43_Rpl22l1              0.9974120
#> 44_Eif3b                0.9974120
#> 45_Med26                0.9974120
#> 46_Mepce                0.9974120
#> 47_Pdcd11               0.9974120
#> 48_Twf1                 0.9974120
#> 49_LOC100759640         1.0000000
#> 50_Wrnip1               0.9974120
#> 51_Poldip3              0.9974120
#> 52_Ampd2                0.9974120
#> 53_Mea1                 0.9974120
#> 54_Dbn1                 0.9974120
#> 55_Snip1                0.9974120
#> 56_Srsf6                0.9974120
#> 57_LOC113834282         0.9974120
#> 58_Map9                 0.9974120
#> 59_Cdc42ep1             0.9974120
#> 60_Poldip3              1.0000000
#> 61_LOC100764225         1.0000000
#> 62_Epb41l2              0.9974120
#> 63_H671_4g11480         0.9974120
#> 64_Nbn                  0.9974120
#> 65_U2surp               0.9974120
#> 66_Gigyf1               0.9974120
#> 67_NA                   0.9974120
#> 68_Luc7l3               1.0000000
#> 69_LOC100752363         0.9974120
#> 70_Ampd2                1.0000000
#> 71_LOC100759640         1.0000000
#> 72_Stam                 0.9974120
#> 73_Nsfl1c               0.9974120
#> 74_Pfkfb3               0.9974120
#> 75_Rad23a               0.9974120
#> 76_Elf2                 0.9974120
#> 77_Crem                 0.9974120
#> 78_Rragc                0.9974120
#> 79_Lrrfip2              0.9974120
#> 80_Zyx                  0.9974120
#> 81_Lrrfip2              0.9974120
#> 82_Gatad2b              0.9995772
#> 83_Bcar1                0.9974120
#> 84_Ehd1                 0.9974120
#> 85_LOC113834282         0.9974120
#> 86_Tmem230              0.9974120
#> 87_Ncbp1                0.9974120
#> 88_Mllt1                1.0000000
#> 89_Stk17b               0.9974120
#> 90_Dlgap4               0.9974120
#> 91_Papolg               0.9974120
#> 92_Cyld                 0.9974120
#> 93_Gigyf1               0.9974120
#> 94_Lrrfip2              1.0000000
#> 95_Lrrfip2              0.9974120
#> 96_Rlim                 0.9974120
#> 97_Eif3b                0.9974120
#> 98_Mphosph10            1.0000000
#> 99_Gatad2b              0.9974120
#> 100_Srsf6               0.9974120
#> 101_Zyx                 0.9974120
#> 102_Mphosph10           1.0000000
#> 103_Psip1               0.9974120
#> 104_Fbl                 0.9974120
#> 105_H671_1g2680         0.9974120
#> 106_Sgtb                0.9974120
#> 107_Gnl3                0.9974120
#> 108_Eif3b               0.9974120
#> 109_Serpinb1            1.0000000
#> 110_N4bp1               0.9974120
#> 111_Snip1               1.0000000
#> 112_Psip1               1.0000000
#> 113_Mlh1                0.9974120
#> 114_Bsg                 0.9974120
#> 115_Tnpo1               0.9974120
#> 116_H671_1g2680         1.0000000
#> 117_Cbx8                1.0000000
#> 118_Mideas              1.0000000
#> 119_Mideas              0.9974120
#> 120_Dcun1d3             0.9974120
#> 121_Dlg1                1.0000000
#> 122_Rad23a              0.9974120
#> 123_Srsf6               0.9974120
#> 124_Stx7                0.9995772
#> 125_Pdcd11              0.9974120
#> 126_Kiaa1958            0.9995772
#> 127_Pwp1                0.9974120
#> 128_Txlng               0.9974120
#> 129_Junb                1.0000000
#> 130_LOC100759640        1.0000000
#> 131_Dbn1                0.9974120
#> 132_Top2b               0.9974120
#> 133_Rusc2               0.9974120
#> 134_NA                  1.0000000
#> 135_LOC113837251        0.9974120
#> 136_Fam76b              1.0000000
#> 137_Ptpn14              0.9974120
#> 138_Chmp4b              0.9974120
#> 139_Prpf4b              1.0000000
#> 140_Eif3b               1.0000000
#> 141_Nsfl1c              1.0000000
#> 142_Pdlim7              1.0000000
#> 143_Rnf113a             1.0000000
#> 144_Epb41l2             0.9974120
#> 145_Hnrnpc              1.0000000
#> 146_LOC113834282        0.9974120
#> 147_Plekho2             0.9974120
#> 148_Med26               0.9974120
#> 149_Arhgef40            0.9974120
#> 150_NA                  0.9974120
#> 151_Phf8                0.9974120
#> 152_Minar1              1.0000000
#> 153_H671_21690          0.9974120
#> 154_Arhgef40            0.9974120
#> 155_Chaf1b              0.9974120
#> 156_Prpf4b              0.9974120
#> 157_Znf367              0.9974120
#> 158_Luzp1               0.9974120
#> 159_LOC113833882        0.9974120
#> 160_Hnrnpc              0.9974120
#> 161_Mepce               0.9974120
#> 162_Ubxn1               0.9974120
#> 163_Mllt1               0.9974120
#> 164_Chaf1b              0.9974120
#> 165_Raly                0.9974120
#> 166_Gas2l1              0.9974120
#> 167_Dlg1                0.9974120
#> 168_Hoxc10              0.9974120
#> 169_Gigyf1              1.0000000
#> 170_Luzp1               0.9974120
#> 171_Srp72               0.9974120
#> 172_LOC100771461        1.0000000
#> 173_Chaf1b              0.9974120
#> 174_C3H11orf58          1.0000000
#> 175_Pdcd11              0.9974120
#> 176_Psip1               0.9974120
#> 177_Prpf4b              0.9974120
#> 178_Rnf113a             0.9974120
#> 179_Irf3                0.9974120
#> 180_Smim13              0.9974120
#> 181_Gnl3                0.9974120
#> 182_Psma5               0.9974120
#> 183_Ptpn14              0.9974120
#> 184_Prpf4b              1.0000000
#> 185_Top2b               1.0000000
#> 186_Prpf38b             0.9974120
#> 187_Epb41l2             0.9974120
#> 188_Eif3b               1.0000000
#> 189_Hnrnpc              0.9974120
#> 190_LOC100758278        0.9974120
#> 191_Prpf4b              0.9974120
#> 192_Caskin2             0.9974120
#> 193_LOC100752363        0.9974120
#> 194_Septin6             0.9974120
#> 195_Max                 0.9974120
#> 196_Mid1ip1             0.9974120
#> 197_NA                  0.9974120
#> 198_Hsph1               1.0000000
#> 199_Nol7                0.9974120
#> 200_Raly                0.9974120
#> 201_Smim13              0.9974120
#> 202_LOC100757535        0.9974120
#> 203_Net1                0.9974120
#> 204_LOC100754077        1.0000000
#> 205_Snip1               0.9974120
#> 206_Hnrnpc              0.9974120
#> 207_Ldlrap1             0.9974120
#> 208_Luzp1               0.9974120
#> 209_Rpl26               0.9974120
#> 210_Epb41l2             0.9974120
#> 211_Znf367              0.9974120
#> 212_Dlgap4              1.0000000
#> 213_Plekho2             1.0000000
#> 214_Zpr1                0.9974120
#> 215_Dlgap4              1.0000000
#> 216_Def6                0.9974120
#> 217_Eif4ebp2            0.9974120
#> 218_Eef1b2              0.9974120
#> 219_Rad23a              0.9974120
#> 220_Morf4l2             0.9974120
#> 221_Arhgef40            0.9974120
#> 222_NA                  0.9974120
#> 223_LOC100773565        0.9974120
#> 224_Dus2                0.9974120
#> 225_Pip4p2              0.9974120
#> 226_Top2b               0.9974120
#> 227_Znf280b             0.9974120
#> 228_Pdcd11              0.9974120
#> 229_Bckdk               1.0000000
#> 230_Arhgef40            1.0000000
#> 231_Mepce               1.0000000
#> 232_Ccnd3               0.9974120
#> 233_Phf8                0.9974120
#> 234_H671_1g2680         0.9974120
#> 235_Ell                 0.9974120
#> 236_U2surp              0.9974120
#> 237_Rps10               0.9974120
#> 238_Ctdspl2             0.9974120
#> 239_Top2b               0.9974120
#> 240_Msantd3             0.9974120
#> 241_Fam76b              0.9974120
#> 242_Ppp4r3a             0.9974120
#> 243_Gpatch4             0.9974120
#> 244_Nudc                0.9974120
#> 245_Nol7                0.9974120
#> 246_Plekho2             0.9974120
#> 247_Prpf4b              0.9974120
#> 248_Mta2                0.9974120
#> 249_U2surp              0.9974120
#> 250_Ubxn1               1.0000000
#> 251_Rlim                0.9974120
#> 252_Atat1               0.9974120
#> 253_Ubxn1               0.9995772
#> 254_H671_1g2680         0.9974120
#> 255_eIF2aK2             0.9974120
#> 256_Skiv2l              1.0000000
#> 257_Rpl28               1.0000000
#> 258_LOC100759640        1.0000000
#> 259_Gatad2b             0.9974120
#> 260_NA                  0.9974120
#> 261_Gprasp1             0.9974120
#> 262_Luzp1               1.0000000
#> 263_Slc1a5              1.0000000
#> 264_LOC113834282        0.9974120
#> 265_Srsf6               0.9974120
#> 266_Cdc42ep1            0.9974120
#> 267_Net1                0.9974120
#> 268_Caskin2             0.9974120
#> 269_LOC100759640        0.9974120
#> 270_Mideas              0.9974120
#> 271_Luzp1               0.9974120
#> 272_Emd                 1.0000000
#> 273_Plpp6               0.9974120
#> 274_LOC100759640        1.0000000
#> 275_Rps7                0.9974120
#> 276_Fkbp1a              0.9974120
#> 277_Gatad2b             0.9974120
#> 278_Znf385a             0.9974120
#> 279_Arhgef6             0.9974120
#> 280_Slirp               0.9974120
#> 281_Skiv2l              0.9974120
#> 282_H671_21690          0.9974120
#> 283_Kat8                0.9974120
#> 284_Nkap                0.9974120
#> 285_Gsk3b               0.9974120
#> 286_Ints1               0.9974120
#> 287_Gas2l1              1.0000000
#> 288_LOC100759640        1.0000000
#> 289_Top2b               0.9974120
#> 290_Kif20b              0.9974120
#> 291_Phf8                1.0000000
#> 292_Snip1               0.9974120
#> 293_Gsk3b               1.0000000
#> 294_Caskin2             1.0000000
#> 295_C3H11orf58          0.9974120
#> 296_Lrch4               0.9974120
#> 297_LOC113834282        0.9974120
#> 298_LOC100750407        1.0000000
#> 299_LOC113833392        1.0000000
#> 300_LOC113833882        0.9974120
#> 301_Ldlrap1             0.9974120
#> 302_Wee1                0.9974120
#> 303_Caap1               1.0000000
#> 304_Eif4ebp2            0.9974120
#> 305_Ripk2               0.9974120
#> 306_Srp72               0.9974120
#> 307_Taok2               0.9974120
#> 308_Nr2f6               0.9974120
#> 309_Arhgef40            0.9974120
#> 310_Gys1                1.0000000
#> 311_Dlg1                1.0000000
#> 312_Vapb                0.9974120
#> 313_LOC100757535        0.9974120
#> 314_Mkrn2               0.9974120
#> 315_Eif3b               0.9974120
#> 316_Isyna1              1.0000000
#> 317_Prpf4b              0.9974120
#> 318_LOC113833882        0.9974120
#> 319_Lrch4               0.9974120
#> 320_Dbn1                0.9974120
#> 321_Abcf1               0.9974120
#> 322_Ints1               0.9974120
#> 323_C3H11orf58          1.0000000
#> 324_Psma5               0.9974120
#> 325_Fundc1              0.9974120
#> 326_Papolg              1.0000000
#> 327_Mideas              0.9974120
#> 328_Ubxn1               0.9974120
#> 329_Synm                0.9974120
#> 330_Arhgef6             1.0000000
#> 331_Ptpn14              0.9974120
#> 332_Pgrmc1              0.9974120
#> 333_Myh9                0.9974120
#> 334_Etv3                1.0000000
#> 335_Ip6k1               1.0000000
#> 336_Luzp1               1.0000000
#> 337_Ptpn14              1.0000000
#> 338_Caskin2             1.0000000
#> 339_Chaf1b              0.9974120
#> 340_Ubxn1               1.0000000
#> 341_Ube2c               1.0000000
#> 342_Gins2               0.9974120
#> 343_Nlgn2               0.9974120
#> 344_Nf2                 0.9974120
#> 345_Pip4p2              0.9974120
#> 346_Emd                 0.9974120
#> 347_Top2b               0.9974120
#> 348_Trim35              0.9974120
#> 349_NA                  0.9974120
#> 350_NA                  1.0000000
#> 351_Mideas              0.9974120
#> 352_Gas2l1              0.9974120
#> 353_Ampd2               0.9974120
#> 354_Calu                1.0000000
#> 355_Fam76b              1.0000000
#> 356_Dlg1                0.9995772
#> 357_Srsf6               0.9974120
#> 358_Chaf1b              0.9974120
#> 359_Dbn1                0.9974120
#> 360_Tcf25               0.9974120
#> 361_Psip1               0.9974120
#> 362_Cnpy3               0.9974120
#> 363_LOC100759640        0.9974120
#> 364_Zyx                 0.9974120
#> 365_Lrch4               0.9974120
#> 366_Bola1               1.0000000
#> 367_Znf385a             0.9974120
#> 368_Kif20b              0.9974120
#> 369_Ell                 0.9974120
#> 370_Ell                 0.9974120
#> 371_Srsf6               0.9974120
#> 372_Pwp1                0.9974120
#> 373_Def6                0.9974120
#> 374_Cbx8                1.0000000
#> 375_Ddx51               0.9974120
#> 376_Psip1               0.9974120
#> 377_Arhgef40            0.9974120
#> 378_Raly                0.9974120
#> 379_NA                  0.9974120
#> 380_Lrrfip2             0.9974120
#> 381_Gnl3                0.9974120
#> 382_Caskin2             0.9974120
#> 383_Rragc               0.9974120
#> 384_Caskin2             0.9974120
#> 385_Bcar1               1.0000000
#> 386_Homer3              0.9974120
#> 387_Luzp1               0.9974120
#> 388_N4bp1               1.0000000
#> 389_Ppp4r3a             0.9974120
#> 390_H671_1g2680         0.9974120
#> 391_Gnl3                0.9974120
#> 392_Top2b               0.9974120
#> 393_Oser1               1.0000000
#> 394_Snrk                0.9974120
#> 395_Kat8                0.9974120
#> 396_Raver1              1.0000000
#> 397_Pdcd11              0.9974120
#> 398_Rps20               0.9974120
#> 399_Bsg                 0.9974120
#> 400_Raly                0.9974120
#> 401_Pdcd2               0.9974120
#> 402_Caskin2             0.9974120
#> 403_LOC100773571        0.9974120
#> 404_Papolg              1.0000000
#> 405_LOC100757535        0.9974120
#> 406_Caap1               0.9974120
#> 407_Psip1               0.9974120
#> 408_Dbn1                0.9974120
#> 409_Mta2                0.9995772
#> 410_Abcf1               1.0000000
#> 411_LOC100754108        0.9974120
#> 412_Slirp               0.9974120
#> 413_Nelfa               0.9974120
#> 414_Aggf1               0.9974120
#> 415_Bap1                0.9974120
#> 416_Luc7l3              0.9974120
#> 417_Rrp1                0.9974120
#> 418_Wrnip1              0.9974120
#> 419_NA                  0.9974120
#> 420_Abcf1               0.9974120
#> 421_Cluap1              0.9974120
#> 422_Hnrnpc              0.9974120
#> 423_Ptpn1               0.9974120
#> 424_Myh9                0.9974120
#> 425_U2surp              0.9974120
#> 426_NA                  0.9974120
#> 427_Arhgef40            1.0000000
#> 428_Chaf1b              0.9974120
#> 429_Prpf4b              0.9974120
#> 430_Epb41l2             0.9974120
#> 431_Eif3b               1.0000000
#> 432_Isyna1              0.9974120
#> 433_U2surp              0.9974120
#> 434_LOC100765020        0.9974120
#> 435_Arhgef6             0.9974120
#> 436_Ptpn1               0.9974120
#> 437_Prpf4b              0.9974120
#> 438_Rpl35a              0.9974120
#> 439_Prpf4b              1.0000000
#> 440_Zyx                 0.9974120
#> 441_Dbn1                0.9974120
#> 442_Chaf1b              1.0000000
#> 443_LOC113834282        0.9974120
#> 444_Gpsm2               0.9974120
#> 445_LOC100757535        1.0000000
#> 446_Cfap410             0.9974120
#> 447_Epb41l2             0.9974120
#> 448_Ncbp1               1.0000000
#> 449_Pacsin1             1.0000000
#> 450_Cstf2               0.9974120
#> 451_LOC100769437        0.9974120
#> 452_eIF2aK2             0.9974120
#> 453_Kiaa1191            0.9974120
#> 454_Mepce               0.9974120
#> 455_Cbx8                0.9974120
#> 456_Eed                 0.9974120
#> 457_Cdc42ep1            0.9974120
#> 458_Lrrfip2             1.0000000
#> 459_Pacsin1             0.9974120
#> 460_Gpatch4             0.9974120
#> 461_Plin4               1.0000000
#> 462_NA                  0.9974120
#> 463_Snip1               0.9974120
#> 464_Cyld                0.9974120
#> 465_Plin4               0.9974120
#> 466_Twf1                1.0000000
#> 467_LOC113834282        1.0000000
#> 468_Snip1               0.9974120
#> 469_Ppp4r3a             0.9974120
#> 470_Psip1               0.9974120
#> 471_Dnajc5              1.0000000
#> 472_Phf8                0.9974120
#> 473_Bola1               0.9974120
#> 474_Cdc42ep1            0.9974120
#> 475_Eif4ebp2            0.9974120
#> 476_Prpf38b             0.9974120
#> 477_Klhl26              0.9974120
#> 478_Hsph1               1.0000000
#> 479_Snip1               0.9974120
#> 480_Caskin2             0.9974120
#> 481_Plpp6               0.9974120
#> 482_NA                  0.9974120
#> 483_Mlh1                0.9974120
#> 484_Gys1                0.9995772
#> 485_Tfg                 0.9974120
#> 486_Arhgef6             1.0000000
#> 487_Mphosph10           1.0000000
#> 488_Hoxc10              0.9974120
#> 489_LOC100759640        1.0000000
#> 490_Arhgef40            0.9974120
#> 491_Dnajc5              1.0000000
#> 492_Tbc1d23             1.0000000
#> 493_Ubxn1               0.9974120
#> 494_Rab1a               1.0000000
#> 495_Eif3b               1.0000000
#> 496_Tceal8              0.9974120
#> 497_Dlgap4              0.9974120
#> 498_Smim13              0.9974120
#> 499_NA                  1.0000000
#> 500_Lrch4               1.0000000
#> 501_Bola1               1.0000000
#> 502_NA                  0.9974120
#> 503_Ptpn14              0.9974120
#> 504_LOC100759640        0.9974120
#> 505_Rps10               0.9974120
#> 506_Top2b               1.0000000
#> 507_Ssr3                0.9974120
#> 508_Homer3              1.0000000
#> 509_Phf8                0.9974120
#> 510_LOC100767716        0.9974120
#> 511_Xpa                 0.9974120
#> 512_H671_21690          0.9974120
#> 513_LOC100769471        0.9974120
#> 514_Gas2l1              0.9974120
#> 515_Luzp1               0.9974120
#> 516_Gpbp1               0.9974120
#> 517_Gatad2b             1.0000000
#> 518_Gys1                0.9974120
#> 519_Top2b               0.9974120
#> 520_LOC100757535        0.9974120
#> 521_Lpcat4              0.9974120
#> 522_Arhgef6             1.0000000
#> 523_Cavin3              0.9974120
#> 524_Gpatch4             0.9974120
#> 525_Prpf38b             0.9974120
#> 526_Timm8a              0.9974120
#> 527_Cavin3              0.9974120
#> 528_Mkrn2               0.9974120
#> 529_Oser1               0.9974120
#> 530_Gsk3b               0.9974120
#> 531_Eef1b2              1.0000000
#> 532_Ampd2               0.9974120
#> 533_Lrrfip2             0.9974120
#> 534_Ring1               0.9974120
#> 535_Rlim                0.9974120
#> 536_LOC100759640        0.9974120
#> 537_LOC100759640        0.9974120
#> 538_Atp5pf              0.9974120
#> 539_Max                 0.9974120
#> 540_Bap1                1.0000000
#> 541_Nsfl1c              0.9974120
#> 542_Prpf4b              1.0000000
#> 543_LOC100757535        0.9974120
#> 544_Mtmr10              0.9974120
#> 545_Hoxc10              0.9974120
#> 546_Trim35              1.0000000
#> 547_Eif4ebp2            0.9974120
#> 548_Dlgap4              1.0000000
#> 549_Gys1                0.9974120
#> 550_Sgtb                0.9974120
#> 551_Eri2                0.9974120
#> 552_Ccnd3               1.0000000
#> 553_Smim13              0.9974120
#> 554_Snrk                0.9974120
#> 555_Caskin2             1.0000000
#> 556_Pdcd11              0.9974120
#> 557_Pgam5               0.9974120
#> 558_Mphosph10           1.0000000
#> 559_Mideas              1.0000000
#> 560_Top2b               0.9974120
#> 561_LOC100763014        0.9974120
#> 562_Snip1               1.0000000
#> 563_Ubxn1               0.9974120
#> 564_LOC100750407        0.9974120
#> 565_Morf4l2             1.0000000
#> 566_Ctdspl2             1.0000000
#> 567_Cwf19l1             1.0000000
#> 568_Eef1b2              1.0000000
#> 569_C1H12orf45          1.0000000
#> 570_Znf367              0.9974120
#> 571_Ankrd34a            0.9974120
#> 572_Mllt11              0.9974120
#> 573_LOC100774792        0.9974120
#> 574_NA                  0.9974120
#> 575_Cbx8                0.9974120
#> 576_Bckdk               0.9974120
#> 577_Snip1               1.0000000
#> 578_Nsfl1c              0.9974120
#> 579_Gas2l1              1.0000000
#> 580_Nudc                0.9974120
#> 581_Epb41l2             0.9974120
#> 582_Mtmr6               1.0000000
#> 583_Znf668              0.9974120
#> 584_Hsph1               0.9974120
#> 585_LOC113834282        0.9974120
#> 586_Ctdspl2             0.9974120
#> 587_Foxf1               0.9974120
#> 588_Luzp1               0.9974120
#> 589_Xpa                 1.0000000
#> 590_Psip1               0.9974120
#> 591_Rbm7                0.9974120
#> 592_Mtrex               0.9974120
#> 593_Arhgef40            0.9974120
#> 594_Plekho2             0.9974120
#> 595_Bckdk               0.9974120
#> 596_Dut                 1.0000000
#> 597_Abcf1               1.0000000
#> 598_Txnl1               1.0000000
#> 599_Nudc                1.0000000
#> 600_Sh3gl1              0.9974120
#> 601_Gatad2b             0.9974120
#> 602_Homer3              0.9974120
#> 603_Septin6             0.9974120
#> 604_Smim13              0.9974120
#> 605_Arhgef40            0.9974120
#> 606_Rpl32               1.0000000
#> 607_Tomm34              0.9974120
#> 608_Mlh1                0.9974120
#> 609_Tbcc                1.0000000
#> 610_Eif3d               0.9974120
#> 611_Snrk                0.9974120
#> 612_Bckdk               1.0000000
#> 613_Wdr3                0.9974120
#> 614_LOC100757535        1.0000000
#> 615_Dlg1                0.9974120
#> 616_LOC100767716        0.9974120
#> 617_Hnrnpc              0.9974120
#> 618_Mphosph10           0.9974120
#> 619_Eif3b               0.9974120
#> 620_Emd                 1.0000000
#> 621_Txlng               0.9974120
#> 622_Prpf4b              1.0000000
#> 623_Rlim                1.0000000
#> 624_Eef1b2              0.9974120
#> 625_Def6                0.9974120
#> 626_LOC100765020        0.9974120
#> 627_U2surp              0.9974120
#> 628_Elf2                0.9974120
#> 629_Slc1a5              0.9974120
#> 630_NA                  0.9974120
#> 631_Tfg                 1.0000000
#> 632_Top2b               0.9974120
#> 633_Pip4p2              0.9974120
#> 634_Cdc42ep1            1.0000000
#> 635_Hsph1               1.0000000
#> 636_Twf1                1.0000000
#> 637_Nbn                 0.9974120
#> 638_Psmd4               0.9974120
#> 639_Bap1                0.9974120
#> 640_Mepce               0.9974120
#> 641_Mideas              0.9974120
#> 642_LOC100759640        0.9974120
#> 643_Epb41l2             1.0000000
#> 644_Sav1                0.9974120
#> 645_Prpf4b              0.9974120
#> 646_Gnas                0.9974120
#> 647_Mllt1               1.0000000
#> 648_Poldip3             0.9974120
#> 649_Aldoa               0.9974120
#> 650_Rbbp8               0.9974120
#> 651_LOC113834282        0.9974120
#> 652_Gys1                0.9974120
#> 653_Hnrnpc              0.9974120
#> 654_Vps35               0.9974120
#> 655_Miga2               0.9974120
#> 656_Epb41l2             0.9974120
#> 657_Tob2                1.0000000
#> 658_Lamtor1             1.0000000
#> 659_LOC100759640        1.0000000
#> 660_Epb41l2             0.9974120
#> 661_Rlim                0.9974120
#> 662_Gys1                1.0000000
#> 663_LOC100750437        0.9974120
#> 664_NA                  0.9974120
#> 665_Nbn                 0.9974120
#> 666_Tyw3                1.0000000
#> 667_Gas2l1              0.9974120
#> 668_Fus                 1.0000000
#> 669_Prpf38b             0.9974120
#> 670_Calu                1.0000000
#> 671_Rras2               0.9974120
#> 672_Prpf4b              0.9974120
#> 673_Nelfa               0.9974120
#> 674_LOC100754077        0.9974120
#> 675_Rbm28               0.9974120
#> 676_Nsfl1c              0.9974120
#> 677_Rnf126              0.9974120
#> 678_Eme1                0.9974120
#> 679_Nbn                 0.9974120
#> 680_Eif4ebp2            0.9974120
#> 681_Wee1                1.0000000
#> 682_Prpf38b             1.0000000
#> 683_Luzp1               1.0000000
#> 684_Gas2l1              1.0000000
#> 685_Pdcd11              0.9974120
#> 686_Chaf1b              0.9974120
#> 687_Pycr1               0.9974120
#> 688_Phf8                0.9974120
#> 689_Raver1              1.0000000
#> 690_Dbn1                0.9974120
#> 691_Dut                 0.9974120
#> 692_Prpf4b              1.0000000
#> 693_Prpf4b              0.9974120
#> 694_Efs                 0.9974120
#> 695_NA                  0.9974120
#> 696_Ppp2r5b             1.0000000
#> 697_Caskin2             0.9974120
#> 698_Arhgef40            0.9974120
#> 699_Zyx                 1.0000000
#> 700_Mphosph10           0.9974120
#> 701_LOC113833392        0.9974120
#> 702_Cdc42ep1            0.9974120
#> 703_Snrpa1              1.0000000
#> 704_Ncbp1               0.9974120
#> 705_Gas2l1              1.0000000
#> 706_Gas2l1              0.9974120
#> 707_Bap1                0.9974120
#> 708_LOC100759640        0.9974120
#> 709_Cherp               1.0000000
#> 710_Nbn                 0.9974120
#> 711_LOC100759640        1.0000000
#> 712_NA                  1.0000000
#> 713_Eif3b               0.9974120
#> 714_Miga2               1.0000000
#> 715_Prpf4b              0.9974120
#> 716_Dbn1                0.9974120
#> 717_Ppp2r5b             0.9974120
#> 718_Exosc9              1.0000000
#> 719_Eif3b               0.9974120
#> 720_Ripk2               1.0000000
#> 721_Dlg1                0.9974120
#> 722_N4bp1               0.9974120
#> 723_Nudc                0.9974120
#> 724_Znf367              0.9974120
#> 725_Ring1               1.0000000
#> 726_Snrpa1              0.9974120
#> 727_U2surp              0.9974120
#> 728_LOC100764225        0.9974120
#> 729_Cdc42ep1            1.0000000
#> 730_Znf385a             1.0000000
#> 731_Ints1               1.0000000
#> 732_LOC113833392        1.0000000
#> 733_Lrch4               0.9974120
#> 734_Ctdspl2             0.9974120
#> 735_Prpf4b              0.9974120
#> 736_Luzp1               0.9974120
#> 737_Eif3b               1.0000000
#> 738_Ptpn14              1.0000000
#> 739_Rrp1                0.9974120
#> 740_Lrrfip2             1.0000000
#> 741_Nsfl1c              1.0000000
#> 742_Ddx51               0.9974120
#> 743_Prpf38b             0.9974120
#> 744_Eef1b2              0.9974120
#> 745_Znf385a             1.0000000
#> 746_Map9                1.0000000
#> 747_Rflnb               1.0000000
#> 748_NA                  0.9974120
#> 749_C1H12orf45          0.9974120
#> 750_U2surp              0.9974120
#> 751_Caskin2             0.9974120
#> 752_Eri1                0.9974120
#> 753_Gsk3b               0.9974120
#> 754_LOC100766946        0.9974120
#> 755_Cnpy3               0.9974120
#> 756_Hnrnpc              0.9974120
#> 757_Ptpn14              0.9974120
#> 758_Slc7a11             0.9974120
#> 759_Hnrnpc              0.9974120
#> 760_Cdc37l1             0.9974120
#> 761_LOC100768405        0.9974120
#> 762_Rragc               0.9974120
#> 763_LOC113834282        0.9974120
#> 764_Fus                 0.9974120
#> 765_Ubxn1               0.9974120
#> 766_Mmut                0.9974120
#> 767_Pdcd11              0.9974120
#> 768_LOC100757535        0.9974120
#> 769_Eif3b               0.9974120
#> 770_Rnf113a             1.0000000
#> 771_Sytl4               1.0000000
#> 772_Tlnrd1              0.9974120
#> 773_H671_1g1131         0.9974120
#> 774_Neurl1              0.9974120
#> 775_Zyx                 0.9974120
#> 776_Ctdspl2             0.9974120
#> 777_Chaf1b              0.9995772
#> 778_Rragc               1.0000000
#> 779_Srfbp1              0.9974120
#> 780_Gys1                0.9974120
#> 781_Usp15               0.9974120
#> 782_Arhgef40            0.9974120
#> 783_Gigyf1              0.9974120
#> 784_Minar1              0.9974120
#> 785_Dus2                1.0000000
#> 786_Gatad2b             0.9974120
#> 787_Eif5                0.9974120
#> 788_Epb41l2             0.9974120
#> 789_Arl6ip4             0.9974120
#> 790_Plin4               1.0000000
#> 791_Elf2                0.9974120
#> 792_Plin4               1.0000000
#> 793_Snip1               1.0000000
#> 794_Txlng               1.0000000
#> 795_LOC100769437        0.9974120
#> 796_Caskin2             0.9974120
#> 797_NA                  0.9974120
#> 798_Synm                1.0000000
#> 799_Synm                0.9974120
#> 800_Ube2c               0.9974120
#> 801_Sgtb                0.9974120
#> 802_Prpf4b              0.9974120
#> 803_Epb41l2             0.9974120
#> 804_Mllt1               0.9974120
#> 805_LOC100759640        1.0000000
#> 806_Epb41l2             0.9974120
#> 807_Znf280b             0.9974120
#> 808_Kiaa1143            1.0000000
#> 809_Gas2l1              1.0000000
#> 810_Srp72               0.9974120
#> 811_Tomm22              0.9974120
#> 812_Psip1               0.9974120
#> 813_Arhgef37            0.9974120
#> 814_Bckdk               0.9974120
#> 815_Strip1              0.9974120
#> 816_Usp15               1.0000000
#> 817_Ssr3                0.9974120
#> 818_Strip1              0.9974120
#> 819_Eif3b               0.9974120
#> 820_U2surp              0.9974120
#> 821_Bend3               0.9974120
#> 822_Rps10               0.9974120
#> 823_Rpl23a              0.9974120
#> 824_Nbn                 0.9974120
#> 825_Rpap3               0.9974120
#> 826_LOC100759640        1.0000000
#> 827_Ric8a               0.9974120
#> 828_Hsph1               0.9974120
#> 829_LOC100759640        0.9974120
#> 830_LOC100757535        0.9974120
#> 831_Gigyf1              0.9974120
#> 832_Dbn1                0.9974120
#> 833_Snrk                1.0000000
#> 834_Prpf38b             0.9974120
#> 835_LOC100766868        0.9974120
#> 836_LOC100766868        0.9974120
#> 837_Wbp11               0.9974120
#> 838_Rusc2               0.6510735
#> 839_Eif3b               0.9974120
#> 840_Ptpn14              1.0000000
#> 841_Rlim                0.9974120
#> 842_Ints1               1.0000000
#> 843_Chaf1b              0.9974120
#> 844_Dlg1                1.0000000
#> 845_Lamtor1             0.9974120
#> 846_Tab1                0.9974120
#> 847_Dbn1                0.9974120
#> 848_Psip1               0.9974120
#> 849_Dbn1                0.9974120
#> 850_Pabpc1              0.9974120
#> 851_Hnrnpc              0.9974120
#> 852_Emd                 0.9974120
#> 853_LOC100764225        0.9974120
#> 854_Nup50               1.0000000
#> 855_Ctcf                0.9974120
#> 856_Raly                0.9974120
#> 857_Bard1               0.9974120
#> 858_Ptpn14              0.9974120
#> 859_LOC100757535        0.9974120
#> 860_Psmd2               1.0000000
#> 861_Junb                0.9974120
#> 862_C1qbp               0.9974120
#> 863_Lrch4               0.9974120
#> 864_CUNH14orf93         1.0000000
#> 865_U2surp              0.9974120
#> 866_Raly                1.0000000
#> 867_LOC100774417        0.9974120
#> 868_Srp72               0.9974120
#> 869_LOC100764225        1.0000000
#> 870_Morf4l2             0.9974120
#> 871_CUNH9orf40          0.9974120
#> 872_Gas2l1              0.9974120
#> 873_Atp5pf              0.9974120
#> 874_Lrrfip2             0.9974120
#> 875_Prpf4b              0.9974120
#> 876_Top2b               1.0000000
#> 877_Mepce               0.9974120
#> 878_Ptpn14              0.9974120
#> 879_Dnajc25             0.9974120
#> 880_Cbx8                0.9974120
#> 881_Synm                0.9974120
#> 882_Def6                0.9974120
#> 883_Gys1                0.9974120
#> 884_Luzp1               0.9974120
#> 885_Synm                0.9974120
#> 886_Snip1               1.0000000
#> 887_Top2b               0.9974120
#> 888_NA                  0.9974120
#> 889_Trim35              0.9974120
#> 890_Znf385a             0.9974120
#> 891_Chaf1b              0.9974120
#> 892_Abcf1               0.9974120
#> 893_Pdcd11              1.0000000
#> 894_Dlg1                0.9974120
#> 895_Dbn1                1.0000000
#> 896_LOC100752363        0.9974120
#> 897_Ppp4r3a             1.0000000
#> 898_Gas2l1              0.9974120
#> 899_Mtmr10              1.0000000
#> 900_Cyld                0.9974120
#> 901_NA                  0.9974120
#> 902_Rnf113a             1.0000000
#> 903_Nelfa               1.0000000
#> 904_Zkscan1             1.0000000
#> 905_Chaf1b              0.9974120
#> 906_Eif3b               0.9974120
#> 907_Top2b               0.9995772
#> 908_Chaf1b              0.9974120
#> 909_Epb41l2             0.9974120
#> 910_C3H11orf58          0.9974120
#> 911_Top2b               0.9974120
#> 912_Wee1                0.9974120
#> 913_Raly                0.9974120
#> 914_H671_1g2680         0.9974120
#> 915_Eef1b2              1.0000000
#> 916_Gas2l1              1.0000000
#> 917_Epb41l2             0.9974120
#> 918_Rpl23a              0.9974120
#> 919_Chmp2b              0.9974120
#> 920_Lrrfip2             0.9974120
#> 921_Aldoa               0.9974120
#> 922_Cby1                0.9974120
#> 923_LOC100759640        0.9974120
#> 924_Rbm28               0.9974120
#> 925_Skiv2l              0.9974120
#> 926_Ints1               0.9974120
#> 927_Ehd1                0.9974120
#> 928_Nr2f6               0.9974120
#> 929_Top2b               0.9974120
#> 930_Lrrfip2             0.9974120
#> 931_Pip4p2              0.9974120
#> 932_Srp72               0.9974120
#> 933_Mtmr9               0.9974120
#> 934_Gigyf1              1.0000000
#> 935_Rbm7                0.9974120
#> 936_LOC100773565        0.9974120
#> 937_Trim35              0.9974120
#> 938_Cbx8                0.9974120
#> 939_Rplp0               0.9974120
#> 940_Aldoa               0.9974120
#> 941_NA                  0.9974120
#> 942_Zyx                 0.6510735
#> 943_Psip1               0.9974120
#> 944_Slc7a11             0.9974120
#> 945_Miga2               0.9974120
#> 946_Arhgef6             1.0000000
#> 947_Dlgap4              0.9974120
#> 948_Ampd2               1.0000000
#> 949_Luzp1               0.9974120
#> 950_Camlg               0.9974120
#> 951_Pfkfb3              1.0000000
#> 952_NA                  0.9974120
#> 953_Raly                0.9974120
#> 954_Kiaa1143            0.9974120
#> 955_Bcar1               0.9974120
#> 956_Gatad2b             0.9974120
#> 957_Eif4ebp2            0.9974120
#> 958_Fam76b              1.0000000
#> 959_Camlg               0.9974120
#> 960_LOC100754077        0.9974120
#> 961_NA                  0.9974120
#> 962_Epb41l2             0.9974120
#> 963_Ankrd34a            0.9995772
#> 964_Zc3h15              0.9974120
#> 965_Def6                0.9974120
#> 966_Srsf6               1.0000000
#> 967_H671_4g11480        0.9974120
#> 968_Top2b               0.9974120
#> 969_LOC100769471        0.9974120
#> 970_Raver1              0.9974120
#> 971_Etv3                0.9974120
#> 972_Psd                 0.9974120
#> 973_Usp15               1.0000000
#> 974_Nol7                0.9974120
#> 975_Stk38               0.9974120
#> 976_Smim13              0.9974120
#> 977_Etv3                0.9974120
#> 978_Synm                0.9974120
#> 979_Pwp1                0.9974120
#> 980_Fus                 0.9974120
#> 981_Junb                0.9974120
#> 982_Phf8                1.0000000
#> 983_Nelfa               0.9974120
#> 984_Prpf4b              0.9974120
#> 985_Abraxas1            0.9974120
#> 986_Prpf4b              0.9974120
#> 987_Raver1              0.9974120
#> 988_Caap1               0.9974120
#> 989_Rpap3               1.0000000
#> 990_Hsph1               0.9974120
#> 991_LOC100750437        1.0000000
#> 992_Mepce               0.9974120
#> 993_Efs                 0.9974120
#> 994_Epb41l2             0.9974120
#> 995_Abcf1               0.9974120
#> 996_NA                  0.9974120
#> 997_Eif4ebp2            0.9974120
#> 998_Pfkfb3              0.9974120
#> 999_Hnrnpc              0.9974120
#> 1000_Psmd2              0.9974120
#> 
#> $Stationary$pvc_pattern_summary
#>   -60 15 60 90 120 240
#> p   0  0 40  0   0   0
#> v   0  0 22  0   0   0
#> b   0  0  1  0   0   0
#> t   0  0  1  0   0   0
#> 
#> $Stationary$plots
#> list()
```

The terminal output of the function gives information about the total
amount of excursions found in each level, and also a breakdown by
timepoint and type of excursion (p = peak, v = valley, b = bottom of
cliff, t = top of cliff, always for the respective timepoint that is
mentioned in that row).
[Here](https://csbg.github.io/SplineOmics_html_reports/pvc_report_PPTX.html)
you can see the resulting HTML report.

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
    #> [1] knitr_1.51        dplyr_1.2.0       SplineOmics_0.4.3
    #> 
    #> loaded via a namespace (and not attached):
    #>   [1] Rdpack_2.6.5             bitops_1.0-9             pbapply_1.7-4           
    #>   [4] writexl_1.5.4            rlang_1.1.7              magrittr_2.0.4          
    #>   [7] clue_0.3-66              GetoptLong_1.1.0         otel_0.2.0              
    #>  [10] matrixStats_1.5.0        compiler_4.5.2           reshape2_1.4.5          
    #>  [13] png_0.1-8                systemfonts_1.3.1        vctrs_0.7.1             
    #>  [16] stringr_1.6.0            pkgconfig_2.0.3          shape_1.4.6.1           
    #>  [19] crayon_1.5.3             fastmap_1.2.0            backports_1.5.0         
    #>  [22] caTools_1.18.3           rmarkdown_2.30           nloptr_2.2.1            
    #>  [25] ragg_1.5.0               purrr_1.2.1              xfun_0.56               
    #>  [28] cachem_1.1.0             jsonlite_2.0.0           progress_1.2.3          
    #>  [31] EnvStats_3.1.0           remaCor_0.0.20           gmp_0.7-5               
    #>  [34] BiocParallel_1.42.2      broom_1.0.12             parallel_4.5.2          
    #>  [37] prettyunits_1.2.0        cluster_2.1.8.1          R6_2.6.1                
    #>  [40] stringi_1.8.7            bslib_0.10.0             RColorBrewer_1.1-3      
    #>  [43] limma_3.64.3             boot_1.3-32              car_3.1-5               
    #>  [46] ClusterR_1.3.6           numDeriv_2016.8-1.1      jquerylib_0.1.4         
    #>  [49] Rcpp_1.1.1               iterators_1.0.14         base64enc_0.1-6         
    #>  [52] IRanges_2.42.0           Matrix_1.7-4             splines_4.5.2           
    #>  [55] tidyselect_1.2.1         rstudioapi_0.18.0        abind_1.4-8             
    #>  [58] yaml_2.3.12              doParallel_1.0.17        gplots_3.3.0            
    #>  [61] codetools_0.2-19         plyr_1.8.9               lmerTest_3.2-0          
    #>  [64] lattice_0.22-5           tibble_3.3.1             withr_3.0.2             
    #>  [67] Biobase_2.68.0           S7_0.2.1                 evaluate_1.0.5          
    #>  [70] desc_1.4.3               zip_2.3.3                circlize_0.4.17         
    #>  [73] pillar_1.11.1            BiocManager_1.30.27      carData_3.0-6           
    #>  [76] KernSmooth_2.23-26       checkmate_2.3.4          renv_1.1.7              
    #>  [79] foreach_1.5.2            stats4_4.5.2             reformulas_0.4.4        
    #>  [82] generics_0.1.4           S4Vectors_0.46.0         hms_1.1.4               
    #>  [85] ggplot2_4.0.2            scales_1.4.0             aod_1.3.3               
    #>  [88] minqa_1.2.8              gtools_3.9.5             RhpcBLASctl_0.23-42     
    #>  [91] glue_1.8.0               tools_4.5.2              fANCOVA_0.6-1           
    #>  [94] variancePartition_1.38.1 lme4_1.1-38              mvtnorm_1.3-3           
    #>  [97] fs_1.6.6                 grid_4.5.2               tidyr_1.3.2             
    #> [100] rbibutils_2.4.1          colorspace_2.1-2         nlme_3.1-168            
    #> [103] Formula_1.2-5            cli_3.6.5                textshaping_1.0.4       
    #> [106] svglite_2.2.2            ComplexHeatmap_2.24.1    corpcor_1.6.10          
    #> [109] gtable_0.3.6             sass_0.4.10              digest_0.6.39           
    #> [112] BiocGenerics_0.54.1      pbkrtest_0.5.5           ggrepel_0.9.6           
    #> [115] rjson_0.2.23             htmlwidgets_1.6.4        farver_2.1.2            
    #> [118] htmltools_0.5.9          pkgdown_2.2.0            lifecycle_1.0.5         
    #> [121] GlobalOptions_0.1.3      statmod_1.5.1            MASS_7.3-65
