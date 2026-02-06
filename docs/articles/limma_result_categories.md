# limma_result_categories

## limma Result Categories

limma analysis results can be divided into three categories, which are
defined in this document:

1.  **Time Effect**: This category focuses on the changes of the feature
    (e.g. protein) value that occur over time within a single condition.

2.  **Average Difference Between Conditions**: This category compares
    the average feature values for levels within a condition, regardless
    of time.

3.  **Interaction Between Condition and Time**: This category examines
    the interaction between time and condition. It identifies features
    whose value changes differently over time depending on the
    condition.

------------------------------------------------------------------------

#### Legend:

- **A hit** is a feature (e.g. a protein) that is significantly changed
  over time.
- **Levels** are the different factors of a condition of the experiment.
  For example, bioreactor phase is a condition, and exponential and
  stationary are levels within that condition.

### Category 1 (time effect)

Temporal pattern within level for a given feature → Hit

##### Example of a Hit

A hit is a feature that shows a clear temporal pattern over time.

![Clear temporal pattern over
time](limma_result_categories_files/figure-html/unnamed-chunk-1-1.png)

Clear temporal pattern over time

##### Example of No Hit

![No clear temporal
pattern](limma_result_categories_files/figure-html/unnamed-chunk-2-1.png)

No clear temporal pattern

### Category 2 (average difference conditions)

Overall mean difference between levels for a given feature → Hit

##### Example of a Hit

![No clear temporal pattern for both levels but overall mean difference
between
them](limma_result_categories_files/figure-html/unnamed-chunk-3-1.png)

No clear temporal pattern for both levels but overall mean difference
between them

##### Example of No Hit

![Clear temporal pattern for both levels but no overall mean difference
in feature value between
them.](limma_result_categories_files/figure-html/unnamed-chunk-4-1.png)

Clear temporal pattern for both levels but no overall mean difference in
feature value between them.

### Category 3 (interaction condition & time)

Treatment interacting with time for a feature (time effect changing with
treatment, the feature must have different temporal patterns in both
conditions/levels) → Hit

##### Examples of Hits

Different temporal patterns are observed for each level –\> hit in
category 3.

![Different temporal patterns of the feature for both
levels.](limma_result_categories_files/figure-html/unnamed-chunk-5-1.png)

Different temporal patterns of the feature for both levels.

![Different temporal patterns of the feature for both
levels.](limma_result_categories_files/figure-html/unnamed-chunk-6-1.png)

Different temporal patterns of the feature for both levels.

##### Example of No Hit

![Overall the same temporal pattern of a feature for both
levels.](limma_result_categories_files/figure-html/unnamed-chunk-7-1.png)

Overall the same temporal pattern of a feature for both levels.

## Session Info

    ## R version 4.5.2 (2025-10-31)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 22.04.5 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0  LAPACK version 3.10.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=de_AT.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=de_AT.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=de_AT.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=de_AT.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: Europe/Vienna
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices datasets  utils     methods   base     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] desc_1.4.3          digest_0.6.39       R6_2.6.1           
    ##  [4] fastmap_1.2.0       xfun_0.56           cachem_1.1.0       
    ##  [7] knitr_1.51          htmltools_0.5.9     rmarkdown_2.30     
    ## [10] lifecycle_1.0.5     cli_3.6.5           sass_0.4.10        
    ## [13] pkgdown_2.2.0       textshaping_1.0.4   jquerylib_0.1.4    
    ## [16] renv_1.1.7          systemfonts_1.3.1   compiler_4.5.2     
    ## [19] rstudioapi_0.18.0   tools_4.5.2         ragg_1.5.0         
    ## [22] bslib_0.10.0        evaluate_1.0.5      yaml_2.3.12        
    ## [25] otel_0.2.0          BiocManager_1.30.27 jsonlite_2.0.0     
    ## [28] htmlwidgets_1.6.4   rlang_1.1.7         fs_1.6.6
