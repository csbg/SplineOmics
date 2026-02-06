# Designing a Limma Design Formula

## Introduction

The `limma` package is a powerful tool for analyzing gene expression
data, particularly in the context of microarrays and RNA-seq. A critical
part of any `limma` analysis is the design formula, which specifies the
experimental conditions and contrasts that you are interested in. This
vignette provides a guide on how to construct a `limma` design formula
correctly, with examples and best practices.

## Understanding the Design Matrix

The design matrix is a crucial component in any differential expression
analysis using `limma`. It defines the relationships between the samples
and the experimental conditions (factors) under investigation. A
well-constructed design matrix allows `limma` to correctly model the
effects of these factors and estimate the differential expression.

### Basic Design Formula

In its simplest form, the design formula includes only one factor, such
as treatment vs. control. If your experiment involves comparing two
conditions (e.g., `treated` vs. `untreated`), you can create a design
formula like this:

``` r
design <- model.matrix(~ 0 + condition, data = meta)
```

Here, condition is a factor variable in your metadata (meta) that
represents the experimental groups.

Important Points: The ~ 0 + condition syntax tells R to create a design
matrix without an intercept (i.e., a matrix where each level of the
factor condition is represented by its own column). This approach is
helpful when you want to make direct comparisons between conditions.
Including Multiple Factors If your experiment includes more than one
factor, such as time points and treatments, you can include both in the
design formula:

``` r
design <- model.matrix(~ 0 + treatment + time, data = meta)
```

This formula assumes that the effects of treatment and time are additive
(no interaction). If you suspect that the interaction between treatment
and time might be important, you can include an interaction term:

``` r
design <- model.matrix(~ 0 + treatment * time, data = meta)
```

Interaction Term: The treatment \* time term includes both the main
effects of treatment and time and their interaction. Blocking Factors In
some experiments, you might have technical or biological replicates, or
other blocking factors (e.g., batch effects). You can include these
blocking factors in your design formula:

``` r
design <- model.matrix(~ 0 + treatment + batch, data = meta)
```

This formula accounts for both treatment and batch effects, ensuring
that your analysis is not confounded by batch effects.

## Creating Contrasts

After defining your design matrix, you will likely want to make specific
comparisons between conditions. This is where contrasts come in. For
example, to compare treated vs. untreated, you can define a contrast
matrix:

``` r
contrast <- makeContrasts(
  treated_vs_untreated = treatmenttreated - treatmentuntreated,
  levels = design
)
```

## Practical Example

Let’s say you have an experiment with two treatments (A and B) and two
time points (early and late). Your metadata might look like this:

``` r
meta <- data.frame(
  sample = c("S1", "S2", "S3", "S4"),
  treatment = factor(c("A", "A", "B", "B")),
  time = factor(c("early", "late", "early", "late"))
)
```

Your design formula would be:

``` r
design <- model.matrix(~ 0 + treatment * time, data = meta)
```

And a contrast to compare treatment A at the early time point against
treatment B at the late time point would be:

``` r
contrast <- makeContrasts(
  A_early_vs_B_late = (treatmentA:timeearly) - (treatmentB:timelate),
  levels = design
)
```

## Summary

1.  **Intercept or No Intercept:**
    - Starting with `~ 0` means no intercept (i.e., not including a
      baseline group in the model).
    - Starting with `~ 1` (or just `~`) includes an intercept (baseline
      group).
2.  **Additive Factors:**
    - Factors separated by `+` indicate additive effects. For example,
      `~ 0 + factor1 + factor2` means you are modeling the effects of
      `factor1` and `factor2` additively, without considering
      interactions.
3.  **Interaction Terms:**
    - The `*` symbol is used to model interactions between factors. For
      example, `~ 0 + factor1 * factor2` will include `factor1`,
      `factor2`, and their interaction (`factor1:factor2`).
    - Alternatively, you can specify the interaction explicitly with
      `:`. For example, `~ 0 + factor1 + factor2 + factor1:factor2` is
      equivalent to `~ 0 + factor1 * factor2`.

### Some examples:

- `~ 0 + factor1 + factor2`: Additive model without an intercept.
- `~ 1 + factor1 + factor2`: Additive model with an intercept.
- `~ 0 + factor1 * factor2`: Model with main effects and their
  interaction, no intercept.
- `~ 1 + factor1 * factor2`: Model with intercept, main effects, and
  their interaction.

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
