# FAQ

## Frequently Asked Questions

------------------------------------------------------------------------

### FAQ on data pre-processing

#### Does SplineOmics handle missing values?

Yes. SplineOmics can handle missing values gracefully because it uses
`limma` as the statistical engine for model fitting and differential
analysis. The `limma` package is designed to work with incomplete
expression matrices — missing values are simply ignored when estimating
model parameters and computing contrasts, without causing errors.

This means you do not necessarily need to impute or remove features with
a few missing entries before running SplineOmics. However, if a feature
has too many missing values (for example, across most samples), consider
removing it.

------------------------------------------------------------------------

### FAQ on parameter selection

#### Which design formula should I use?

The choice of design formula depends on your experimental setup and on
whether you want to model your conditions independently or jointly.

- Isolated designs (no interaction between condition and time)  
  Use this when you want to analyse each condition completely
  separately, for example, to get independent spline fits for each
  treatment.  
  In this case, `SplineOmics` will simply process each dataset
  sequentially for convenience, so you do not need to run the code
  twice.  
  This approach is ideal when you do not expect conditions to share
  information or influence one another.

- Integrated designs (conditions modelled jointly with a time
  interaction)  
  Use this when you want to fit all conditions in a single model and
  allow their time trends to differ through an interaction term between
  condition and time.  
  With this approach, datasets are modelled jointly and can borrow
  statistical power from each other, leading to more stable variance
  estimates and more sensitive detection of shared temporal patterns.  
  This also enables access to the full range of `limma` result
  categories, including category 2 and category 3 results (see the
  vignette on [limma result
  categories](https://csbg.github.io/SplineOmics/articles/limma_result_categories.md)).

#### Should I use B-splines or natural cubic splines?

Both spline types can model smooth trends across time or another
continuous variable, but they differ in how local their flexibility is.

- B-splines are locally adaptive: changing one knot affects the fitted
  curve only in a small neighbourhood around that knot.  
  This makes them ideal when you expect local variations (for example,
  short-term biological responses) and want the rest of the curve to
  remain stable. The trade-off is that B-splines typically use more
  degrees of freedom, so the model can become more complex.

- Natural cubic splines enforce global smoothness: each basis function
  extends across the entire range, so adjusting one part of the curve
  slightly influences the whole fit. They use fewer degrees of freedom
  and can be more stable when you expect overall smooth behaviour, but
  they are less suitable if local detail is important.

In short:  
\> Choose B-splines when local flexibility matters, and natural cubic \>
splines when you want a smoother, more global trend.

#### How many degrees of freedom should I use for the splines?

There is no single optimal choice — the best number of degrees of
freedom (dof) depends on the smoothness and complexity of your data.  
In practice, 2 or 3 degrees of freedom work well in most cases.

Using more degrees of freedom makes the spline wigglier and increases
the risk of overfitting, while using only 1 degree of freedom is usually
too restrictive to capture meaningful trends.

If you set the degrees of freedom to 0, `SplineOmics` will automatically
determine the optimal value using leave-one-out cross-validation,
selecting the dof that provides the best predictive performance for your
dataset.

#### Should I use array weights?

Yes, in most cases this is recommended. Since `SplineOmics` builds on
`limma`, it inherits support for array weights, which help to correct
for heteroskedasticity (unequal variances) often present in time-series
data.

By estimating a quality weight for each sample, `SplineOmics` can give
less influence to noisier samples and more weight to consistent ones,
resulting in smoother spline fits and higher statistical power.

#### How many clusters should I use?

There is no universal rule for choosing the optimal number of
clusters.  
Using more clusters makes each cluster purer — its centroid represents
its members more precisely — but it also makes downstream interpretation
more complex and fragmented.

To help with this choice, `SplineOmics` provides several aids:

- For each cluster, it reports the variance explained by the cluster
  centroid, both as a mean value and as a distribution histogram, with
  guidance on how to interpret these values.  

- For each individual feature, it shows how well it is represented by
  its cluster centroid.  

- You can specify a range or custom set of cluster numbers (for example,
  2–10 or {2, 5, 6}), and `SplineOmics` automatically selects the best
  one using the Bayesian Information Criterion (BIC):

  ``` math

  \text{BIC} = n_{\text{obs}} 
  \log\left(\frac{\text{tot\_within}}{n_{\text{obs}}}\right) 
  + k \log(n_{\text{obs}}) \times p
  ```

The model with the lowest BIC is chosen as the optimal clustering
configuration.

------------------------------------------------------------------------

### FAQ on interpreting the results

#### How is FDR handled in `run_limma_splines()` when conditions are modelled

jointly (\`mode = “integrated”)?

The joint spline model is fitted once across all conditions, but FDR is
controlled separately for each logical family of tests. There are three
types of results, each with its own FDR correction:

1.  Per-condition time effects (`extract_within_level_time_effects()`):
    For each condition level, spline (and interaction) contrasts are
    tested using an F-test (or t-test if `dof = 1`). P-values are
    adjusted across all genes within that condition. FDR is controlled
    per condition-specific time effect, not across conditions.

2.  Pairwise average differences between conditions (`condition_only`
    results): For each condition pair, the main condition effect plus
    its time-related terms are tested jointly. P-values are adjusted
    across all genes for that pairwise contrast. FDR is controlled per
    average-difference contrast per pair, not across all pairs.

3.  Pairwise condition×time interactions (`condition_time` results): For
    each condition pair, only the spline interaction terms are tested,
    via F- or t-tests depending on the degrees of freedom. P-values are
    adjusted across all genes for that interaction contrast. FDR is
    controlled per time-interaction contrast per pair, not jointly with
    average-difference contrasts or other pairs.

#### How is FDR handled in `run_ora()`?

FDR is controlled in a structured way across three levels of testing:

1.  Limma spline results (highest level):  
    Each limma spline result (e.g., time effect for one condition,
    interaction effects, etc.) represents a separate biological
    question.  
    Therefore, FDR is not shared or adjusted across different limma
    results.

2.  Clusters within a limma result:  
    A limma result produces a set of hit genes that are clustered.  
    All clusters belong to the same biological comparison, so they form
    a shared multiple-testing family.

3.  Databases within a limma result:  
    Each gene-set database (GO BP, KEGG, Reactome, etc.) is analyzed
    independently because they differ in size, structure, and hypothesis
    space.  
    Therefore, FDR is controlled per database, not across databases.

run_ora() retains all tested terms (no internal filtering), collects all
p-values for each database across all clusters of a given limma result,
and applies a single multiple-testing correction within that database.
The resulting adjusted p-values are stored in `p_adj_by_db` and are used
for all downstream filtering and visualization.

#### How do I interpret category 2 (avrg diff between conditions) ORA results?

Category 2 evaluates average differences between conditions after
accounting  
for time. Each pairwise contrast creates one
`cluster_cat2_<cond1>_vs_<cond2>`  
column in the cluster table.

For each feature, this column contains:

- `"cond1_higher"` — significantly higher on average in cond1  
- `"cond2_higher"` — significantly higher on average in cond2  
- `NA` — not significant

These labels define two gene sets: genes higher in cond1, and genes
higher in  
cond2. SplineOmics runs over-representation analysis separately for
these  
two sets, using all modeled features as the background universe.

Thus, Category 2 enrichment answers:

- Which biological processes are higher in condition A than B?  
- Which processes are higher in condition B than A?

#### How do I interpret Category 3 (condition–time interaction) ORA results?

Category 3 evaluates whether the time trajectories differ between two
conditions. For each pairwise contrast, SplineOmics creates one
`cluster_cat3_<cond1>_vs_<cond2>` column in the cluster table.

For each feature, this column contains:

- `"<k1>_<k2>"` — the combo cluster, built from the time-effect cluster
  of the feature in `cond1` (`k1`) and in `cond2` (`k2`)
- `"ns"` replaces a missing cluster if the feature is not significant in
  that condition (e.g. `"2_ns"`)
- `NA` — the feature is not significant for the interaction contrast

Each unique combo cluster (e.g. `2_3`, `1_4`, `2_ns`) defines a gene
set. SplineOmics performs over-representation analysis for each such
set, using all modeled features as the background universe.

Thus, Category 3 enrichment answers questions such as:

- Which biological processes are enriched among genes that rise in
  condition A but fall in condition B?
- Which processes are specific to genes that change only in one
  condition?

This highlights biological pathways whose time-course behavior differs
between the compared conditions.

------------------------------------------------------------------------

### General FAQ

#### Can SplineOmics handle datasets with more than two conditions?

Not directly. `SplineOmics` is designed for pairwise comparisons between
conditions. You can, however, analyse experiments with more than two
conditions by either:

- using an isolated design, where each condition is analysed separately,
  or  
- performing all pairwise comparisons among the conditions of interest.

These approaches allow you to explore multiple conditions, but they do
not integrate them into a single joint model. Integrated designs are
currently limited to comparisons between two conditions at a time.

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
