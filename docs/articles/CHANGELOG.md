# Changelog

## Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a
Changelog](https://keepachangelog.com/en/1.0.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

### Template for Changelog Entries

Each release section should follow the structure below:

#### \[Version\] - YYYY-MM-DD

##### Added

- New features or functionality introduced in this release.

##### Changed

- Updates or modifications to existing features.

##### Fixed

- Bugs that have been resolved.

##### Deprecated

- Features that are still functional but are slated for removal in the
  future.

##### Removed

- Features or functionality that have been removed.

##### Security

- Any security-related updates or patches.

------------------------------------------------------------------------

#### \[v0.4.3\] (in preparation)

##### Added

- condition_colours support to plot options. Users can now provide a
  named list mapping condition levels to colours. When supplied, these
  colours override the default condition colours in plots, enabling
  consistent and user-defined visual styling across figures.

- min_cluster_r2 pruning with explicit ‘other’ cluster. Introduce a  
  min_cluster_r2 parameter to enforce a minimum r²-to-centroid
  constraint after k-means clustering.

##### Changed

- Quality control reporting and modality normalization in
  cluster_genes_multiomics().

- Separated statistical computations and plot/HTML-report generation
  for: cluster_hits(), run_ora(), and find_pvc(). Therefore, the new
  functions: create_clustering_report(), create_ora_report(), and
  create_pvc_report(); are now containing the plotting/HTML-report
  generation code for the mentioned statistical functions.

#### \[v0.4.2\]

##### Added

- [`cluster_genes_multiomics()`](https://csbg.github.io/SplineOmics/reference/cluster_genes_multiomics.md),
  a gene-centric multi-omics clustering function that integrates
  multiple data modalities within analytical blocks via weighted
  gene–gene distance aggregation. The function returns both per-gene
  cluster assignments and modality-specific cluster centroid summaries,
  enabling interpretation of cluster structure, coverage, and
  within-cluster coherence across omics layers.

- [`run_ora()`](https://csbg.github.io/SplineOmics/reference/run_ora.md):
  Gene names additional safety checks

##### Fixed

- [`cluster_hits()`](https://csbg.github.io/SplineOmics/reference/cluster_hits.md):
  cT and cDT results sheet bug

- 1 condition bug in
  [`explore_data()`](https://csbg.github.io/SplineOmics/reference/explore_data.md).

- `extract_effect()` internal function, made it much more robustly (had
  clear weaknesses before)

- Part of the documentation

#### \[v0.4.1\]

##### Added

- Effect size thresholds are now reported in the HTML report generated
  by the cluster_hits() function.

- Double spline plots now better inform about non significant (ns) and
  not relevant (nr) results, and filter better based on this.

- cluster_genes_multiomics(): Gene-centric multi-omics clustering
  function that integrates multiple layers per block, supports
  many-to-one mappings, computes weighted distance matrices, assigns
  gene clusters, and returns centroid-level summaries for visualization.
  However, this is still a beta addition, so it might not work as
  intendet! Also, no vignette was added yet. This will be done in a
  future version.

##### Fixed

- use_array_weights behavior for variancePartition linear mixed models.

- Non-significant category 2 and 3 results plots shown in HTML report.
  This was now fixed (now they are not shown anymore).

- Multi condition bug in the contrasts in run_limma_splines(). For the
  average diff condition and interaction condition time results, it did
  not correctly create the intended pair contrasts. Now, the code
  overall was refactored, and this issue fixed in this refactoring. Now,
  a joint contrast matrix L is created before the model is fit, and is
  used in the fit for variancePartion, and in the limma case after model
  fitting with contrast.fit.

#### \[v0.4.0\]

##### Changed

- SplineOmics now supports any amount of conditions (values in the meta
  condition column).

- adj_pthresholds argument of cluster_hits() renamed to
  adj_pthresh_time_effect. Now, it instead requires a single numeric
  input specifying the adj. p-value treshold for all conditions (e.g.
  0.05), and does not support anymore setting a separate value for each
  condition. This is because using separate values break the uniformity
  of inference and introduces uncontrolled type I error inflation.

- plot_info argument of cluster_hits() now does not accept anymore the
  value ‘double_spline_plots’ in the treatment_labels and
  treatment_timepoints lists. Instead, the treatment lines specified for
  the conditions will automatically be applied to all pairwise
  comparisons where this condition appears.

- Updated ORA workflow to disable internal filtering, apply
  database-wise multiple-testing correction across all clusters via a
  new `p_adj_by_db` column, and use this value for all downstream
  significance filtering to ensure consistent FDR control.

#### Removed

- qvalueCutoff value in clusterProfiler_params argument of run_ora().
  The reason is that the qvalueCutoff is no longer statistically
  meaningfull with the new way the raw p-values are adjusted (see
  above).

#### \[v0.3.4\]

##### Changed

- Renamed the function download_bioconductor_database() to
  extract_gene_sets().

- Stylistic changes to cluster_hits HTML report, such as font sizes for
  plots.

#### \[v0.3.3\]

##### Fixed

- Bugs occuring in the functions extract_data(), explore_data(), and
  cluster_hits()

##### Removed

- Removed the report argument of the cluster_hits() and run_ora()
  function. Instead, now, the report_dir argument is per default NULL,
  which means that no HTML report is generated. Only when a directory
  path is provided with the report_dir argument, a HTML report is
  generated.

##### Added

- The signed r² as the quality score for the cluster members is now also
  written in the downloadable topTables inside the cluster_hits() HTML
  report.

#### \[v0.3.2\]

##### Fixed

- The cluster_hits() function now does not modify the gene symbols
  anymore. Before, it attempted to clean them by removing characters
  that are not valid and making them uppercase. This is implicit
  behavior (even though the user was informed) and was removed in the
  favor of the user now having to explicitly cleaning them himself and
  he is responsible for that.

##### Changed

- Changed cluster quality metric from Silhouette score to signed r^2
  (variance explained). This is because Silhouette score focuses on
  determining how well a cluster member fits in this cluster compared to
  other clusters, but does not quantify how well the cluster centroid
  represents each member (which is what is desired).

- Identifying the best amount of clusters in the case the user specified
  a range of clusters is now again done with the Bayesian Information
  Criterion. This is because in order to use the Silhouette score, a new
  dependency is needed, and one aim is to minimize dependencies.

##### Removed

- Removed the double spline clustering and cat3 now forms the clusters
  by combining the cluster assignments of the time effect clusters. This
  was done to simplify the task, because one now needs to only focus on
  nicely clustering the time effect hits.

- Gene check inside run_ora(). Now, the user is responsible of providing
  appropriate symbols.

- report argument of the cluster_hits() and run_ora() functions. This
  was a Boolean argument specifying if a report should be generated.
  Now, this was removed, and the functions only generate a report when a
  directory path to write was provided. Now also the default for the
  directory path is NULL, instead of the project dir. This means that
  per default, no HMTL report is generated.

#### \[v0.3.1\]

##### Added

- Effect size computations and thresholding for the splines. Because
  with enough statistical power, also tiny differences can be
  significant, it can make sense to determine the hits by not just
  filtering with the adj.p-value threshold, but also with an effect size
  threshold. Therefore, the cluster_hits() function now computes the
  spline effect sizes via cumulative travel and accepts a threshold
  argument that is used for the cutoffs. For category 2 with the average
  difference between the conditions, the effect size is already in the
  topTable from limma.

- Double spline clustering for the hits for category 3. Those clusters
  are assigned to those hits and serve for the following clustered
  enrichment. To control this, a new argument of the cluster_hits()
  function, nr_clusters_interaction, was defined.

- Cluster quality metric (Silhouette score). Now, in the HTML reports of
  the cluster_hits() function, for every hit, there is the Silhouette
  score reported. Additionally, in the beginning of the time effect
  results sections, histograms show the distribution of Silhouette
  scores in all clusters.

- The functions preprocess_rna_seq_data, run_limma_splines, and
  cluster_hits now report the amount of time it took for them to run.

##### Fixed

- Some documentation was outdated or too vague. Improved that.

- Some smaller things were fixed, such as remnants of the
  screen_limma_hyperparams() function in the code in other places, or
  input control functions that still expected older argument names.

#### \[0.3.0\]

##### Changed

- cluster_hits() returns. Now it returns only the plots and a large
  table, summarizing all clustering results, including the results from
  the limma result categories 2 and 3.

- run_ora() now takes the cluster_hits() summary table return as input,
  and performs enrichment analysis with all three limma result
  categories.

- More info in run_ora() reports: The HTML reports generated by the
  run_ora() function contained too little information about the inputs
  and the results to which the ora report is connected. To adress this,
  the ora report now contains downloadable files with the foreground and
  background genes, reports the clusterProfiler parameters used, and the
  run_ora() function has a new argument which allows to pass a string
  which is the name of the corresponding cluster_hits() HTML report,
  that contains the results that are the basis for the respective
  run_ora() report.

- Separated isolated and integrated mode computations: Before, there was
  a part in the code that generated the limma time_effect category
  results. This was used both for the isolated and integrated mode,
  because they share those results. However, this meant that the
  analysis for the integrated mode was run more often then necessary,
  because it was run once for each condition to the the time_effect
  results, and then once to get the average diff condition and
  interaction results. This was unnecessary, since all the integrated
  results can be derived with a single global model (with the help of
  contrasts). Therefore, I now completely separated the compuations for
  the isolated and integrated approach, and for the integrated approach,
  just one global model is run and all the three limma result categories
  extracted from there with the help of contrasts.

- Made cluster_hits() usable for many hits: Introduced a new function
  argument, called max_hit_number. This controls the plotting logic of
  the report generation, so it limits the amount of plots and amount of
  detail shown in the plots to a maximum number of hits per cluster.
  This allows to also generate clustering reports for analyses with many
  hits, because nowbody anyways looks at thousands of spline plots.

- Switched from hierarchical clustering to kmeans. This was done because
  kmeans is much faster (more efficient in terms of computations). Also,
  now the number of clusters is always specified as a range (min:max,
  e.g. 2:4) and all of those cluster numbers are tried out and the best
  one identified with the help of the silhuette score. Also if you want
  to have just one definite amount of clusters, you still specify it as
  a range due to technical reasons. Additionally, when the datasets is
  large, it uses the minibatch version of kmeans, which is a lot faster
  (while sacrificing a bit of the cluster quality).

##### Added

- Auto-dof with loocv: Introduced new capability, that automatically
  finds the best spline dof based on the error on loocv. For every
  linear model fit, it performs the analytical loocv and then looks at
  the error sum of all features. It does this with every dof from 2 to
  min(timepoints, 10). Autodof is activated when the the dof = 0L is
  selected.

- Added cross-species gene mapping: Added optional support for automatic
  gene symbol conversion between species via gprofiler2 or orthogene.
  Controlled through the new mapping_cfg argument in run_ora(), this
  enables mapping from from_species to to_species using either method.
  The conversion applies across all cluster hit tables, preserves order,
  and retains unmapped genes. No suffix stripping is done, but users are
  warned if trailing \_digits are detected. The qvalue package is now
  optionally required to enable stricter FDR filtering in enrichment
  results.

- Cross-platform and memory-safe parallelization support for dream() via
  a new helper function bp_setup(). This function sets up a BiocParallel
  backend using either MulticoreParam (Unix) or SnowParam (Windows),
  based on user-supplied configuration (bp_cfg). BLAS threads per worker
  are throttled using RhpcBLASctl to avoid oversubscription. Cleanup
  with bpstop() is conditionally applied for SnowParam backends.
  Includes validation of bp_cfg input. Enables efficient multi-core
  execution of dream() and eBayes() while preventing system resource
  exhaustion.

- download_bioconductor_database() function. This one allows to
  conveniently download a database from a specific organism from
  BioConductor (see function documentation).

- Introduced clusterProfiler::enrichGO() along with
  clusterProfiler::simplify into SplineOmics::run_ora(). This allows to
  run the enrichment of the BioConductor GO databases with the enrichGO
  function, followed by the simplify function, to collapse hierarchical
  terms. This behavior is controlled over a new argument of run_ora(),
  namely enrichGO_cfg (see documentation).

- Better documentation in the limma and cluster_hits HTML reports about
  the heteroscedasticity arguments and the result of the Levene test.

##### Fixed

- Red spline curves: Before, it was just generated with the help of the
  spline coefficients, and the intercept. But this is not the full
  linear model that is used by limma. Now, it generates the spline also
  with the help of the other variables used in the linear model. The
  predicted splines are now generated at the very beginning of the
  cluster_hits() function, and used for clustering AND plotting.

- Wrong coefficients passes to topTable() function. For the limma result
  category average difference conditions, not all necessary coefficients
  were passed. Before, it just checked if there is a difference between
  the conditions for the reference time point (the first timepoint). Now
  it passes all the necessary coefficients, and then performs a joint
  F-test to get the average difference between the conditions across the
  whole timecourse.

- The use_array_weights argument internally was always set to FALSE when
  it was set to TRUE by the user. This is why setting it to TRUE
  seemingly had no effect. This was now fixed.

------------------------------------------------------------------------

#### \[0.2.0\] - 2025-05-08

##### Fixed

- Bug that prevented the generation of the double spline plots when the
  condition column of meta was a factor instead of a string. This is now
  solved (it internally in the code converts it to a string), so having
  the condition column as a factor now does not cause problems anymore
  (as having it as factor is perfectly valid).

- Bug that occured in the internal clean_gene_symbols() function when
  one of the gene names was NA. Now, this does not cause a problem
  anymore and instead of the gene name on the respective spline plot it
  will write NA.

- The error message when more clusters where specified than hits in a
  condition was very cryptic. Now, it gives a message that is easy to
  understand.

- A bug in input control that threw an error indicating that there are
  rows in the dataset that contains just zeros, even though this was not
  true.

- When running cluster_hits() with report = FALSE, but there where many
  hits, it anyways warned the user that the report generation will take
  long. This was fixed, because when report = FALSE, no report is
  generated, which means the user does not have to be warned then (that
  would just confuse the user).

##### Changed

- Replicates in the double spline plots (limma result category 2 and

  3.  are now shown with shapes instead of numbers besides the data
      points. This is now possible because the imputed values are not
      anymore shown with shapes in those plots, like in the single
      spline plots (limma result category 1), but with colors (red =
      imputed values for condition 1, dodgerblue = imputed values for
      condition 2).

- The images in the HTML reports are now considerably smaller, but
  zoomable, with the option to drag within the zoomed image. A single
  left-click with the mouse starts the zoom process. Then, you can zoom
  with the mouse wheel. If you reach the max zoom, you can move in the
  image by further zooming with the mouse wheel but placing the mouse in
  one of the corners. Further, one can zoom by holding left-click
  pressed and dragging (you can then also let go).

- SplineOmics can now handle datasets with NA values. limma can do that,
  and since SplineOmics is based on limma, it can in principle do it to.
  The issue before was just that there were checks in place that
  prevented NA values from entering the pipeline, and also other steps
  in the code caused problems, raised warnings, or stopped the execution
  when facing NA values.

- Changed extract_data() function so that it can handle missing values
  correctly, and so that it now is fully explicit and does not do any
  internal magic.

- When there where less than two hits in all levels for the
  cluster_hits() function, then this function provided a message that
  informed about this and threw an error. This was fine when just
  running the function to generate one result. However, when using this
  function in an outside loop to generate multiple results, this
  behavior is suboptimal, because it then must be captured with a
  tryCatch clause. Now, if there are less than two hits in all levels,
  it simply returns NULL, which is a cleaner behavior that does not
  break loops if one result cannot be generated.

- extract_data() function, so that now it does not work with an implicit
  logic, but instead with an explicit logic. Before, it tried to
  automatically find the numeric matrix field. Now, it requires to
  specific the top- and bottom row, and the left and right column. This
  then marks the field where the matrix is. The rest of the logic (the
  logic about the row name columns) is still the same.

- make_scatter_plot_html() (now renamed to make_scatter_plots_html()).

- renamed run_gsea() to run_ora(), because the function uses
  clusterProfiler, which performs ora (overrepresentation analysis) and
  not gsea (gene set enrichment analysis) in the stricter sense.

#### Added

- find_pvc() function. This function performs a compount contrast for
  every timepoint triple in the data (adjacent timepoints) to find
  peaks, valley, and cliffs as a form of local temporal patterns. It
  generates an HTML report containing all findings in files and also
  containing the plots.

- New vignette showcasing the find_pvc() function.

- Implicit (default) and explicit heteroscedasticity handling of the
  linear model input. Controlled by the `robust_fit` Boolean flag
  indicating whether the robust modeling strategy should be used to
  account for heteroscedasticity. If omitted (i.e.,
  `robust_fit = NULL`), the decision is made implicitly based on the
  result of the Breusch–Pagan test, which is applied independently to
  each feature. This test fits a linear model of the form
  `expression ~ time` and checks whether the residual variance
  systematically depends on the fitted values. If a significant
  violation of homoscedasticity is detected in a sufficient fraction of
  features (default: 10%), the robust strategy is applied. For RNA-seq
  data, this means using `voomWithQualityWeights()` instead of the
  standard `voom()` function. For other, non-count-based data,
  [`limma::arrayWeights()`](https://rdrr.io/pkg/limma/man/arrayWeights.html)
  is used during model fitting, combined with `robust = TRUE` in
  [`limma::eBayes()`](https://rdrr.io/pkg/limma/man/ebayes.html). These
  strategies downweight samples with higher residual variance to prevent
  bias in the linear modeling step, as violations of the
  homoscedasticity assumption can lead to misleading p-values and
  unreliable inference.

- New note in the HTML report of the cluster_hits() function, that
  informs the user that the splines shown in the plots can appear to
  have the wrong intercept. This can occur when a batch effect and/or
  random effect is modeled with limma or the linear mixed models, but
  the plotting data is batch corrected only with the dedicated limma
  batch correction function. For a reason, there is a gap. The results
  are fine! Just the plotting is off!

- New user-available function called compare_results() that allows to
  correlate the topTable results from two SplineOmics results (for
  example intergrated vs. isolated approach with the same data).

------------------------------------------------------------------------

#### \[0.1.2\] - 2025-02-11

##### Added

- Level name in title of overall spline plots.

##### Fixed

- Small bug that made it impossible having no treatment label for a
  condition.
- Small bug that did not allow to specify two treatment labels.
- mode == isolated for RNA-seq data. Before, it caused an error, because
  it splits up the meta into the different conditions, but does not do
  so for the data. Now, it informs the user that it cannot do this, and
  the user has to split up the data himself outside (just for RNA-seq
  data).

##### Removed

- open_template() function. This was a function that opened a template
  for your own analysis. I removed it, because otherwise, whenever I
  change the code, I have to update the get-started vignette (tutorial)
  and this, which is twice the maintainance effort.

------------------------------------------------------------------------

#### \[0.1.1\] - 2025-01-29

#### Changed

- The design formula must now contain the string ‘Time’ rather than ‘X’
  like it was before. X from before stood for the time. This change is
  intended to make the design formula more explicit and self
  explanatory.
- Random effects can now be directly be specified in the design formula,
  rather than being passed as part of the dream_params.

#### Added

- Added linear mixed models for modeling random effects. The
  variancePartition package is used for that –\> voomWithDreamWeights
  for RNA-seq data processing, and dream as the replacement for limma.
  For example, if you Reactor is your random effect, you can write the
  design formula like this: design \<- “~ 1 + Condition\*Time + Plate +
  (1\|Reactor)” and SplineOmics will automatically run the
  variancePartition functions voomWithDreamWeights() and dream() instead
  of the limma::voom and lmfit. dream() has additional parameters, such
  as the method and the degree of freedom (different than the degree of
  freedom used for the splines in this package) and you can pass these
  with the dream_params argument. See RNA-seq analysis vignette or the
  respective function references for more info.
- Raw data plotting function –\> make_scatter_plot_html() –\> see
  references.
- Imputed values are marked with triangle symbols in cluster_hits()
  spline plots if raw data is passed.
- Package version is written in each HTML report (the tag) and the
  session info is added as an embedded text file.
- Standard error written on top of “double spline” plots (limma result
  category 2 and 3).
- Used analysis script is embedded as a text file in the reports.
- The mode (integrated or isolated) is written on top of the reports in
  a separate field.
- Now there are 4 significance stars.
- Treatment lines for the double spline plots (limma result category 2
  and 3).

##### Fixed

- A few smaller visual things for the HTML reports.

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
