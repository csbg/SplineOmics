# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Template for Changelog Entries
Each release section should follow the structure below:

### [Version] - YYYY-MM-DD
#### Added
- New features or functionality introduced in this release.

#### Changed
- Updates or modifications to existing features.

#### Fixed
- Bugs that have been resolved.

#### Deprecated
- Features that are still functional but are slated for removal in the future.

#### Removed
- Features or functionality that have been removed.

#### Security
- Any security-related updates or patches.

Examples:
- **Added**: Introduced a new plotting function `plotTimeSeries()`.
- **Fixed**: Resolved an issue causing crashes when input data was missing values.
- **Changed**: Modified the default parameters for `normalizeData()`.

---

### [0.1.1] - 2025-MM-DD

### Changed
- The design formula must now contain the string 'Time' rather than 'X' like it was before. X from
  before stood for the time. This change is intended to make the design formula more explicit and 
  self explanatory.
- Random effects can now be directly be specified in the design formula, rather
  than being passed as part of the dream_params.
  
### Added
- Added linear mixed models for modeling random effects. The variancePartition
  package is used for that --> voomWithDreamWeights for RNA-seq data processing,
  and dream as the replacement for limma. For example, if you Reactor is your
  random effect, you can write the design formula like this:
  design <- "~ 1 + Condition*Time + Plate + (1|Reactor)" and SplineOmics will
  automatically run the variancePartition functions voomWithDreamWeights() and
  dream() instead of the limma::voom and lmfit. dream() has additional
  parameters, such as the method and the degree of freedom (different than the
  degree of freedom used for the splines in this package) and you can pass these
  with the dream_params argument. See RNA-seq analysis vignette or the
  respective function references for more info.
- Raw data plotting function --> make_scatter_plot_html()  --> see references.
- Imputed values are marked with triangle symbols in cluster_hits() spline
  plots if raw data is passed.
- Package version is written in each HTML report (the tag) and the session info
  is added as an embedded text file.
- Standard error written on top of "double spline" plots (limma result category
  2 and 3). 
- Used analysis script is embedded as a text file in the reports.
- The mode (integrated or isolated) is written on top of the reports in a 
  separate field.
- Now there are 4 significance stars.
- Treatment lines for the double spline plots (limma result category 2 and 3).

#### Fixed
- A few smaller visual things for the HTML reports.