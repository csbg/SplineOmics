% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_general.R
\name{check_homoscedasticity_violation}
\alias{check_homoscedasticity_violation}
\title{Check for Violation of Homoscedasticity in Linear Model Inputs}
\usage{
check_homoscedasticity_violation(
  data,
  meta,
  condition,
  data_type = "other-omics",
  p_threshold = 0.05
)
}
\arguments{
\item{data}{A numeric matrix of expression values (rows = features, 
columns = samples). Should already be log-transformed or 
otherwise variance-stabilized if required.}

\item{meta}{A data frame of sample metadata. Must have one row per column 
in `data`.}

\item{condition}{A character string indicating the name of the column in
`meta` representing the condition of interest.}

\item{data_type}{String specifying if rna-seq data is passed, or other omics
data. Based on this, the message displayed in case of a 
significant violation of homoscedasticity informing about
the strategy to mitigate the issue will change (different
strategy for rna-seq and other omics data).}

\item{p_threshold}{Numeric. Significance threshold for the Wilcoxon test
(default is 0.05).}
}
\value{
A logical value indicating whether a statistically significant 
        difference in variance was detected (`TRUE` = assumption violated,
        `FALSE` = no violation). A summary of the test result is also printed
        to the console.
}
\description{
This internal helper function tests whether the assumption of 
homoscedasticity (equal variances) is violated across two levels of a given
experimental condition. It performs a paired Wilcoxon signed-rank test on 
per-feature (e.g., gene or protein) sample variances across the two groups.
If the test is significant, it suggests that variance differs systematically 
between conditions, which may bias linear model fits.
}
