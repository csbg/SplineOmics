% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/screen_limma_hyperparams.R
\name{get_limma_combos_results}
\alias{get_limma_combos_results}
\title{Generate LIMMA Combination Results}
\usage{
get_limma_combos_results(
  datas,
  rna_seq_datas,
  metas,
  designs,
  modes,
  condition,
  spline_test_configs,
  feature_names,
  adj_pthresholds,
  padjust_method
)
}
\arguments{
\item{datas}{A list of matrices.}

\item{rna_seq_datas}{A list of RNA-seq data objects, such as the voom object
derived from the limma::voom function.}

\item{metas}{A list of metadata corresponding to the data matrices.}

\item{designs}{A list of design matrices.}

\item{modes}{A character vector containing 'isolated' or 'integrated'.}

\item{condition}{A single character string specifying the condition.}

\item{spline_test_configs}{A configuration object for spline tests.}

\item{feature_names}{A character vector of feature names.}

\item{adj_pthresholds}{A numeric vector with elements > 0 and < 1.}

\item{padjust_method}{A single character string specifying the p-adjustment
method.}
}
\value{
A list of results for each combination of data, design, and spline
configuration.
}
\description{
Computes results for various combinations of data, design matrices, and
spline configurations using the LIMMA method.
}
