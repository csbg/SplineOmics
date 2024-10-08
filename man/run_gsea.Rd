% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_gsea.R
\name{run_gsea}
\alias{run_gsea}
\title{Generate a GSEA Report}
\usage{
run_gsea(
  levels_clustered_hits,
  databases,
  report_info,
  clusterProfiler_params = NA,
  plot_titles = NA,
  universe = NULL,
  report_dir = here::here()
)
}
\arguments{
\item{levels_clustered_hits}{A list of clustered hits at different levels.}

\item{databases}{A list of databases for the gene set enrichment analysis.}

\item{report_info}{A list containing information for the report generation.}

\item{clusterProfiler_params}{Additional parameters for the GSEA analysis,
default is NA. Those include adj_p_value,
pAdjustMethod, etc (see clusterProfiler
documentation).}

\item{plot_titles}{Titles for the plots, default is NA.}

\item{universe}{Enrichment background data, default is NULL.}

\item{report_dir}{Directory where the report will be saved, default is
`here::here()`.}
}
\value{
A list of plots generated for the GSEA report.
}
\description{
This function generates a Gene Set Enrichment Analysis (GSEA) report based
on clustered hit levels, gene data, and specified databases. It processes
the input data, manages GSEA levels, and produces an HTML report with plots.
}
