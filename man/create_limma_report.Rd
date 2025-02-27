% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_limma_report.R
\name{create_limma_report}
\alias{create_limma_report}
\title{create_limma_report()}
\usage{
create_limma_report(splineomics, adj_pthresh = 0.05, report_dir = here::here())
}
\arguments{
\item{splineomics}{An S3 object of class `SplineOmics` that contains all the
necessary data and parameters for the analysis, including:
\itemize{
  \item \code{limma_splines_result}: A list containing top tables from
   differential expression analysis for the three different limma results.
  \item \code{meta}: A data frame with sample metadata. Must contain a column
                     "Time".
  \item \code{condition}: A character string specifying the column name in
                         the metadata (\code{meta}) that defines groups
                         for analysis. This column contains levels such as
                          "exponential" and "stationary" for phases, or
                          "drug" and "no_drug" for treatments.
  \item \code{annotation}: A data frame containing feature information,
                           such as gene and protein names, associated with
                           the expression data.
  \item \code{report_info}: A list containing metadata about the analysis
                            for reporting purposes.
}}

\item{adj_pthresh}{A numeric value specifying the adjusted p-value threshold
for significance. Default is 0.05. Must be > 0 and < 1.}

\item{report_dir}{A string specifying the directory where the report should
be saved. Default is the current working directory.}
}
\value{
A list of plots included in the generated HTML report.
}
\description{
Generates an HTML report based on the results of a limma analysis with
splines.
The report includes various plots and sections summarizing the analysis
results for time effects, average differences between conditions,
and interaction effects between condition and time.
}
