% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_limma_report.R
\name{generate_avrg_diff_plots}
\alias{generate_avrg_diff_plots}
\title{Generate Plots for Average Difference Conditions}
\usage{
generate_avrg_diff_plots(avrg_diff_conditions, adj_pthresh)
}
\arguments{
\item{avrg_diff_conditions}{A list of top tables from the LIMMA analysis
representing the average difference between conditions.}

\item{adj_pthresh}{A numeric value specifying the adjusted p-value threshold
for significance.}
}
\value{
A list containing the plots and their sizes, as well as the
section header information.
}
\description{
Creates p-value histograms and volcano plots for each condition in the
average difference conditions. This function is used internally in the
`create_limma_report` function.
}
