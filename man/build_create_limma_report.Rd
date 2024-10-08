% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_limma_report.R
\name{build_create_limma_report}
\alias{build_create_limma_report}
\title{Build Cluster Hits Report}
\usage{
build_create_limma_report(
  header_section,
  plots,
  plots_sizes,
  level_headers_info,
  report_info,
  output_file_path = here::here()
)
}
\arguments{
\item{header_section}{A character string containing the HTML header section.}

\item{plots}{A list of ggplot2 plot objects.}

\item{plots_sizes}{A list of integers specifying the size of each plot.}

\item{level_headers_info}{A list of header information for each level.}

\item{report_info}{A named list containg the report info fields. Here used
for the email hotkey functionality.}

\item{output_file_path}{A character string specifying the path to save the
HTML report.}
}
\value{
No return value, called for side effects.
}
\description{
Generates an HTML report for clustered hits, including plots and
spline parameter details, with a table of contents.
}
\seealso{
\code{\link{plot2base64}}, \code{\link{create_progress_bar}}
}
