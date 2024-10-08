% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_report_generation.R
\name{create_enrichr_zip}
\alias{create_enrichr_zip}
\title{Create a ZIP File for Enrichr Gene Lists}
\usage{
create_enrichr_zip(enrichr_format)
}
\arguments{
\item{enrichr_format}{A list with the formatted gene lists and background
gene list, typically the output of `prepare_gene_lists_for_enrichr`.}
}
\value{
A base64-encoded string representing the ZIP file.
}
\description{
This function creates a ZIP file containing directories for each level of
gene lists. Each directory contains text files for each cluster. The ZIP file
is then encoded to base64 for easy download.
}
\details{
The function creates a temporary directory to store the files. For each level
in the `enrichr_format$gene_lists`, it creates a directory named after the
level. Within each level directory, it creates a text file for each cluster,
containing the genes in that cluster. The directories and files are added
to a ZIP file, which is then encoded to base64.
}
