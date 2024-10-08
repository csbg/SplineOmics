% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download_enrichr_databases.R
\name{enrichr_get_genesets}
\alias{enrichr_get_genesets}
\title{Get Enrichr Gene Sets}
\usage{
enrichr_get_genesets(databases)
}
\arguments{
\item{databases}{A character vector of database names to download from
Enrichr.}
}
\value{
A named list of gene sets from the specified Enrichr databases. Each
        database is represented as a list, with gene set names as list
        names and vectors of human gene symbols as list elements.
}
\description{
This function downloads gene sets from specified Enrichr databases.
It returns a list where each element is a list corresponding to a database,
with each element containing a vector of human gene symbols for a gene set.
}
