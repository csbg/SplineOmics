% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_gsea.R
\name{dbs_to_term2genes}
\alias{dbs_to_term2genes}
\title{Convert Database File to TERM2GENE List}
\usage{
dbs_to_term2genes(databases)
}
\arguments{
\item{databases}{A dataframe, containing the three columns DB, Geneset, and
gene. This dataframe contains the databases downloaded from
Enrichr with the SplineOmics package function:
download_enrichr_databases.}
}
\value{
A nested list where the first level of names corresponds to database
names ('DB'),
        the second level to gene sets ('Geneset'), and the innermost lists
        contain gene names ('Gene') associated with each gene set.
}
\description{
Reads a specified .tsv file containing information about databases,
gene sets, and genes. The file should have three columns: 'DB' for database
names, Geneset' for gene set identifiers, and 'Gene' for gene names. This
function organizes this information into a nested list. Each top-level
element corresponds to a unique database, and within each, gene sets map to
lists of associated genes.
}
