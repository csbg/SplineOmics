% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/screen_limma_hyperparams.R
\name{remove_batch_effect}
\alias{remove_batch_effect}
\title{Remove Batch Effect}
\usage{
remove_batch_effect(
  datas,
  metas,
  meta_batch_column,
  meta_batch2_column,
  condition
)
}
\arguments{
\item{datas}{A list of matrices.}

\item{metas}{A list of metadata corresponding to the data matrices.}

\item{meta_batch_column}{A character string specifying the meta batch column.}

\item{meta_batch2_column}{A character string specifying the second meta batch
column.}

\item{condition}{A character vector of length 1, specifying the column name
of the meta dataframe, that contains the levels that
separate the experiment.}
}
\value{
A list of matrices with batch effects removed where applicable.
}
\description{
Removes batch effects from the data matrices using the specified batch
column in the metadata.
}
\seealso{
\link[limma]{removeBatchEffect}
}
