% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/explore_data.R
\name{explore_data}
\alias{explore_data}
\title{explore_data()}
\usage{
explore_data(splineomics, report_dir = here::here(), report = TRUE)
}
\arguments{
\item{splineomics}{A SplineOmics object, containing the data, meta,
condition, report_info, meta_batch_column, and
meta_batch2_column.}

\item{report_dir}{A non-empty string specifying the report directory. The
output HTML reports will be placed there. Default is the
current working dir, determined with the here library:
here::here().}

\item{report}{A Boolean TRUE or FALSE value, specifying if a report should
be generated or not. A report is generated per default, but
when only the plots as plot objects inside R are desired, this
argument can be set to FALSE.}
}
\value{
A list of ggplot objects representing various exploratory plots.
}
\description{
This function automatically generates exploratory data analysis (EDA) plots 
from the provided data. These include density plots, boxplots, PCA plots, 
MDS plots, variance explained plots, violin plots, mean correlation with 
time, first lag autocorrelation, lag1 differences, and coefficient of 
variation. The function returns all EDA plots as a list and, by default, 
creates an HTML report containing the plots, saving it to the specified 
report directory.
}
