% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/explore_data.R
\name{make_violin_box_plots}
\alias{make_violin_box_plots}
\title{Generate Violin Box Plot}
\usage{
make_violin_box_plots(data, meta, condition)
}
\arguments{
\item{data}{A numeric matrix containing the data.}

\item{meta}{A dataframe containing the column meta data of data}

\item{condition}{The name of the factor column of meta for the experiment}
}
\value{
A ggplot object representing the violin plot.
}
\description{
This function generates a violin plot for a given data matrix. The violin
plot shows the distribution of the values in the data matrix across different
variables, with each variable's distribution displayed as a separate violin.
}
