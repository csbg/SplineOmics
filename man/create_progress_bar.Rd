% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_general.R
\name{create_progress_bar}
\alias{create_progress_bar}
\title{utils scripts contains shared functions that are used by at least two package
functions of the SplineOmics package.
Create Progress Bar}
\usage{
create_progress_bar(iterable, message = "Processing")
}
\arguments{
\item{iterable}{An iterable object (e.g., list or vector) whose length
determines the total number of steps.}

\item{message}{A message to display with the progress bar
(default is "Processing").}
}
\value{
A progress bar object from the 'progress' package.
}
\description{
Creates a progress bar for tracking the progress of an iterable task.
}
\seealso{
\code{\link{progress_bar}}
}
