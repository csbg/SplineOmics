% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/open_tutorial.R
\name{open_tutorial}
\alias{open_tutorial}
\title{open_tutorial()}
\usage{
open_tutorial()
}
\value{
If successful, opens the `tutorial.Rmd` file in RStudio for the user to
interact with.
If `rstudioapi` is not installed or available, or the tutorial file is
not found,
an error is thrown with a corresponding message.
}
\description{
This function opens the `tutorial.Rmd` file in RStudio for
interactive use. Users can then run each code chunk step by step.
}
