% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster_hits.R
\name{plot_spline_comparisons}
\alias{plot_spline_comparisons}
\title{Create spline comparison plots for two conditions}
\usage{
plot_spline_comparisons(
  time_effect_1,
  condition_1,
  time_effect_2,
  condition_2,
  avrg_diff_conditions,
  interaction_condition_time,
  data,
  meta,
  condition,
  X_1,
  X_2,
  plot_info,
  adj_pthresh_avrg_diff_conditions,
  adj_pthresh_interaction
)
}
\arguments{
\item{time_effect_1}{A data frame containing the time effects for the first
condition.}

\item{condition_1}{The name of the first condition.}

\item{time_effect_2}{A data frame containing the time effects for the second
condition.}

\item{condition_2}{The name of the second condition.}

\item{avrg_diff_conditions}{A data frame with the adjusted p-values for the
 average difference
between conditions.}

\item{interaction_condition_time}{A data frame with the adjusted p-values
 for the interaction between
condition and time.}

\item{data}{The data matrix containing the measurements.}

\item{meta}{The metadata associated with the measurements.}

\item{condition}{Column name of meta that contains the levels of the
experiment.}

\item{X_1}{A matrix of spline basis values for the first condition.}

\item{X_2}{A matrix of spline basis values for the second condition.}

\item{plot_info}{A list containing plotting information such as time unit
and axis labels.}

\item{adj_pthresh_avrg_diff_conditions}{The adjusted p-value threshold for
the average difference
between conditions.}

\item{adj_pthresh_interaction}{The adjusted p-value threshold for the
interaction between
condition and time.}
}
\value{
A list containing:
\describe{
  \item{plots}{A list of ggplot2 plots comparing the two conditions.}
  \item{feature_names}{A list of feature names for the plotted features.}
}
}
\description{
This function generates comparison plots for spline fits of two conditions
over time. It compares the time effects of two conditions, plots the data
points, and overlays the fitted spline curves. The function checks if the
adjusted p-values for the average difference between conditions and the
interaction between condition and time are below the specified thresholds
before generating plots.
}
