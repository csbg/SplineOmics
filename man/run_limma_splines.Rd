% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_limma_splines.R
\name{run_limma_splines}
\alias{run_limma_splines}
\title{run_limma_splines()}
\usage{
run_limma_splines(splineomics)
}
\arguments{
\item{splineomics}{An S3 object of class `SplineOmics` that contains the
following elements:
\itemize{
  \item \code{data}: The matrix of the omics dataset, with the feature
  names optionally as row headers.
  \item \code{rna_seq_data}: An object containing the preprocessed
  RNA-seq data,
  such as the output from `limma::voom` or a similar preprocessing pipeline.
  \item \code{meta}: A dataframe containing metadata corresponding to the
  \code{data}, must include a 'Time' column and the column specified by
  \code{condition}.
  \item \code{design}: A character string representing the limma design
  formula.
  \item\code{dream_params}: #' A named list or NULL. When not NULL, it can 
                               contain the following named 
     elements:
     - `dof`: An integer greater than 1, specifying the degrees of freedom 
     for  the dream topTable. If set to 0, then the best dof is automatically
     found with the help of leave-one-out-crossvalidation (loocv). The dof 
     with the lowest error on the loocv is chosen.
     - `KenwardRoger`: A boolean indicating whether to use the Kenward-Roger 
     approximation for mixed models.
     Note that random effects are now directly specified in the design 
     formula and not in `dream_params`.
  \item \code{condition}: A character string specifying the column name
  in \code{meta} used to define groups for analysis.
  \item \code{spline_params}: A list of spline parameters used in the
  analysis, including:
    \itemize{
      \item \code{spline_type}: The type of spline (e.g., "n" for natural
      splines or "b" for B-splines).
      \item \code{dof}: Degrees of freedom for the spline.
      \item \code{knots}: Positions of the internal knots (for B-splines).
      \item \code{bknots}: Boundary knots (for B-splines).
      \item \code{degree}: Degree of the spline (for B-splines only).
    }
}}
}
\value{
The SplineOmics object, updated with a list with three elements:
        - `time_effect`: A list of top tables for each level with the time
                         effect.
        - `avrg_diff_conditions`: A list of top tables for each comparison
                                 between the levels. The comparison is the
                                 average difference of the values.
        - `interaction_condition_time`: A list of top tables for each
                                        comparison between levels. The
                                        comparison is the interaction between
                                        the condition and the time.
}
\description{
This function performs a limma spline analysis to identify significant
time-dependent changes in features (e.g., proteins) within an omics
time-series dataset. It evaluates features within each condition level
and between levels by comparing average differences and interactions
between time and condition.
}
