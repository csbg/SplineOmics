% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_input_validation.R
\name{Level2Functions}
\alias{Level2Functions}
\title{Level2Functions: A class providing level 2 functionalities}
\description{
Level2Functions: A class providing level 2 functionalities

Level2Functions: A class providing level 2 functionalities
}
\details{
This class provides various level 2 functionalities, including
methods to check dataframes and spline parameters.
}
\seealso{
\code{\link{InputControl}}

Check Data Matrix
}
\keyword{Check}
\keyword{Column}
\keyword{Columns}
\keyword{Dataframe}
\keyword{Degrees}
\keyword{Freedom}
\keyword{Maximum}
\keyword{Minimum}
\keyword{Parameters}
\keyword{Spline}
\keyword{This}
\keyword{Type}
\keyword{a}
\keyword{and}
\keyword{checks}
\keyword{column}
\keyword{columns}
\keyword{data}
\keyword{dataframe}
\keyword{expected}
\keyword{function}
\keyword{if}
\keyword{internal}
\keyword{match}
\keyword{names}
\keyword{of}
\keyword{respective}
\keyword{the}
\keyword{their}
\keyword{types.}
\section{Super classes}{
\code{\link[SplineOmics:Level4Functions]{SplineOmics::Level4Functions}} -> \code{\link[SplineOmics:Level3Functions]{SplineOmics::Level3Functions}} -> \code{Level2Functions}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-Level2Functions-check_data}{\code{Level2Functions$check_data()}}
\item \href{#method-Level2Functions-check_meta}{\code{Level2Functions$check_meta()}}
\item \href{#method-Level2Functions-check_dataframe}{\code{Level2Functions$check_dataframe()}}
\item \href{#method-Level2Functions-check_spline_params_generally}{\code{Level2Functions$check_spline_params_generally()}}
\item \href{#method-Level2Functions-check_spline_params_mode_dependent}{\code{Level2Functions$check_spline_params_mode_dependent()}}
\item \href{#method-Level2Functions-check_columns_spline_test_configs}{\code{Level2Functions$check_columns_spline_test_configs()}}
\item \href{#method-Level2Functions-check_spline_type_column}{\code{Level2Functions$check_spline_type_column()}}
\item \href{#method-Level2Functions-check_spline_type_params}{\code{Level2Functions$check_spline_type_params()}}
\item \href{#method-Level2Functions-check_max_and_min_dof}{\code{Level2Functions$check_max_and_min_dof()}}
\item \href{#method-Level2Functions-check_columns}{\code{Level2Functions$check_columns()}}
\item \href{#method-Level2Functions-clone}{\code{Level2Functions$clone()}}
}
}
\if{html}{\out{
<details open><summary>Inherited methods</summary>
<ul>
<li><span class="pkg-link" data-pkg="SplineOmics" data-topic="Level4Functions" data-id="create_error_message"><a href='../../SplineOmics/html/Level4Functions.html#method-Level4Functions-create_error_message'><code>SplineOmics::Level4Functions$create_error_message()</code></a></span></li>
<li><span class="pkg-link" data-pkg="SplineOmics" data-topic="Level3Functions" data-id="check_batch_column"><a href='../../SplineOmics/html/Level3Functions.html#method-Level3Functions-check_batch_column'><code>SplineOmics::Level3Functions$check_batch_column()</code></a></span></li>
<li><span class="pkg-link" data-pkg="SplineOmics" data-topic="Level3Functions" data-id="check_condition_time_consistency"><a href='../../SplineOmics/html/Level3Functions.html#method-Level3Functions-check_condition_time_consistency'><code>SplineOmics::Level3Functions$check_condition_time_consistency()</code></a></span></li>
<li><span class="pkg-link" data-pkg="SplineOmics" data-topic="Level3Functions" data-id="check_voom_structure"><a href='../../SplineOmics/html/Level3Functions.html#method-Level3Functions-check_voom_structure'><code>SplineOmics::Level3Functions$check_voom_structure()</code></a></span></li>
</ul>
</details>
}}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Level2Functions-check_data"></a>}}
\if{latex}{\out{\hypertarget{method-Level2Functions-check_data}{}}}
\subsection{Method \code{check_data()}}{
This function checks the validity of the data matrix, ensuring that it
is a
matrix, contains only numeric values,
has no missing values, and all elements are non-negative. Additionally, it
verifies that no rows or columns are
entirely zeros.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Level2Functions$check_data(data, data_meta_index = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{data}}{A dataframe containing numeric values.}

\item{\code{data_meta_index}}{An optional parameter specifying the index of the
 data
for error messages. Default is NA.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
Returns TRUE if all checks pass. Stops execution and returns an
error
message if any check fails.

Check Metadata
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Level2Functions-check_meta"></a>}}
\if{latex}{\out{\hypertarget{method-Level2Functions-check_meta}{}}}
\subsection{Method \code{check_meta()}}{
This function checks the validity of the metadata dataframe, ensuring it
contains the 'Time' column,
does not contain missing values, and that the specified condition column
 is
valid and of the appropriate type.
Additionally, it checks for an optional batch effect column and prints
messages regarding its use.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Level2Functions$check_meta(
  meta,
  condition,
  meta_batch_column = NULL,
  meta_batch2_column = NULL,
  data_meta_index = NULL
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{meta}}{A dataframe containing the metadata, including the 'Time'
column.}

\item{\code{condition}}{A single character string specifying the column name
in the
meta dataframe to be checked.}

\item{\code{meta_batch_column}}{An optional parameter specifying the column
name in
 the meta dataframe used to remove the batch effect. Default is NA.}

\item{\code{meta_batch2_column}}{An optional parameter specifying the column
name in
the meta dataframe used to remove the batch effect. Default is NA.}

\item{\code{data_meta_index}}{An optional parameter specifying the index of the
data/meta pair for error messages. Default is NA.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
Returns TRUE if all checks pass. Stops execution and returns an
error message if any check fails.

Check Dataframe
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Level2Functions-check_dataframe"></a>}}
\if{latex}{\out{\hypertarget{method-Level2Functions-check_dataframe}{}}}
\subsection{Method \code{check_dataframe()}}{
Validates that the dataframe contains all required columns with the
correct data types.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Level2Functions$check_dataframe(df)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{df}}{A dataframe to check.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
TRUE if the dataframe is valid, otherwise an error is thrown.

Check Spline Parameters Generally
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Level2Functions-check_spline_params_generally"></a>}}
\if{latex}{\out{\hypertarget{method-Level2Functions-check_spline_params_generally}{}}}
\subsection{Method \code{check_spline_params_generally()}}{
Validates the general structure and contents of spline parameters.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Level2Functions$check_spline_params_generally(spline_params)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{spline_params}}{A list of spline parameters.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
No return value, called for side effects.

Check Spline Parameters Mode Dependent
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Level2Functions-check_spline_params_mode_dependent"></a>}}
\if{latex}{\out{\hypertarget{method-Level2Functions-check_spline_params_mode_dependent}{}}}
\subsection{Method \code{check_spline_params_mode_dependent()}}{
Validates the spline parameters depending on the specified mode.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Level2Functions$check_spline_params_mode_dependent(
  spline_params,
  mode,
  meta,
  condition
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{spline_params}}{A list of spline parameters.}

\item{\code{mode}}{A character string specifying the mode
('integrated' or 'isolated').}

\item{\code{meta}}{A dataframe containing metadata.}

\item{\code{condition}}{A character string specifying the condition.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
No return value, called for side effects.

Check Columns in Spline Test Configurations
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Level2Functions-check_columns_spline_test_configs"></a>}}
\if{latex}{\out{\hypertarget{method-Level2Functions-check_columns_spline_test_configs}{}}}
\subsection{Method \code{check_columns_spline_test_configs()}}{
Validates that the spline test configurations contain the required columns
in the correct order.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Level2Functions$check_columns_spline_test_configs(spline_test_configs)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{spline_test_configs}}{A dataframe containing spline test
configurations.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
No return value, called for side effects.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Level2Functions-check_spline_type_column"></a>}}
\if{latex}{\out{\hypertarget{method-Level2Functions-check_spline_type_column}{}}}
\subsection{Method \code{check_spline_type_column()}}{
Validates that the 'spline_type' column in the spline test configurations
contains only 'n' or 'b'.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Level2Functions$check_spline_type_column(spline_test_configs)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{spline_test_configs}}{A dataframe containing spline test
configurations.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
No return value, called for side effects.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Level2Functions-check_spline_type_params"></a>}}
\if{latex}{\out{\hypertarget{method-Level2Functions-check_spline_type_params}{}}}
\subsection{Method \code{check_spline_type_params()}}{
Validates the parameters for each row in the spline test configurations
based on the spline type.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Level2Functions$check_spline_type_params(spline_test_configs)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{spline_test_configs}}{A dataframe containing spline test
configurations.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
TRUE if all checks pass, otherwise an error is thrown.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Level2Functions-check_max_and_min_dof"></a>}}
\if{latex}{\out{\hypertarget{method-Level2Functions-check_max_and_min_dof}{}}}
\subsection{Method \code{check_max_and_min_dof()}}{
Validates the degrees of freedom (DoF) for each row in the spline test
configurations based on the metadata.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Level2Functions$check_max_and_min_dof(spline_test_configs, metas)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{spline_test_configs}}{A dataframe containing spline test
configurations.}

\item{\code{metas}}{A list of metadata corresponding to the data matrices.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
No return value, called for side effects.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Level2Functions-check_columns"></a>}}
\if{latex}{\out{\hypertarget{method-Level2Functions-check_columns}{}}}
\subsection{Method \code{check_columns()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Level2Functions$check_columns(df, expected_cols)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{df}}{A dataframe to check.}

\item{\code{expected_cols}}{A character vector of expected column names.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
This function does not return a value. It stops execution if the
dataframe columns or their classes do not match the expected structure.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Level2Functions-clone"></a>}}
\if{latex}{\out{\hypertarget{method-Level2Functions-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Level2Functions$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
