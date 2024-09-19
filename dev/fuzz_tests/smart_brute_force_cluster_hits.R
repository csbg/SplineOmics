# Load necessary libraries
library(testthat)
library(here)
library(ggplot2)

# Function to generate random dataframes
generate_random_dataframe <- function() {
  rows <- sample(5:20, 1)
  cols <- sample(5:10, 1)
  data.frame(matrix(runif(rows * cols), nrow = rows, ncol = cols))
}

# Function to generate random strings
generate_random_string <- function() {
  paste0(sample(letters, sample(1:10, 1), replace = TRUE), collapse = "")
}

# Function to generate random lists
generate_random_list <- function() {
  list(a = generate_random_string(), b = runif(1))
}

# Function to generate varied inputs for cluster_hits parameters
generate_cluster_hits_inputs <- function() {
  list(
    top_tables = list(data.frame(matrix(runif(sample(25:100, 1)), 
                                        ncol = sample(5:10, 1)))),
    data = if (runif(1) > 0.2) generate_random_dataframe() else generate_random_string(),
    meta = if (runif(1) > 0.2) generate_random_dataframe() else generate_random_string(),
    design = generate_random_string(),
    condition = generate_random_string(),
    report_info = generate_random_string(),
    spline_params = list(spline_type = c(generate_random_string()), 
                         dof = as.integer(runif(1, 1, 10))),
    adj_pthresholds = runif(sample(1:5, 1), 0, 0.1),
    clusters = ifelse(runif(1) > 0.5, "auto", as.integer(runif(1, 1, 10))),
    meta_batch_column = if (runif(1) > 0.5) generate_random_string() else NA,
    meta_batch2_column = if (runif(1) > 0.5) generate_random_string() else NA,
    time_unit = generate_random_string(),
    report_dir = example_report_dir,
    report = sample(c(TRUE, FALSE), 1)
  )
}

# Wrapper to check if the function stops gracefully
test_function <- function(inputs) {
  tryCatch({
    cluster_hits(
      top_tables = inputs$top_tables,
      data = inputs$data,
      meta = inputs$meta,
      design = inputs$design,
      condition = inputs$condition,
      report_info = inputs$report_info,
      spline_params = inputs$spline_params,
      adj_pthresholds = inputs$adj_pthresholds,
      clusters = inputs$clusters,
      meta_batch_column = inputs$meta_batch_column,
      meta_batch2_column = inputs$meta_batch2_column,
      time_unit = inputs$time_unit,
      report_dir = inputs$report_dir,
      report = inputs$report
    )
    return(list(success = FALSE, error = NULL, inputs = inputs))
  }, error = function(e) {
    if (inherits(e, "simpleError")) {
      return(list(success = TRUE, error = NULL, inputs = inputs))  # Expected stop error
    } else {
      return(list(success = FALSE, error = e, inputs = inputs))  # Unexpected error
    }
  })
}

# Run fuzz testing
set.seed(123)  # For reproducibility
results <- replicate(1000, {
  inputs <- generate_cluster_hits_inputs()
  test_function(inputs)
}, simplify = FALSE)

# Analyze the results
passed_tests <- sapply(results, function(x) x$success)
failed_tests <- results[!passed_tests]

cat("Number of tests passed: ", sum(passed_tests), "\n")
cat("Number of tests failed: ", length(failed_tests), "\n")

# Print details of failed tests
if (length(failed_tests) > 0) {
  for (i in seq_along(failed_tests)) {
    cat("Failed Test ", i, ":\n")
    print(failed_tests[[i]]$inputs)
    cat("Error: ", failed_tests[[i]]$error$message, "\n\n")
  }
}
