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

# Function to generate varied inputs for limma_report parameters
generate_limma_report_inputs <- function() {
  list(
    run_limma_splines_result = list(
      time_effect = list(top_table1 = if (runif(1) > 0.2) generate_random_dataframe() else generate_random_string()),
      avrg_diff_conditions = list(top_table2 = if (runif(1) > 0.2) generate_random_dataframe() else generate_random_string()),
      interaction_condition_time = list(top_table3 = if (runif(1) > 0.2) generate_random_dataframe() else generate_random_string())
    ),
    report_info = generate_random_list(),
    adj_pthresh = runif(1, 0, 0.1),
    report_dir = example_report_dir
  )
}

# Wrapper to check if the function stops gracefully
test_function <- function(inputs) {
  tryCatch({
    limma_report(
      run_limma_splines_result = inputs$run_limma_splines_result,
      report_info = inputs$report_info,
      adj_pthresh = inputs$adj_pthresh,
      report_dir = inputs$report_dir
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
  inputs <- generate_limma_report_inputs()
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
