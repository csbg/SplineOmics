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
  list(spline_type = generate_random_string(), dof = as.integer(runif(1, 1, 10)))
}

# Function to generate varied inputs for run_limma_splines parameters
generate_run_limma_splines_inputs <- function() {
  list(
    data = if (runif(1) > 0.2) generate_random_dataframe() else generate_random_string(),
    meta = if (runif(1) > 0.2) generate_random_dataframe() else generate_random_string(),
    design = generate_random_string(),
    condition = generate_random_string(),
    spline_params = generate_random_list(),
    padjust_method = generate_random_string()
  )
}

# Wrapper to check if the function stops gracefully
test_function <- function(inputs) {
  tryCatch({
    browser()
    run_limma_splines(
      data = inputs$data,
      meta = inputs$meta,
      design = inputs$design,
      condition = inputs$condition,
      spline_params = inputs$spline_params,
      padjust_method = inputs$padjust_method
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
  inputs <- generate_run_limma_splines_inputs()
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
