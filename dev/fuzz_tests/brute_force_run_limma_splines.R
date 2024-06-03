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

# Function to generate random values of various types
generate_random_value <- function() {
  types <- list(
    logical = sample(c(TRUE, FALSE), 1),
    integer = sample(-100:100, 1),
    double = runif(1, -100, 100),
    character = generate_random_string(),
    factor = factor(generate_random_string()),
    date = Sys.Date() + sample(-100:100, 1),
    dataframe = generate_random_dataframe(),
    list = generate_random_list()
  )
  return(sample(types, 1)[[1]])
}

# Function to generate varied inputs for run_limma_splines parameters
generate_run_limma_splines_inputs <- function() {
  list(
    data = generate_random_value(),
    meta = generate_random_value(),
    design = generate_random_value(),
    condition = generate_random_value(),
    spline_params = generate_random_value(),
    padjust_method = generate_random_value()
  )
}

# Wrapper to check if the function stops gracefully
test_function <- function(inputs) {
  tryCatch({
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
