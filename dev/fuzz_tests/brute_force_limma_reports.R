# Load necessary libraries
library(devtools)
devtools::load_all()

library(testthat)
library(here)

# Function to generate random inputs of various types
generate_random_input <- function() {
  types <- list(
    logical = sample(c(TRUE, FALSE), 1),
    integer = sample(-100:100, 1),
    double = runif(1, -100, 100),
    character = paste0(sample(letters, sample(1:10, 1), replace = TRUE), collapse = ""),
    factor = factor(sample(letters, sample(1:10, 1), replace = TRUE)),
    date = Sys.Date() + sample(-100:100, 1),
    dataframe = data.frame(matrix(runif(100), nrow = 10)),
    list = list(sample(letters, 1))
  )
  return(sample(types, 1)[[1]])
}

# Function to generate random inputs for limma_report parameters
generate_limma_report_inputs <- function() {
  list(
    run_limma_splines_result = list(
      time_effect = generate_random_input(),
      avrg_diff_conditions = generate_random_input(),
      interaction_condition_time = generate_random_input()
    ),
    report_info = generate_random_input(),
    adj_pthresh = runif(1, 0, 1),
    report_dir = generate_random_input()
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
results <- replicate(10000, {
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
