library(testthat)

# Simple function definition
add_two <- function(x) {
  if (!is.numeric(x)) {
    stop("Input must be numeric")
  }
  return(x + 2)
}
# Unit test for the add_two function
library(testthat)
test_that("add_two works correctly", {
  # Test that add_two returns the expected output
  expect_equal(add_two(2), 4)
  # Test that add_two handles non-numeric input
  expect_error(add_two("a"), "Input must be numeric")
})