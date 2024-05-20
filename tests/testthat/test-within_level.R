library(testthat)

# tests/testthat/test-within_level.R
library(limma)
library(splines)

# Define the test function
test_that("within_level handles input parameters correctly", {
  # Mock data input
  set.seed(123)
  data <- matrix(rnorm(100), nrow = 10, ncol = 10)
  
  # Mock meta input
  meta <- data.frame(
    Time = rep(1:5, each = 2),
    Condition = rep(c("A", "B"), each = 5)
  )
  
  # Function parameters
  design <- "~ Condition + X"
  factor <- "Condition"
  level <- "A"
  level_index <- 1
  padjust_method <- "BH"
  spline_params <- list(
    spline_type = c("n"),
    dof = c(2L)
  )
  
  # Run the function
  result <- within_level(data, 
                         meta, 
                         design, 
                         factor, 
                         level, 
                         spline_params, 
                         level_index, 
                         padjust_method)
  
  # Check that the result is a list
  expect_type(result, "list")
  
  # Check that the list contains two elements
  expect_length(result, 2)
  
  # Check the names of the elements in the list
  expect_named(result, c("top_table", "fit"))
  
  # Check if top_table is a data frame
  expect_true(is.data.frame(result$top_table))
  
  # Check if fit is an MArrayLM object
  expect_true(inherits(result$fit, "MArrayLM"))
  
  # Verify that the p-value adjustment method is used in the top table
  top_table_adjust_method <- attr(result$top_table, "adjust.method")
  expect_equal(top_table_adjust_method, padjust_method)
  
  # Check the dimensions of top_table
  expect_equal(nrow(result$top_table), nrow(data))
})

test_that("within_level handles different spline_params correctly", {
  # Mock data input
  set.seed(123)
  data <- matrix(rnorm(100), nrow = 10, ncol = 10)
  
  # Mock meta input
  meta <- data.frame(
    Time = rep(1:5, each = 2),
    Condition = rep(c("A", "B"), each = 5)
  )
  
  # Function parameters
  design <- "~ Condition + X"
  factor <- "Condition"
  level <- "A"
  padjust_method <- "BH"
  level_index <- 2  # Test with the second set of spline params
  
  # Test with different spline_params
  spline_params <- list(
    spline_type = c("n"),
    dof = c(2L)
  )
  
  result <- within_level(data, 
                         meta, 
                         design, 
                         factor,
                         level, 
                         spline_params, 
                         level_index,
                         padjust_method)
  
  # Check that the correct degrees of freedom or knots are used
  if (!is.null(spline_params$dof)) {
    expect_equal(attr(meta$X, "df"), spline_params$dof[level_index])
  } else {
    expect_equal(attr(meta$X, "knots"), spline_params$knots[[level_index]])
  }
})

# test_that("within_level handles different padjust_methods correctly", {
#   # Mock data input
#   set.seed(123)
#   data <- matrix(rnorm(100), nrow = 10, ncol = 10)
#   
#   # Mock meta input
#   meta <- data.frame(
#     Time = rep(1:5, each = 2),
#     Condition = rep(c("A", "B"), each = 5)
#   )
#   
#   # Function parameters
#   design <- "~ Condition + X"
#   factor <- "Condition"
#   level <- "A"
#   level_index <- 1
#   spline_params <- list(
#     spline_type = c("n"),
#     dof = c(2L)
#   )
#   
#   supported_methods <- c("holm", "hochberg", "hommel", "bonferroni", 
#                          "BH", "BY", "fdr", "none")
#   
#   for (method in supported_methods) {
#     result <- within_level(data, meta, design, factor, level, spline_params, 
#                            level_index, method)
#     expect_equal(result$fit$adjust.method, method)
#   }
# })
