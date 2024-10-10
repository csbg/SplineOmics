#' Interactive Tutorial for Getting Started
#'
#' @description
#' This function opens the `tutorial.Rmd` file in RStudio for
#' interactive use. Users can then run each code chunk step by step.
#'
#' @return
#' If successful, opens the `tutorial.Rmd` file in RStudio for the user to
#' interact with.
#' If `rstudioapi` is not installed or available, or the tutorial file is
#' not found,
#' an error is thrown with a corresponding message.
#'
#' @export
#'
open_tutorial <- function() {
  # Check if rstudioapi is installed
  if (!requireNamespace("rstudioapi", quietly = TRUE)) {
    stop_call_false(
      "The 'rstudioapi' package is not installed. ",
      "Please install it manually with: install.packages('rstudioapi')."
    )
  }

  # Find the tutorial file
  file <- system.file("tutorial", "tutorial.Rmd", package = "SplineOmics")

  if (file != "") {
    if (rstudioapi::isAvailable()) {
      rstudioapi::navigateToFile(file)
    } else {
      stop_call_false("RStudio API not available. Cannot open tutorial.")
    }
  } else {
    stop_call_false("tutorial.Rmd file not found under inst/tutorial")
  }
}


#' Open Template for Quick Setup
#'
#' @description
#' This function opens the `template.Rmd` file in RStudio for
#' interactive use. The template file provides a structure for users
#' to quickly set up their personal analysis.
#'
#' @return
#' If successful, opens the `template.Rmd` file in RStudio for the user to
#' interact with.
#' If `rstudioapi` is not installed or available, or the template file is
#' not found,
#' an error is thrown with a corresponding message.
#'
#' @export
#'
open_template <- function() {
  # Check if rstudioapi is installed
  if (!requireNamespace("rstudioapi", quietly = TRUE)) {
    stop(
      "The 'rstudioapi' package is not installed. ",
      "Please install it manually with: install.packages('rstudioapi').",
      call. = FALSE
    )
  }

  # Find the template file
  file <- system.file("template", "template.Rmd", package = "SplineOmics")

  if (file != "") {
    if (rstudioapi::isAvailable()) {
      rstudioapi::navigateToFile(file)
    } else {
      stop("RStudio API not available. Cannot open template.", call. = FALSE)
    }
  } else {
    stop("template.Rmd file not found under inst/template", call. = FALSE)
  }
}
