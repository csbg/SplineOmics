# This function is exported, but not part of the functionality of the 
# SplineOmics package. Rather, it provides a convenient way of opening the 
# tutorial for the SplineOmics package in a R Markdown file, which
# provides an interactive experience.


#' Interactive Tutorial for Getting Started
#'
#' @description
#' This function opens the `tutorial.Rmd` file in RStudio for
#' interactive use. Users can then run each code chunk step by step.
#' 
#' @export
#' 
open_tutorial <- function() {
  
  # Check if rstudioapi is installed
  if (!requireNamespace(
    "rstudioapi",
    quietly = TRUE
  )) {
    repeat {
      # Prompt the user for action
      cat("The 'rstudioapi' package is not installed.\n")
      cat("1: Install 'rstudioapi'\n")
      cat("2: Do not install and quit\n")
      cat("3: Resolve manually and retry\n")
      choice <- readline(prompt = "Please enter your choice (1, 2, or 3): ")
      
      # Check user input and take appropriate action
      if (choice == "1") {
        install.packages("rstudioapi")
        break
      } else if (choice == "2") {
        stop(
          "User chose not to install 'rstudioapi'. Exiting function.",
          call. = FALSE
          )
      } else if (choice == "3") {
        stop(
          "Please install 'rstudioapi' manually and retry.",
          call. = FALSE)
      } else {
        cat("Invalid choice. Please enter 1, 2, or 3.\n")
      }
    }
  }
  
  file <- system.file(
    "tutorial",
    "tutorial.Rmd",
    package = "SplineOmics"
    )
  if (file != "") {
    if (rstudioapi::isAvailable()) {
      rstudioapi::navigateToFile(file)
    } else {
      stop("RStudio API not available. Cannot open tutorial.")
    }
  } else {
    stop("tutorial.Rmd file not found under inst/tutorial")
  }
}


#' Open Template for Quick Setup
#'
#' @description
#' This function opens the `template.Rmd` file in RStudio for
#' interactive use. The template file provides a structure for users 
#' to quickly set up their personal analysis.
#' 
#' @export
#' 
open_template <- function() {
  
  # Check if rstudioapi is installed
  if (!requireNamespace(
    "rstudioapi",
    quietly = TRUE
    )) {
    repeat {
      # Prompt the user for action
      cat("The 'rstudioapi' package is not installed.\n")
      cat("1: Install 'rstudioapi'\n")
      cat("2: Do not install and quit\n")
      cat("3: Resolve manually and retry\n")
      choice <- readline(prompt = "Please enter your choice (1, 2, or 3): ")
      
      # Check user input and take appropriate action
      if (choice == "1") {
        install.packages("rstudioapi")
        break
      } else if (choice == "2") {
        stop(
          "User chose not to install 'rstudioapi'. Exiting function.",
          call. = FALSE
        )
      } else if (choice == "3") {
        stop(
          "Please install 'rstudioapi' manually and retry.",
          call. = FALSE)
      } else {
        cat("Invalid choice. Please enter 1, 2, or 3.\n")
      }
    }
  }
  
  file <- system.file(
    "template",
    "template.Rmd",
    package = "SplineOmics"
  )
  if (file != "") {
    if (rstudioapi::isAvailable()) {
      rstudioapi::navigateToFile(file)
    } else {
      stop("RStudio API not available. Cannot open template.")
    }
  } else {
    stop("template.Rmd file not found under inst/template")
  }
}


