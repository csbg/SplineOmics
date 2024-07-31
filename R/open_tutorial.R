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
#' @importFrom rstudioapi navigateToFile
#' 
#' @export
#' 
open_tutorial <- function() {
  
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
