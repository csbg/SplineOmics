# This function is exported, but not part of the functionality of the 
# SplineOmics package. Rather, it provides a convenient way of opening the 
# demo (tutorial) for the SplineOmics package in a R Markdown file, which
# provides an interactive experience.


#' Interactive Demo for Getting Started
#'
#' @description
#' This function opens the `demo.Rmd` file in RStudio for
#' interactive use. Users can then run each code chunk step by step.
#'
#' @importFrom rstudioapi navigateToFile
#' 
#' @export
#' 
splineomics_demo <- function() {
  
  file <- system.file(
    "demo",
    "demo.Rmd",
    package = "SplineOmics"
    )
  if (file != "") {
    if (rstudioapi::isAvailable()) {
      rstudioapi::navigateToFile(file)
    } else {
      stop("RStudio API not available. Please open the file manually: ", file)
    }
  } else {
    stop("Demo file not found: demo.Rmd file not found under inst/demo.")
  }
}
