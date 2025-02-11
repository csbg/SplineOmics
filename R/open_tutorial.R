#' open_tutorial()
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
      "The 'rstudioapi' package is not installed.\n",
      "Please install it manually (to your custom_lib_path) using the command", 
      "below and re-run the function:\n\n",
      "  install.packages('rstudioapi')\n\n",
      "This is an optional dependency of the SplineOmics package, ",
      "only needed for optional functionality and not part of the core", 
      "package, which is why it must be installed manually if this function", 
      "is used."
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
