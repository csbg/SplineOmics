#' utils scripts contains shared functions that are used by at least two package 
#' functions of the SplineOmics package.

# Level 1 internal functions ---------------------------------------------------


#' Create Progress Bar
#'
#' @description
#' Creates a progress bar for tracking the progress of an iterable task.
#'
#' @param iterable An iterable object (e.g., list or vector) whose length 
#' determines the total number of steps.
#'
#' @return A progress bar object from the 'progress' package.
#'
#' @examples
#' items <- 1:10
#' pb <- create_progress_bar(items)
#' for (i in items) {
#'   pb$tick()
#'   Sys.sleep(0.1)
#' }
#'
#' @seealso
#' \code{\link{progress_bar}}
#' 
create_progress_bar <- function(iterable) {
  
  library(progress)
  
  # Create and return the progress bar
  pb <- progress_bar$new(
    format = "  Processing [:bar] :percent :elapsed",
    total = length(iterable),
    width = 60
  )
  
  return(pb)
}