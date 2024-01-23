#' Helper function for displaying path to data
#'
#' @return print the table to console
#' @export
#'
path_to_example <- function(
){
  list.files(system.file("extdata", "", package = "mspredictr"), full.names = T, pattern = "\\.mzML")[1]
}
