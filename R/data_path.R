#' Helper function for displaying path to data
#'
#' @export
#'
path_to_example <- function(
){
  list.files(system.file("extdata", "", package = "mspredictr"),
             full.names = T, pattern = "\\.csv")[1]
}
