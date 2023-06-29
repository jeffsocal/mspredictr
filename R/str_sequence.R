#' iteratable function to clean up comet sequence strings
#'
#' @description
#' `str_sequence()` get the mass of a poly amino acid
#'
#' @param sequence a vector of character string
#'
#' @return a string vector
#' @export
#'
str_sequence <- function(
    sequences = NULL
){
  sequences |>
    stringr::str_extract_all("[A-Z]") |>
    lapply(function(x){paste(x, collapse = "")}) |>
    unlist()
}
