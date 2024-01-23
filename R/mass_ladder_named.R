#' Convert a peptide string to a named variable
#'
#' @description
#' `mass_ladder()` Generates a named variable
#'
#' @param sequence as character string
#'
#' @return a vector
#'
mass_ladder_named <- function(
    sequence = NULL
    ){

  if(!is.character(sequence)) { cli::cli_abort("`sequence` must be a character string")}

  seq_regex <- "[a-zA-Z]|\\[.+?\\]";
  mod_regex <- "[a-zA-Z]|[\\-\\+]*[0-9]+\\.*[0-9]*";
  seq_array <- stringr::str_extract_all(sequence, seq_regex)[[1]]
  mass_array <- mass_ladder(sequence)
  names(mass_array) <- seq_array

  return(mass_array)
}
