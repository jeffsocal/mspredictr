#' Convert a peptide string to a named variable
#'
#' @description
#' `mass_ladder()` Generates a named variable
#'
#' @param sequence as character string
#'
#' @return a vector
#'
mass_ladder <- function(
    sequence = NULL
    ){

  if(!is.character(sequence)) { cli::cli_abort("`sequence` must be a character string")}

  seq_regex <- "[a-zA-Z]|\\[.+?\\]";
  mod_regex <- "[a-zA-Z]|[\\-\\+]*[0-9]+\\.*[0-9]*";
  seq_array <- stringr::str_extract_all(sequence, seq_regex)[[1]]
  mass_array <- rep(0, length(seq_array))

  for(i in 1:length(seq_array)){
    residue <- stringr::str_extract_all(seq_array[i], mod_regex)[[1]]
    if(grepl("[0-9]+", residue[1])){
      mass_array[i] <- as.numeric(residue[1])
    } else {
      mass_array[i] <- mass_residue(residue[1])
    }
    if(length(residue) == 2) {
      mass_array[i] <- mass_array[i] + as.numeric(residue[2])
      }
  }
  names(mass_array) <- seq_array

  return(mass_array)
}
