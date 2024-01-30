#' Generate a named mass ladder vector.
#'
#' @description
#' `mass_ladder()` Generates the lass ladder from a peptide sequence.
#'
#' @param sequence
#' The character string representing a peptide, or poly amino acid. The canonical
#' 20 amino acids are encoded in and chemical modifications can be represented by
#' and floating point numerical value enclosed by square brackets. If a canonical
#' amino acid is also enclosed in the square brackets `[M15.99]` it is assumed that
#' the numerical value is in addition to the mass of the residue, and thus represents
#' a post-translational modification (PTM).
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
