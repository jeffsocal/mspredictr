#' A function to extract the plain amino acid sequence
#'
#' @description
#' `str_sequence()` Get just the amino acid sequence of a peptide
#'
#' @param sequences
#' The character string representing a peptide, or poly amino acid. The canonical
#' 20 amino acids are encoded in and chemical modifications can be represented by
#' and floating point numerical value enclosed by square brackets. If a canonical
#' amino acid is also enclosed in the square brackets `[M15.99]` it is assumed that
#' the numerical value is in addition to the mass of the residue, and thus represents
#' a post-translational modification (PTM).
#'
#' @export
#'
#' @examples
#' str_sequence('SA[M15.99]PLER')
#'
str_sequence <- function(
    sequences = NULL
){
  sequences |>
    stringr::str_extract_all("[A-Z]") |>
    lapply(function(x){paste(x, collapse = "")}) |>
    unlist()
}
