#' A function to extract the plain amino acid sequence
#'
#' @description
#' `specdist()` Get just the amino acid sequence of a peptide
#'
#' @param seq1
#' The character string representing a peptide, or poly amino acid. The canonical
#' 20 amino acids are encoded in and chemical modifications can be represented by
#' and floating point numerical value enclosed by square brackets. If a canonical
#' amino acid is also enclosed in the square brackets `[M15.99]` it is assumed that
#' the numerical value is in addition to the mass of the residue, and thus represents
#' a post-translational modification (PTM).

#' @param seq2
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
#' specdist('SA[M15.99]PLER', 'SA[M15.99]PLQR')
#'
specdist <- function(
    seq1 = NULL,
    seq2 = NULL
){
  indf1 <- indef_fragments(seq1)
  indf2 <- indef_fragments(seq2)
  length(setdiff(indf1, indf2)) + length(setdiff(indf2, indf1))
}
