#' A non-iterable helper function to clean up comet sequence strings
#'
#' @description
#' `munge_seq()` Get a normalized peptide sequence string
#'
#' @param sequence
#' The character string representing a peptide, or poly amino acid. The canonical
#' 20 amino acids are encoded in and chemical modifications can be represented by
#' and floating point numerical value enclosed by square brackets. This function
#' takes a string input and capitalizes all letters and creates a square bracket
#' containg any mass values along with the immediate-left letter.
#'
#'
munge_seq <- function(
    sequence = NULL
){

  if(!is.character(sequence)) { cli::cli_abort("`sequence` must be a character string")}

  seq_regex <- "[a-zA-Z]|\\[.+?\\]";
  mod_regex <- "[a-zA-Z]|[\\-\\+]*[0-9]+\\.*[0-9]*";
  seq_array <- stringr::str_extract_all(sequence, seq_regex)[[1]]

  w_ptm <- which(grepl("\\[|[0-9]", seq_array))

  if(length(intersect(c(1,3), w_ptm)) == 2) {
    mass_1 <- gsub("\\[|\\]", "", seq_array[1]) |> as.numeric()
    mass_2 <- gsub("\\[|\\]", "", seq_array[3]) |> as.numeric()

    seq_array[3] <- paste0("[", seq_array[2], (mass_1 + mass_2), "]")
    seq_array[1] <- ''
    seq_array[2] <- ''
    w_ptm <- setdiff(w_ptm, c(1,3))
  }

  if(length(w_ptm) != 0) {
    for(i_ptm in w_ptm){
      seq_array[i_ptm] <- gsub("\\[|\\]", "", seq_array[i_ptm]) |> as.numeric() |> num_trunc()

      if(i_ptm == 1){
        seq_array[i_ptm + 1] <- paste0("[", seq_array[i_ptm + 1], seq_array[i_ptm], "]")
      } else {
        seq_array[i_ptm - 1] <- paste0("[", seq_array[i_ptm - 1], seq_array[i_ptm], "]")
      }
      seq_array[i_ptm] <- ''
    }
  }

  return(paste(seq_array, collapse = '') |> stringr::str_to_upper())
}

#' An iterable function to clean up peptide sequence strings.
#'
#' @description
#' `str_peptide()` Get a normalized peptide sequence string
#'
#' @param sequences
#' The character string representing a peptide, or poly amino acid. The canonical
#' 20 amino acids are encoded in and chemical modifications can be represented by
#' and floating point numerical value enclosed by square brackets. This function
#' takes a string input and capitalizes all letters and creates a square bracket
#' containg any mass values along with the immediate-left letter.
#'
#' @export
#'
#' @examples
#' str_peptide('sam[15.99]ple')
#'
#' str_peptide('[47.89]sample')
#'
#' str_peptide('n[47.89]SAMPLE')
#'
str_peptide <- function(
    sequences = NULL
){
  sequences |> str_clean() |> lapply(munge_seq) |> unlist()
}

#' An iterable helper function to clean up comet sequence strings
#'
#' @description
#' `str_peptide()` Get a normalized peptide sequence string
#'
#' @param sequences
#' The character string representing a peptide, or poly amino acid. The canonical
#' 20 amino acids are encoded in and chemical modifications can be represented by
#' and floating point numerical value enclosed by square brackets. This function
#' takes a string input and capitalizes all letters and creates a square bracket
#' containg any mass values along with the immediate-left letter.
#'
str_clean <- function(
    sequences = NULL
){
  sequences |> lapply(
    function(sequence){
      sequence <- sub("^.\\.", "", sequence)
      sequence <- sub("^n", "", sequence)
      sequence <- sub("\\..$", "", sequence)
      sequence <- sequence |> stringr::str_to_upper()
      return(sequence)
    }
  ) |> unlist()
}
