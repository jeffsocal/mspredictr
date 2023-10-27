#' non-iteratable helper function to clean up comet sequence strings
#'
#' @description
#' `munge_seq()` get the mass of a poly amino acid
#'
#' @param sequence as character string
#'
#' @return a string
#'
munge_seq <- function(
    sequence = NULL
){

  if(!is.character(sequence)) { cli::cli_abort("`sequence` must be a character string")}

  sequence <- sub("^.\\.", "", sequence)
  sequence <- sub("^n", "", sequence)
  sequence <- sub("\\..$", "", sequence)

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

#' iteratable function to clean up comet sequence strings
#'
#' @description
#' `str_peptide()` get the mass of a poly amino acid
#'
#' @param sequence as character string
#'
#' @return a string vector
#' @export
#'
str_peptide <- function(
    sequences = NULL
){

  sequences |> lapply(munge_seq) |> unlist()

}
