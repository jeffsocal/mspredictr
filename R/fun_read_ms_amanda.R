#' helper function to read in platform specific results
#'
#' @param x location of file to parse
#'
#' @return a tibble
#'
read_ms_amanda <- function(
    x,
    cpus = 1
){
  # Scan Number is the experiment scan level (eg ms1, ms2 included)
  out <- x |>
    readr::read_tsv(
      skip=1,
      show_col_types = FALSE
    ) |>
    dplyr::mutate(ms_event = ifelse(`Scan Number` == 0, Title, `Scan Number`),
           psm_dp = `Nr of matched peaks` / `number of considered fragment ions`) |>
    dplyr::rename(
      psm_score = `Amanda Score`,
      psm_rank = Rank,
      psm_sequence = Sequence,
      psm_protein = `Protein Accessions`
    )

  # normalize peptide
  out$psm_sequence <- out$psm_sequence |> lapply(str_clean) |> unlist()
  out <- out |>
    dplyr::mutate(psm_peptide = purrr::map2(psm_sequence, `Modifications`, ms_amanda_peptide) |> unlist()) |>
    dplyr::select(!c('Filename','Scan Number'))
  # compute psm mass
  out$psm_nmass <- out$psm_peptide |> lapply(peptide_mass) |> unlist()

  return(out)
}

#' helper function to read in platform specific results
#'
#' @param sequence string
#' @param modifications string
#'
#' @return string
#'
ms_amanda_peptide <- function(
    sequence = NULL,
    modifications = NULL
){

  sequence <- sequence |> stringr::str_extract_all('[A-Z]') |> unlist()
  modifications <- modifications |> stringr::str_split(";") |> unlist()

  masses <- modifications |> stringr::str_extract("[0-9]+\\.[0-9]+") |> as.numeric()
  local <- modifications |> stringr::str_extract("(?<=[A-Z])[0-9]+")

  out <- ''
  for(i in 1:length(sequence)){
    if(i %in% local){
      sequence[i] <- paste0("[", sequence[i], num_trunc(masses[which(i == local)],2), "]")
    }
    out <- paste0(out, sequence[i])
  }
  return(out)
}
