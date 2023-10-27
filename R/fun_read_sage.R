#' helper function to read in platform specific results
#'
#' @param x location of file to parse
#'
#' @return a tibble
#'
read_sage <- function(
    x,
    cpus = 1
){
  # scannr contains experiment scan level (eg ms1, ms2 included)
  out <- x |> readr::read_tsv(show_col_types = FALSE) |>
    dplyr::rename(
      ms_event = scannr,
      psm_rank = rank,
      psm_nmass = calcmass,
      psm_score = sage_discriminant_score,
      psm_peptide = peptide,
      psm_protein = proteins
    ) |>
    dplyr::mutate(ms_event = stringr::str_remove(ms_event, ".+scan\\=") |> as.numeric(),
                  psm_dp = matched_peaks / (peptide_len * 2 * charge)) |>
    dplyr::select(dplyr::matches('file|ms_event|origin|decoy|psm')) |>
    dplyr::select(!filename)

  out$psm_peptide <- out$psm_peptide |> lapply(str_peptide) |> unlist()
  out$psm_sequence <- out$psm_peptide |> lapply(str_sequence) |> unlist()

  return(out)
}
