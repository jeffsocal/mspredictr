#' helper function to read in platform specific results
#'
#' @param x location of file to parse
#'
#' @return a tibble
#'
read_comet <- function(
    x,
    cpus = 1
){

  proton_mass <- mass_proton()

  # scan is experiment scan level (eg ms1, ms2 included)
  out <- x |> readr::read_tsv(
    skip=1,
    show_col_types = FALSE
  ) |>
    dplyr::mutate(
      psm_score = -log10(`e-value`),
      psm_dp = ions_matched / ions_total
    ) |>
    dplyr::rename(
      ms_event = scan,
      # 1Th correction to get [M+H]+
      psm_mh = calc_neutral_mass + proton_mass,
      psm_rank = num,
      psm_peptide = modified_peptide,
      psm_sequence = plain_peptide,
      psm_protein = protein
    )

  out$psm_peptide <- out$psm_peptide |> lapply(str_peptide) |> unlist()

  return(out)
}
