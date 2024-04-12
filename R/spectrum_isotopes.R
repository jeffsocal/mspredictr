#' Define the isotopes within a spectrum
#'
#' @description
#' `spectrum_isotopes()` assigns the predicted fragment masses of a given peptide
#' sequence with the mass spectrum.
#'
#' @param spectrum
#' A spectrum data object.
#'
#' @export
#'
#' @examples
#'  # using the supplied spectrum from the msreadr package
#'  library(msreadr)
#'  mzml <- path_to_example() |>
#'          read_spectra()
#'  mzml |>
#'    subset(spectrum_num == 1) |>
#'    spectrum_extract()
#'
spectrum_isotopes <- function(
    spectrum = NULL){

  # visible bindings
  n <- NULL
  mz <- NULL
  intensity <- NULL
  isotope_id <- NULL
  isotope_num <- NULL

  if(is.null(spectrum)) { cli::cli_abort("`spectrum` must not be null")}

  # convert to a dataframe for tidyverse
  if(is.matrix(spectrum)) {
    spectrum <- spectrum |> as.data.frame()
  }

  spectrum <- spectrum |>
    dplyr::mutate(isotope_id = which_isotopes(mz, intensity),
                  isotope_id = isotope_id |> stringr::str_pad(width = ceiling(log10(nrow(spectrum) + 1)), pad = "0"),
                  isotope_id = paste0("i", isotope_id)) |>
    dplyr::group_by(isotope_id) |>
    dplyr::mutate(isotope_num = dplyr::row_number() - 1,
                  n = dplyr::n(),
                  isotope_z = round(1 / mean(diff(mz))) ) |>
    dplyr::ungroup() |>
    dplyr::mutate(isotope_id = ifelse(n == 1, NA, isotope_id),
                  isotope_num = ifelse(n == 1, NA, isotope_num)) |>
    dplyr::select(-n)

  return(spectrum)
}
