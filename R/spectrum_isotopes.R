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
#'  tbl <- path_to_example() |>
#'          readr::read_csv()
#'  tbl |>
#'    spectrum_isotopes()
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
    dplyr::mutate(isotope_id = group_isotopes(mz * 1.0, intensity * 1.0),
                  isotope_num = label_isotopes(isotope_id),
                  isotope_id = isotope_id |> stringr::str_pad(width = ceiling(log10(nrow(spectrum) + 1)), pad = "0"),
                  isotope_id = paste0("i", isotope_id)) |>
    dplyr::group_by(isotope_id) |>
    dplyr::mutate(n = dplyr::n(),
                  isotope_z = round(1 / mean(diff(mz))) ) |>
    dplyr::ungroup() |>
    dplyr::mutate(isotope_id = ifelse(n == 1, NA, isotope_id),
                  isotope_num = ifelse(n == 1, NA, isotope_num)) |>
    dplyr::select(-n)

  return(spectrum)
}
