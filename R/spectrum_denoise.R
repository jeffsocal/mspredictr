#' A helper function used to clean up a mass spectrum to expose only the most abundant peaks
#'
#' @description
#' `denoise_spectrum()` Removes peaks from a given mass spectrum filtered to remove residual
#' precursor peaks, isotopes peaks and any low-level noise.
#'
#' @param spectrum
#' A spectrum data object.
#'
#' @param precursor
#' The floating point value of precursor mass to remove.
#'
#' @param hedge
#' The mz space to remove on either side of a tall peak when filtering.
#'
#' @param n
#' The number of peaks to retain, ranked by abundance.
#'
#' @export
#'
spectrum_denoise <- function(
    spectrum = NULL,
    precursor = NULL,
    isotopes = FALSE
){

  if(is.null(spectrum)) { cli::cli_abort("`spectrum` must not be null")}
  if(!is.null(precursor) & !is.numeric(precursor)) { cli::cli_abort("`precursor` must be a numeric")}
  if(!is.logical(isotopes)) { cli::cli_abort("`isotopes` must be a T/F")}

  if(!is.null(precursor)) {
    mz <- spectrum[,1]
    w_trim <- which(abs(mz - (precursor + 1.3)) > 1.5)
    spectrum <- spectrum[w_trim,]
  }

  if(isotopes == TRUE) {
    spectrum <- spectrum |>
      dplyr::filter(!is.na(isotope_num), isotope_num == 0)
  }

  return(spectrum)
}
