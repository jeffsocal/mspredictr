#' A helper function to extract just the spectrum of the spectra data object
#'
#' @param spectra
#' A spectra data object
#'
#' @param precursor
#' A boolean to indicate spectrum filtering should happen
#'
#' @param isotopes
#' A boolean to remove isotopes.
#'
#' @export
#'
spectrum_extract <- function(
    spectra = NULL,
    precursor = FALSE,
    isotopes = FALSE
){

  spectrum <- spectra$ms2$peaks

  if(precursor == TRUE){
    precursor_mz <- spectra$ms2$precursor_mz
  } else {
    precursor_mz <- 0
  }

  # precursor_z <- spectra$precursor_z
  # precursor_nm <- mass_neutral(precursor_mz, precursor_z)
  # n_expect <- round(precursor_nm / 114.35) * 3

  spectrum <- spectrum |> as.data.frame() |> spectrum_isotopes()
  colnames(spectrum)[1:2] <- c('mz', 'intensity')

  spectrum <- spectrum |>
    spectrum_denoise(
      precursor = precursor_mz,
      isotopes = isotopes)

  return(spectrum |> tibble::as_tibble())
}
