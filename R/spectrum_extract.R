#' A helper function to extract just the spectrum of the spectra data object
#'
#' @param spectra
#' A spectra data object
#'
#' @param filter
#' A boolean to indicate spectrum filtering should happen
#'
#' @export
#'
spectrum_extract <- function(
    spectra = NULL,
    filter = FALSE
){
  spectrum <- spectra$peaks
  precursor_mz <- spectra$precursor_mz
  precursor_z <- spectra$precursor_z
  precursor_nm <- mass_neutral(precursor_mz, precursor_z)

  n_expect <- round(precursor_nm / 114.35) * 3
  spectrum <- spectrum |> as.data.frame()
  colnames(spectrum) <- c('mz', 'intensity')

  if(filter == TRUE){
    spectrum <- spectrum |>
      spectrum_denoise(
        precursor = precursor_mz,
        n = n_expect)
  }

  return(spectrum)
}
