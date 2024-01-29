#' helper function to read in platform specific results
#'
#' @param x location of file to parse
#'
#' @return a tibble
#'
spectrum_extract <- function(
    x,
    filter = FALSE
){
  spectra <- x$peaks
  precursor_mz <- x$precursor_mz
  precursor_z <- x$precursor_z
  precursor_nm <- mass_neutral(precursor_mz, precursor_z)

  n_expect <- round(precursor_nm / 114.35) * 3
  spectra <- spectra |> as.data.frame()
  colnames(spectra) <- c('mz', 'intensity')

  if(filter == TRUE){
    spectra <- spectra |>
      spectrum_denoise(
        precursor = precursor_mz,
        n = n_expect)
  }

  return(spectra)
}
