#' Simulate isotope peaks
#'
#' @description
#' `spectrum_isotopes()` is a helper function to generate isotope peaks in
#' synthetic spectra.
#'
#' @param spectrum
#' A dataframe of monoisotopic peaks
#'
#' @param isotope_model
#' A model of isotopes for a given polypeptide
#'
#' @param charge
#' The proton charge to assume for the isotope profile
#'
spectrum_isotopes <- function(
    spectrum = NULL,
    isotope_model = NULL,
    charge = 1
){

  # visible bindings
  mz <- NULL

  averagene <- 'N'
  if(is.null(isotope_model)) {
    isotope_model <- model_isotopes(averagene)
  }
  a_mass <- peptide_mass(averagene) - mass_atomic('H') * 2 - mass_atomic('0')
  n_mass <- mass_neutron()
  p_mass <- mass_proton()

  spectrum_new <- list()
  for(i in 1:nrow(spectrum)){
    index <- ceiling(spectrum[i,1] / a_mass)
    iso <- isotope_model[[index]]

    spectrum_new[[i]] <- data.frame(
      mz = (n_mass * (1:length(iso) - 1) + spectrum[i,1] + (p_mass * (charge - 1))) / charge,
      intensity = iso * spectrum[i,2]
    )
  }

  return(spectrum_new |> dplyr::bind_rows() |> dplyr::arrange(mz))
}

#' An isotope model using the BRAIN algorithm
#'
#' @description
#' `model_isotopes()` a helper function to generate isotopes
#'
#' @param averagene
#' The amino acid to use for a base polypeptide
#'
#' @param max_n
#' The maximum number of isotopes to account for
#'
model_isotopes <- function(
    averagene = c('N','A','R','D','C','E','Q','G','H','I',
                  'K','M','F','P','S','T','W','Y','V'),
    max_n = 42
){

  averagene <- rlang::arg_match(averagene)

  peaks <- list()
  for(i in 1:max_n){
    peaks[[i]] <- paste(rep(averagene, i), collapse = "") |>
      BRAIN::getAtomsFromSeq() |>
      BRAIN::calculateIsotopicProbabilities(nrPeaks = 12)

    peaks[[i]] <- peaks[[i]][which(peaks[[i]] > 0.01)]
  }

  return(peaks)
}
