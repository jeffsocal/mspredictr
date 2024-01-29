#' Convert a peptide string to a named variable
#'
#' @description
#' `spectrum_isotopes()` Generates a named variable
#'
#' @param df a dataframe of monoisotopic peaks
#' @param isotope_model a model of isotopes for a given polypeptide
#' @param charge proton charge to assume for the isotope profile
#'
#' @return data object
#'
spectrum_isotopes <- function(
    df = NULL,
    isotope_model = NULL,
    charge = 1
){

  averagene <- 'N'
  if(is.null(isotope_model)) {
    isotope_model <- model_isotopes(averagene)
  }
  a_mass <- peptide_mass(averagene) - mass_atomic('H') * 2 - mass_atomic('0')
  n_mass <- mass_neutron()
  p_mass <- mass_proton()

  df_new <- list()
  for(i in 1:nrow(df)){
    index <- ceiling(df[i,1] / a_mass)
    iso <- isotope_model[[index]]

    df_new[[i]] <- data.frame(
      mz = (n_mass * (1:length(iso) - 1) + df[i,1] + (p_mass * (charge - 1))) / charge,
      intensity = iso * df[i,2]
    )
  }

  return(df_new |> dplyr::bind_rows() |> dplyr::arrange(mz))
}

#' Convert a peptide string to a named variable
#'
#' @description
#' `model_isotopes()` Generates a named variable
#'
#' @param averagene amino acid to use for a base polypeptide
#'
#' @return data list
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
