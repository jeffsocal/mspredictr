#' A helper function used to simulate the accuracy of synthetic spectra.
#'
#' @description
#' `sim_spectrum_accuracy()` Adds the noise to the accuracy of a synthetic spectrum.
#'
#' @param df
#' A dataframe of monoisotopic peaks.
#'
#' @param accuracy the targeted mass accuracy defined as the mean of the absolute
#' mass differences from predicted. Note, that a median mass measurement
#' accuracy (MMA) of 0.1Th will contain 95% of all peaks within +/-0.3Th, and
#' 99% of all peaks within +/-0.4Th, therefor with any given MMA, tolerances
#' should account for 3x or 4x, to capture 95% or 99% of all peaks, respectively.
#'
#' @param model
#' A the statistical sampling method to use of assessing the match accuracy
#'
sim_spectrum_accuracy <- function(
    df = NULL,
    accuracy = 0.1,
    model = stats::rnorm
){

  ui_mean <- 0
  ui_sd <- accuracy / 0.67

  df[,1] <- df[,1] + model(n = length(df[,1]),
                           mean = ui_mean,
                           sd = ui_sd)

  return(df)
}
