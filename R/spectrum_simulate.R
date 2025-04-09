#' Create a synthetic fragment mass spectrum.
#'
#' @description
#' `spectrum_simulate()` generates a synthetic fragment mass spectrum, based on
#' a peptide string, using normally-distributed random processes.
#'
#' @param peptide
#' The peptide string to simulate a spectrum for
#'
#' @param charge
#' The charges expected in the fragment mass spectrum
#'
#' @param accuracy
#' The floating point mass value to use for calculating mass accuracy
#'
#' @param noise_mean
#' The mean of the normal distribution to sample from
#'
#' @param noise_quantile
#' The Nth quantile used to assess noise
#'
#' @param f_peaks
#' The proportion of peaks to keep. However, if a peak is randomly below the noise
#' threashold then it will be dropped.
#'
#' @export
#'
spectrum_simulate <- function(
    peptide = NULL,
    charge = 1:2,
    accuracy = 0.1,
    noise_mean = NULL,
    noise_quantile = .33,
    f_peaks = 0.67
){

  # visible bindings
  mz <- intensity <- NULL

  peaks_new <- list()
  for(z in charge){
    peaks_new[[z]] <- peptide |>
      fragments(charge = 1) |>
      dplyr::filter(z == 1) |>
      dplyr::slice_sample(prop = f_peaks) |>
      dplyr::mutate(intensity = 10^sample(100:360/100, dplyr::n())) |>
      dplyr::select(mz, intensity) |>
      sim_spectrum_accuracy(accuracy = accuracy) |>
      sim_spectrum_isotopes(charge = z)
  }

  peaks_new <- peaks_new |>
    dplyr::bind_rows() |>
    dplyr::mutate(intensity = (intensity / max(intensity)) * 1000) |>
    dplyr::mutate(intensity = intensity |> round()) |>
    sim_spectrum_noise(noise_mean = noise_mean,
                   noise_quantile = noise_quantile,
                   peak_range = c(100, peptide_mass(peptide)),
                   n_peaks = 600) |>
    filter_topn(300)

  out <-
    list(
      file = 'synthetic',
      spectrum_num = 0,
      ms_event = 0,
      ms_event_level = 2,
      ms_event_info = glue::glue("synthetic spectrum from {}"),
      precursor_rt = 0,
      precursor_mz = peptide_mass(peptide) |> mass_charged(max(charge) + 1),
      precursor_z = max(charge) + 1,
      precursor_intensity = 0,
      collision_energy = 28,
      injection_time = 100,
      precursor_mh = peptide_mass(peptide) |> mass_charged(1),
      peaks = list(peaks_new)
    )

  return(out)
}
