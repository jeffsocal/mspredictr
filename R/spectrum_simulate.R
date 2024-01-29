#' Convert a peptide string to a named variable
#'
#' @description
#' `generate()` Generates a named variable
#'
#' @param data file path
#'
#' @return data object
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

  peaks_new <- list()
  for(z in charge){
    peaks_new[[z]] <- peptide |>
      fragments(charge = 1) |>
      dplyr::slice_head(prop = f_peaks) |>
      dplyr::mutate(intensity = 10^sample(100:360/100, dplyr::n())) |>
      dplyr::select(mz, intensity) |>
      spectrum_accuracy(accuracy = accuracy) |>
      spectrum_isotopes(charge = z)
  }

  peaks_new <- peaks_new |>
    dplyr::bind_rows() |>
    dplyr::mutate(intensity = (intensity / max(intensity)) * 1000) |>
    dplyr::mutate(intensity = intensity |> round()) |>
    spectrum_noise(noise_mean = noise_mean,
                   noise_quantile = noise_quantile,
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
