#' Match peptide predicted fragments with a spectrum.
#'
#' @description
#' `assign_spectrum()` assigns the predicted fragment masses of a given peptide
#' sequence with the mass spectrum.
#'
#' @param spectrum
#' A spectrum data object.
#'
#' @param peptide
#' The character string representing a peptide, or poly amino acid. The canonical
#' 20 amino acids are encoded in and chemical modifications can be represented by
#' and floating point numerical value enclosed by square brackets. If a canonical
#' amino acid is also enclosed in the square brackets `[M15.99]` it is assumed that
#' the numerical value is in addition to the mass of the residue, and thus represents
#' a post-translational modification (PTM).
#'
#' @param tolerance
#' The tolerance in Th to allow predicted fragments to match observed values.
#'
#' @param ...
#' Flow-through parameters to `fragments()`
#'
#' @export
#'
#' @examples
#'  # using the supplied spectrum from the msreadr package
#'  library(msreadr)
#'  mzml <- path_to_example() |>
#'          read_spectra()
#'  mzml |>
#'    subset(spectrum_num == 1) |>
#'    spectrum_extract() |>
#'    spectrum_assign(peptide = 'HAVSEGTK')
#'
spectrum_assign <- function(
    spectrum = NULL,
    peptide = NULL,
    tolerance = 0.1,
    ...){

  # visible bindings
  rid <- NULL
  dif <- NULL
  cluster_id <- NULL
  feature_id <- NULL
  mz <- NULL
  int <- NULL
  ion <- NULL
  type <- NULL
  z <- NULL
  pos <- NULL
  pair <- NULL
  error <- NULL

  intensity <- NULL
  isotope_id <- NULL
  isotope_num <- NULL
  isotope_z <- NULL
  isotope_num <- NULL


  if(is.null(spectrum)) { cli::cli_abort("`spectrum` must not be null")}
  if(is.null(peptide)) { cli::cli_abort("`peptide` must not be null")}

  # convert to a dataframe for tidyverse
  if(is.matrix(spectrum)) {
    spectrum <- spectrum |> as.data.frame()
  }

  table_fragments <- fragments(sequence = peptide, ...)
  # spectrum <- spc
  # table_fragments <- table_predicted |> unite(ion, ion, z, sep="_")

  cn <- colnames(spectrum)
  w_mz <- which(grepl("^m", cn))
  w_int <- which(grepl("^int", cn))
  colnames(spectrum)[w_mz] <- 'mz'
  colnames(spectrum)[w_int] <- 'intensity'

  intensity_median <- median(spectrum$intensity, na.rm = TRUE)

  spectrum <- spectrum |>
    dplyr::mutate(
      rid = dplyr::row_number() |> as.character(),
      ion = rid)

  out <- pairwise_delta(spectrum |> dplyr::select(mz, intensity, rid, ion),
                        table_fragments, 'mz', 'ion') |>
    dplyr::filter(abs(dif) <= tolerance) |>
    dplyr::rename(rid = cluster_id,
                  ion = feature_id) |>
    dplyr::inner_join(spectrum |> dplyr::select(-ion), by='rid') |>
    dplyr::inner_join(table_fragments |> dplyr::select(-mz), by=c('ion')) |>
    dplyr::rename(error = dif) |>
    dplyr::group_by(mz) |>
    dplyr::slice_min(abs(error), n=1, with_ties = F) |>
    dplyr::ungroup() |>
    dplyr::group_by(ion) |>
    dplyr::slice_max(intensity, n=1, with_ties = F) |>
    dplyr::ungroup() |>
    dplyr::select(mz, intensity, isotope_id, isotope_num, isotope_z,
                  seq, ion, z, pair, pos, type, error) |>
    dplyr::mutate(ion_score = (intensity / intensity_median) * (1 - abs(error)/tolerance))

  return(out)
}
