#' The main function for parsing a fasta file
#'
#' @description
#' `read_mzml()` get the current regex
#'
#' @param file_path a character string of the path to the MGF formatted file
#' @return a tibble
#'
read_mzml <- function(
    file_path = NULL
){

  cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))
  if(is.null(file_path)) {cli::cli_abort(c("x" = "file_path is empty"))}
  if(!grepl("\\.mzML", file_path)) {
    file_ext <- stringr::str_extract(file_path, "\\..+$")
    cli::cli_abort(c("x" = "expected a {.emph .mzML} file, got {.emph {file_ext}}"))
  }

  cli::cli_progress_step("Parsing mzML file {basename(file_path)}")

  tryCatch({

    obj_mzml <- file_path |> mzR::openMSfile()

    ## Get the header a data frame of attributes
    tbl_hdr <- mzR::header(obj_mzml) |>
      dplyr::mutate(precursorCharge = ifelse(precursorCharge == 0, 2, precursorCharge)) |>
      tibble::as_tibble() |>
      dplyr::select(
        spectrum_num = seqNum,
        ms_event = acquisitionNum,
        ms_event_level = msLevel,
        ms_event_info = filterString,
        precursor_rt = retentionTime,
        precursor_mz = precursorMZ,
        precursor_z = precursorCharge,
        precursor_intensity = precursorIntensity,
        collision_energy = collisionEnergy,
        injection_time = injectionTime
      ) |>
      ## Get the spectra a list of mz and intensity
      dplyr::mutate(peaks = mzR::spectra(obj_mzml)) |>
      dplyr::filter(ms_event_level == 2) |>
      dplyr::mutate(precursor_mh = purrr::map2(precursor_mz, precursor_z, mass_neutral) |> unlist(),
                    file = sub("\\.mzML", "", basename(file_path))) |>
      dplyr::relocate(precursor_mh, .before = 'peaks') |>
      dplyr::relocate(file)

  }, error = function(err) {
    err = as.character(as.vector(err))
    cli::cli_process_failed()
    cli::cli_abort(err)
  })

  return(tbl_hdr)
}
