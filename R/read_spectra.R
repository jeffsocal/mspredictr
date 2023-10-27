#' Table of spectrum values
#'
#' @description
#' `read_spectra()` get search results from a variety of platforms
#'
#' @param path location of file to parse
#' @param include_spectra boolean to keep/drop the nested MSn spectrum
#'
#' @return a tibble
#'
read_spectra <- function(
    path = NULL,
    include_spectra = TRUE
){

  cli::cli_div(theme = list(span.info = list(color = "#ff4500")))
  if(is.null(path)){ cli::cli_abort("no path to a mzML file given") }
  if(!is.character(path)) { cli::cli_abort("`path` must be a character string")}
  if(!file.exists(path)) { cli::cli_abort("Not Found! `config:` {path}")}
  if(!is.logical(include_spectra)) { cli::cli_abort("`include_spectra` is not a boolean")}

  ######################################################################
  # read in the data file
  # - this allows us to cut peptides out of our search window based on what was collected
  obj_mzml <- path |> mzR::openMSfile()

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
    dplyr::mutate(precursor_nm = purrr::map2(precursor_mz, precursor_z, rmstandem::mass_neutral) |> unlist(),
                  file = sub("\\.mzML", "", basename(path))) |>
    dplyr::relocate(precursor_nm, .before = 'peaks') |>
    dplyr::relocate(file)

  if(include_spectra == FALSE){
    tbl_hdr <- tbl_hdr |> dplyr::select(-peaks)
  }

  return(tbl_hdr)
}
