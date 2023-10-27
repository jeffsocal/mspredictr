#' Table of PSMs results
#'
#' @description
#' `read_psms()` get search results from a variety of platforms
#'
#' @param path location of file to parse
#' @param platform as character string
#'
#' @return a tibble
#' @export
#'
read_psms <- function(
    path_to_psms = NULL,
    path_to_mzml = NULL,
    platform = c('comet','ms_amanda','omssa','sage','xtandem'),
    cpus = 4
){

  if(is.null(path_to_psms)){ cli::cli_abort("no path to PSM file given") }
  if(!file.exists(path_to_psms)){ cli::cli_abort("PSM file does not exist") }

  if(!is.null(path_to_mzml)){
    if(!file.exists(path_to_mzml)){ cli::cli_abort("mzML file does not exist") }
    spec <- path_to_mzml |> read_spectra(include_spectra = FALSE)
  }
  if(!cpus %in% 1:parallel::detectCores()){
    cli::cli_abort("number of cores outside the limit of 1|{parallel::detectCores()}")
  }
  rlang::arg_match(platform)

  switch (platform,
          comet = {
            out <- read_comet(path_to_psms, cpus)
          },
          tide = {
            out <- read_tide(path_to_psms, cpus)
          },
          ms_amanda = {
            out <- read_ms_amanda(path_to_psms, cpus)
          },
          omssa = {
            out <- read_omssa(path_to_psms, cpus)
          },
          sage = {
            out <- read_sage(path_to_psms, cpus)
          },
          xtandem = {
            out <- read_xtandem(path_to_psms, cpus)
          },
          {
            return(NULL)
          }
  )

  # add file meta data and down-select columns
  out <- out |>
    dplyr::mutate(platform = platform,
                  file = path_to_psms |> basename() |> stringr::str_remove("\\..+")) |>
    dplyr::select(dplyr::matches('file|platform|decoy|psm|ms_event|spectrum_num'))

  # order the columns
  out <- out |> dplyr::relocate(sort(colnames(out)))

  # add in scan level information if present
  if(!is.null(path_to_mzml)){
    switch (platform,
            xtandem = {
              out <- spec |> dplyr::left_join(out, by = 'spectrum_num')
            },
            {
              out <- spec |> dplyr::left_join(out, by = 'ms_event')
            }
    )
  }

  return(out)
}
