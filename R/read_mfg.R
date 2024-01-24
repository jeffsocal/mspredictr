#' The main function for parsing a fasta file
#'
#' @description
#' `read_mgf()` get the current regex
#'
#' @param file_path a character string of the path to the MGF formatted file
#' @return a tibble
#'
read_mgf <- function(
    file_path = NULL
){

  cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))
  if(is.null(file_path)) {cli::cli_abort(c("x" = "file_path is empty"))}
  if(!grepl("\\.mgf", file_path)) {
    file_ext <- stringr::str_extract(file_path, "\\..+$")
    cli::cli_abort(c("x" = "expected a {.emph .mgf} file, got {.emph {file_ext}}"))
  }

  cli::cli_progress_step("Parsing MGF file {basename(file_path)}")

  tryCatch({

    l_mgf <- readr::read_file(file_path)
    # read in fasta file
    l_mgf <- unlist(base::strsplit(l_mgf, "BEGIN IONS"))
    l_mgf <- l_mgf[-1] # first in list is blank

    # $ file               : chr [1:8] "example.ms2" "example.ms2" "example.ms2" "example.ms2" ...
    # $ spectrum_num       : int [1:8] 1 2 3 4 5 6 7 8
    # $ ms_event           : int [1:8] 1 2 3 4 5 6 7 8
    # $ ms_event_level     : int [1:8] 2 2 2 2 2 2 2 2
    # $ ms_event_info      : chr [1:8] "FTMS + c NSI d Full ms2 414.7148@hcd28.00 [100.0000-865.0000]" "FTMS + c NSI d Full ms2 384.7073@hcd28.00 [100.0000-805.0000]" "FTMS + c NSI d Full ms2 422.7475@hcd28.00 [100.0000-880.0000]" "FTMS + c NSI d Full ms2 427.2117@hcd28.00 [100.0000-890.0000]" ...
    # $ precursor_rt       : num [1:8] 0.29 0.369 0.447 0.526 0.605 ...
    # $ precursor_mz       : num [1:8] 415 385 423 427 443 ...
    # $ precursor_z        : int [1:8] 2 2 2 2 2 2 2 2
    # $ precursor_intensity: num [1:8] 2.05e+08 2.27e+07 1.18e+07 6.82e+06 1.24e+07 ...
    # $ collision_energy   : num [1:8] 28 28 28 28 28 28 28 28
    # $ injection_time     : num [1:8] 100 100 100 100 100 100 100 100
    # $ precursor_mh       : num [1:8] 828 768 844 853 884 ...
    # $ peaks              :List of 8

    file <- basename(file_path)
    spectrum_num <- 1:length(l_mgf)
    ms_event_level <- 2

    l_mgf <- lapply(l_mgf, function(x){

      out <- c(
        'ms_event_info' = NA,
        'precursor_rt' = NA,
        'precursor_mz' = NA,
        'precursor_z' = NA,
        'precursor_intensity' = NA,
        'collision_energy' = NA,
        'injection_time' = NA,
        'precursor_mh' = NA
      )

      peaks <- matrix(nrow=0, ncol=2)

      lines <- base::strsplit(x, "\n") |> unlist()

      for(i in 1:length(lines)){
        line <- lines[[i]]
        if(grepl("\\=", line)){
          attr <- base::strsplit(line, "=") |> unlist()
          if(attr[1] == 'TITLE'){out['ms_event_info'] <- attr[2]}
          if(attr[1] == 'SCANS'){out['ms_event'] <- attr[2]}
          if(attr[1] == 'RTINSECONDS'){out['precursor_rt'] <- attr[2]}
          if(attr[1] == 'PEPMASS'){
            pep <- stringr::str_split(attr[2], "\\s+") |> unlist()
            out['precursor_mz'] <- pep[1]
            if(length(pep) > 1) {out['precursor_intensity'] <- pep[2]}
          }
          if(attr[1] == 'CHARGE'){
            attr[2] <- gsub("\\+|\\,", "", attr[2]) |> stringr::str_split("\\s+")
            out['precursor_z'] <- attr[2][1]
          }
        } else if(line == '' || grepl("END", line)){
          # do nothing
        } else {
          mzi <- stringr::str_split(line, "\\s+") |> unlist() |> as.numeric()
          peaks <- peaks |> rbind(mzi)
        }
      }
      out['precursor_mh'] <- mass_neutral(out['precursor_mz'] |> as.numeric(),
                                          out['precursor_z'] |> as.numeric())
      # clean up the peaks matrix
      dimnames(peaks) <- NULL
      colnames(peaks) <- c('mz', 'intensity')
      # return the data object
      return(
        out |>
          tibble::as_tibble() |>
          dplyr::mutate(peaks = list(peaks))
      )
    }) |>
      dplyr::bind_rows() |>
      dplyr::mutate(
        file = file,
        spectrum_num = spectrum_num,
        ms_event_level = ms_event_level
      ) |>
      dplyr::mutate(
        ms_event = ms_event |> as.numeric(),
        precursor_rt = precursor_rt |> as.numeric(),
        precursor_mz = precursor_mz |> as.numeric(),
        precursor_z = precursor_z |> as.numeric(),
        precursor_intensity = precursor_intensity |> as.numeric(),
        collision_energy = collision_energy |> as.numeric(),
        injection_time = injection_time |> as.numeric()
      ) |>
      dplyr::select(
        file, spectrum_num, ms_event,
        ms_event_level, ms_event_info,
        precursor_rt, precursor_mz,
        precursor_z, precursor_intensity,
        collision_energy, injection_time,
        precursor_mh, peaks
      )

  }, error = function(err) {
    err = as.character(as.vector(err))
    cli::cli_process_failed()
    cli::cli_abort(err)
  })

  class(l_mgf) <- 'ms2spectra'

  return(l_mgf)
}
