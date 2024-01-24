#' The main function for parsing a fasta file
#'
#' @description
#' `read_mgf()` get the current regex
#'
#' @param path a character string of the path to the MGF formatted file
#' @return a tibble
#' @export
#'
write_mgf <- function(
  data = NULL,
  path = NULL
){

  check_ms2spectra(data)

  cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))
  if(is.null(path)) {cli::cli_abort(c("x" = "path is empty"))}
  if(!grepl("\\.mgf", path)) {
    file_ext <- stringr::str_extract(path, "\\..+$")
    cli::cli_abort(c("x" = "expected a {.emph .mgf} file, got {.emph {file_ext}}"))
  }

  cli::cli_progress_step("Writing MGF file {basename(path)}")

  for(i in 1:nrow(data)){

    append <- FALSE
    if(i > 1) {append <- TRUE}

    out_string <- c("BEGIN IONS")
    out_string <- c(out_string, glue::glue("TITLE={data$ms_event_info[i]}"))
    out_string <- c(out_string, glue::glue("SCANS={data$ms_event[i]}"))
    out_string <- c(out_string, glue::glue("RTINSECONDS={data$precursor_rt[i]}"))
    out_string <- c(out_string, glue::glue("PEPMASS={data$precursor_mz[i]} {data$precursor_intensity[i]}"))
    out_string <- c(out_string, glue::glue("CHARGE=+{data$precursor_z[i]}"))
    for(ii in 1:dim(data$peaks[[i]])[1]){
      out_string <- c(out_string, glue::glue("{data$peaks[[i]][ii,1]}\t{round(data$peaks[[i]][ii,2])}"))
    }
    out_string <- c(out_string, "END IONS")
    out_string <- c(out_string, "")

    out_string |> readr::write_lines(path, append=append)
  }

}
