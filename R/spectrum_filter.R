#' Get the top most abundant peaks from a spectrum
#'
#' @description
#' `filter_topn()` is a helper function to get the top n values from a vector
#'
#' @param spectrum
#' The spectra data object
#'
#' @param n
#' The n number of top peaks to return
#'
filter_topn <- function(
    spectrum = NULL,
    n = 30
) {
  if(!is.data.frame(spectrum)) {cli::cli_abort("`spectrum` not a data.frame object")}
  w_keep <- which(which_top_n(spectrum[,2], n))
  return(spectrum[w_keep,])
}


#' Return the top n values from a vector
#'
#' @description
#' `filter_precursor()` is a helper function to remove residual precursor from
#' the spectrum
#'
#' @param spectrum
#' The spectra data object
#'
#' @param mz
#' The m/z value of the precursor to remove
#'
filter_precursor <- function(
    spectrum = NULL,
    mz = 30
) {
  if(!is.data.frame(spectrum)) {cli::cli_abort("`spectrum` not a data.frame object")}
  w_keep <- which_xprecursor(spectrum[,1], mz)
  return(spectrum[w_keep,])
}

#' Xcorr normalized spectrum
#'
#' @description
#' `filter_xcorr()` get the top n values from a vector
#'
#' @param spectrum
#' The spectra data object
#'
#' @param chunk
#' The number of divisions to normalize across
#'
#' @param norm
#' The value to normalize to
#'
#'
filter_xcorr <- function(
    spectrum = NULL,
    chunk = 10,
    norm = 50
) {
  # visible bindings
  mz <- NULL
  # mz divided into 10 chunks, normalized to 50
  spectrum |>
    dplyr::arrange(mz) |>
    split(1:chunk) |>
    lapply(function(x){ x[,2] <- (x[,2] / max(x[,2])) * norm; return(x) }) |>
    dplyr::bind_rows()
}
