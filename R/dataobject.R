#' FASTA data object
#'
#' @param obj FASTA data list
#'
#' @return FASTA data object
#'
ms2spectra <- function(obj) {
  class(obj) <- "ms2spectra"
  return(obj)
}

#' Check the integrity of a msfastar data object
#'
#' @description
#' `check_fasta()` is a helper function that checks the structure and contents of
#' a msfastar data object
#'
#' @param data msfastar data object
#'
#' @return silent on success, an abort message on fail
#'
check_ms2spectra <- function(
    x = NULL
){

  # fail if x is NULL
  if(is.null(x)) {
    cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))
    cli::cli_abort(c("x" = "Input cannot be {.emph NULL}"))
  }
  if(mode(x) != "list") {
    cli::cli_abort(c("x" = "Input is {.emph mode(x)}}, should be an {.emph list}"))
  }
  if(class(x) != 'ms2spectra') {
    cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))
    cli::cli_abort(c("x" = "Input must be of type {.emph ms2spectra}"))
  }

}
