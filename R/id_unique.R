#' create a unique id value
#'
#' @description
#' `id_unqiue()` converts an amino acid character to a numeric
#'
#' @param n a spectrum data object
#' @param length a fragment data object
#'
#' @return a character vector
#' @export
#'
#'
id_unqiue <- function(n = 1,
                      length = 5){

  if(!is.integer(n)) { cli::cli_abort("`n` must be an integer")}
  if(!is.integer(length)) { cli::cli_abort("`length` must be an integer")}

  pool <- c(letters, LETTERS, 0:9)

  res <- character(n) # pre-allocating vector is much faster than growing it
  for(i in seq(n)){
    this_res <- paste0(sample(pool, length, replace = TRUE), collapse = "")
    while(this_res %in% res){ # if there was a duplicate, redo
      this_res <- paste0(sample(pool, length, replace = TRUE), collapse = "")
    }
    res[i] <- paste0("X",this_res)
  }
  res
}
