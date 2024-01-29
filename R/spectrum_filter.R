#' Return the top n values from a vector
#'
#' @description
#' `filter_topn()` get the top n values from a vector
#'
#' @param df dataframe
#' @param n n indexes
#'
#' @return vector
#'
filter_topn <- function(df, n) {
  if(!is.data.frame(df)) {cli::cli_abort("`df` not a data.frame object")}
  w_keep <- which(which_top_n(df[,2], n))
  return(df[w_keep,])
}


#' Return the top n values from a vector
#'
#' @description
#' `filter_topn()` get the top n values from a vector
#'
#' @param df dataframe
#' @param n n indexes
#'
#' @return vector
#'
filter_precursor <- function(df, mz) {
  if(!is.data.frame(df)) {cli::cli_abort("`df` not a data.frame object")}
  w_keep <- which_xprecursor(df[,1], mz)
  return(df[w_keep,])
}

#' Return a spectrum normalized according to xcorr publication
#'
#' @description
#' `filter_xcorr()` get the top n values from a vector
#'
#' @param df dataframe
#' @param chunk number of divisions to normalize across
#' @param norm the value to normalize to
#'
#' @return vector
#'
filter_xcorr <- function(df, chunk = 10, norm = 50) {
  # mz divided into 10 chunks, normalized to 50
  df |>
    dplyr::arrange(mz) |>
    split(1:chunk) |>
    lapply(function(x){ x[,2] <- (x[,2] / max(x[,2])) * norm; return(x) }) |>
    dplyr::bind_rows()
}
