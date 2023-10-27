#' truncate a number to a given decimal value
#'
#' @description
#' `num_trunc()` truncate a decimal to a given precision
#'
#' @param x the number to truncate
#' @param d decimals to truncate
#'
#' @return a numeric vector
#' @export
#'
num_trunc <- function(x, d = 2){floor(x * 10^d)/(10^d)}
