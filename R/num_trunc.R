#' Truncate a number to a given decimal value
#'
#' @description
#' `num_trunc()` truncates a float without rounding, floor or ceiling
#'
#' @param x
#' The input number, or number vector, to truncate. Values must be smaller than 1000.
#'
#' @param d
#' The number of decimals to truncate to.
#'
#' @export
#'
#' @examples
#' num_trunc(1.234567, 3)
#'
num_trunc <- function(
    x,
    d = 4
){
  if(x > 1000) return(x)
  floor(x * 10^d)/(10^d)
}
