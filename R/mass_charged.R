#' Calculate the m/z of a charged mass
#'
#' @description
#' `charge_mass()` get the mz of a mass with charge
#'
#' @param mass as integer
#' @param z an integer
#'
#' @return a vector
#' @export
#'
mass_charged <- function(mass, z){
  if(z == 0){ return(mass) }
  (mass / z) + mass_proton()
}
