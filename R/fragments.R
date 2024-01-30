#' Get a table of fragments based on a peptide sequence.
#'
#' @description
#' `fragments()` Generates a table of fragments from a given peptide character string.
#'
#' @param sequence
#' The character string representing a peptide, or poly amino acid. The canonical
#' 20 amino acids are encoded in and chemical modifications can be represented by
#' and floating point numerical value enclosed by square brackets. If a canonical
#' amino acid is also enclosed in the square brackets `[M15.99]` it is assumed that
#' the numerical value is in addition to the mass of the residue, and thus represents
#' a post-translational modification (PTM). Any value, in a square bracket is assumed
#' to represent a monomer in the chain `[57.021]`, and reflected in the fragment ion
#' series.
#'
#' Example: `SAMPLER`, `SA[M15.99]PLER`
#'
#' @param type
#' The ion fragmentation type to calculate the mass values for.
#'
#' Allowable: `y`, `b`
#'
#' Default: `c("y", "b")`
#'
#' @param charge
#' The specified fragment ion charges to calculate.
#'
#' Allowable: `1:4`
#'
#' Default: `c(1, 2, 3, 4)`
#'
#' @param loss
#' The specified neutral mass loss experienced in some mass spectrometers.
#'
#' Allowable: `water`, `amine`
#'
#' Default: none
#'
#' @return
#' The output is a tibble of predicted mass values based on the input parameters.
#' The table consists of the columns
#' ```
#'    `ion`: the ion type (y,b) and the charge represented by +
#'
#'    `mz`: the calculated m/z value of the ion
#'
#'    `z`: the charge state of the ion
#'
#'    `seq`: the partial sequence represented by the fragmentation
#'
#'    `pair`: a tag use to identify which ions originated from the y/b ion split
#'
#'    `pos`: the position at which the fragmentation occurred
#'
#'    `type`: the ion type
#' ```
#'
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
#' fragments("SAMPLER", charge = 1)
#'
fragments <- function(
    sequence = 'SAMPLE',
    type = c('y','b'),
    charge = c(1,2,3,4),
    loss = NULL
){

  # visible bindings
  ion <- mz <- NULL

  if(!is.character(sequence)) {cli::cli_abort("`sequence` must be a character string")}
  if(!is.numeric(charge)) { cli::cli_abort("`charge` must be an integer")}
  if(!is.null(loss)) { loss <- rlang::arg_match(loss, c('water', 'amine')) }
  type <- rlang::arg_match(type)
  charge <- match(charge, 1:4)

  ml <- mass_ladder_named(sequence)

  mass_isotope <- 1.0025
  mass_p <- mass_proton()
  mass_H <- mass_atomic('H')
  mass_O <- mass_atomic('O')
  mass_N <- mass_atomic('N')
  mass_C <- mass_atomic('C')
  mass_amo <- mass_H * 3 + mass_N
  mass_wtr <- mass_H * 2 + mass_O

  pep_mass <- sum(ml) + mass_H * 2 + mass_O

  if(is.null(loss)) { loss <- "" }

  t <- list()
  n <- length(ml)

  for(z in c(charge, max(charge) + 1)){
    z_ch <- paste(rep("+", z), collapse="")

    t[[z]] <- base::data.frame(
      ion = paste0('MH', z_ch),
      mz = pep_mass / z + mass_p,
      z = z,
      seq = paste(names(ml), collapse = ""),
      pair = 'p00')
  }

  for(z in charge){
    z_ch <- paste(rep("+", z), collapse="")

    for(i in 1:n){

      pair <- paste0('p', stringr::str_pad(i-1, 2, 'left', '0'))

      if(i > 1){
        y <- (sum(ml[i:n]) + mass_O + mass_H * 2) / z + mass_p
        y_seq <- paste(names(ml[i:n]), collapse = "")
        p <- n-i+1
        t[[paste0('y',i,z)]] <- base::data.frame(
          ion = paste0('y',p, z_ch), mz = y, z = z, seq = y_seq, pos=p, pair = pair)

        if (grepl('amine', loss) && grepl('[RKQN]', y_seq)){
          t[[paste0('y',i,z,'a')]] <- base::data.frame(
            ion = paste0('y',p, '-a', z_ch), mz = y-mass_amo, z = z, seq = y_seq, pos=p)
        }
        if (grepl('water', loss) && grepl('[STED]', y_seq)){
          t[[paste0('y',i,z,'w')]] <- base::data.frame(
            ion = paste0('y',p, '-w', z_ch), mz = y-mass_wtr, z = z, seq = y_seq, pos=p)
        }
      }


      if( i == n) {next}

      b <- sum(ml[1:i]) / z + mass_p
      b_seq <- paste(names(ml[1:i]), collapse = "")
      p <- i
      pair <- paste0('p', stringr::str_pad(i, 2, 'left', '0'))
      t[[paste0('b',i,z)]] <- base::data.frame(
        ion = paste0('b',p, z_ch), mz = b, z = z, seq = b_seq, pos=p, pair = pair)

      if (grepl('amine', loss) && grepl('[RKQN]', b_seq)){
        t[[paste0('b',i,z,'a')]] <- base::data.frame(
          ion = paste0('b',p,'-a', z_ch), mz = b-mass_amo, z = z, seq = b_seq, pos=p)
      }

      if (grepl('water', loss) && grepl('[STED]', b_seq)){
        t[[paste0('b',i,z,'w')]] <- base::data.frame(
          ion = paste0('b',p,'-w', z_ch), mz = b-mass_wtr, z = z, seq = b_seq, pos=p)
      }

    }
  }

  out <- dplyr::bind_rows(t) %>%
    dplyr::arrange(mz) %>%
    dplyr::mutate(type = stringr::str_extract(ion, "^."),
           type = ifelse(type == 'M', 'precursor', type))

  return(out)

}
