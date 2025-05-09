# Generated by extendr: Do not edit by hand

# nolint start

#
# This file was created with the following call:
#   .Call("wrap__make_mspredictr_wrappers", use_symbols = TRUE, package_name = "mspredictr")

#' @docType package
#' @usage NULL
#' @useDynLib mspredictr, .registration = TRUE
NULL

#' Return the mass of an atom
#' @param atom
#' The IUPAC amomic letter designation
#' @export
#' @examples
#' mass_atomic('C')
#' mass_atomic('Na')
mass_atomic <- function(atom) .Call(wrap__mass_atomic, atom)

#' Return the mass of a proton
#' @export
#' @examples
#' mass_proton()
mass_proton <- function() .Call(wrap__mass_proton)

#' Return the mass of a neutron
#' @export
#' @examples
#' mass_neutron()
mass_neutron <- function() .Call(wrap__mass_neutron)

#' Return the mass of a amino acid residue
#' @param s
#' The single letter designation for an amino acid
#' @export
#' @examples
#' mass_residue('G')
#' mass_residue('Y')
mass_residue <- function(s) .Call(wrap__mass_residue, s)

#' Return the charged mass
#' @param mass
#' The floating point mass value
#' @param z
#' The integer charge value
#' @export
#' @examples
#' mass_charged(1234.567, 2)
mass_charged <- function(mass, z) .Call(wrap__mass_charged, mass, z)

#' Return the neutral mass
#' @param mz
#' The floating point mass value
#' @param z
#' The integer charge value
#' @export
#' @examples
#' mass_neutral(1234.567, 2)
mass_neutral <- function(mz, z) .Call(wrap__mass_neutral, mz, z)

#' Return the mass of a peptide sequence
#' @param sequences
#' A vector of peptide sequence strings.
#' @export
#' @examples
#' peptide_mass(c('SAMPLE', 'SILLY'))
peptide_mass <- function(sequences) .Call(wrap__peptide_mass, sequences)

#' Return the unit length of a peptide sequence
#' @param sequences
#' A vector of peptide sequence strings.
#' @export
#' @examples
#' peptide_mass(c('SA[M15.99]PLE', 'SAM[15.99]PLE'))
peptide_length <- function(sequences) .Call(wrap__peptide_length, sequences)

#' Return the peptide sequence reversed
#' @param sequences
#' A vector of peptide sequence strings.
#' @export
#' @examples
#' peptide_mass(c('THINK', 'SAM[15.99]PLER'))
peptide_reverse <- function(sequences) .Call(wrap__peptide_reverse, sequences)

#' Return the vector of mass values for each peptide unit
#' @param seq
#' The peptide sequence string
#' @export
#' @examples
#' mass_ladder('SA[M15.99]PLER')
mass_ladder <- function(seq) .Call(wrap__mass_ladder, seq)

#' Return a simplified fragment mz vector
#' @param seq
#' The peptide sequence string
#' @export
#' @examples
#' mass_fragments('SA[M15.99]PLER')
mass_fragments <- function(seq) .Call(wrap__mass_fragments, seq)

#' Return the fragment indexes
#' @param seq
#' The peptide sequence string
#' @export
#' @examples
#' indef_fragments('SA[M15.99]PLER')
indef_fragments <- function(seq) .Call(wrap__indef_fragments, seq)

#' Return the fragment indexes
#' @param seq
#' The peptide sequence string
#' @param tolerance
#' The numerical float for mass tolerance in Th
#' @export
#' @examples
#' index_fragments('SA[M15.99]PLER', 0.05)
index_fragments <- function(seq, tolerance) .Call(wrap__index_fragments, seq, tolerance)

#' Return a mass indexes.
#' @param masses
#' A vector of mass sequences
#' @param tolerance
#' The numerical float for mass tolerance in Th
#' @export
#' @examples
#' index_mass(c(123.45, 234.56, 345.67), 0.05)
index_mass <- function(masses, tolerance) .Call(wrap__index_mass, masses, tolerance)

#' Return a I/L substituted sequence
#' @param s
#' The peptide sequence string.
#' @export
#' @examples
#' peptide_xleucine('SILLY')
peptide_xleucine <- function(s) .Call(wrap__peptide_xleucine, s)

#' Returns the boolean index of the top N largest values
#' @param f
#' A vector of numerical floats
#' @param n
#' The number of top values to keep
#' @export
#' @examples
#' which_top_n(c(123.45, 234.56, 345.67), 2)
which_top_n <- function(f, n) .Call(wrap__which_top_n, f, n)

#' Returns the boolean index of values that are in proximity to mz
#' @param f
#' A vector of numerical floats
#' @param mz
#' The the value to look for
#' @export
#' @examples
#' which_xprecursor(c(123.45, 234.56, 345.67), 233.61)
which_xprecursor <- function(f, mz) .Call(wrap__which_xprecursor, f, mz)

#' A Helper function that returns the array index of the monoisotopes given
#' group_isotopes()
#' @param vec_iso
#' A vector of numerical floats representing the isotopic groups
#' @export
#' @examples
#' group_isotopes(c(287.171, 288.119, 288.174, 290.161, 291.137,
#'                  291.164, 292.177, 293.124, 296.135, 298.139),
#'                c(218487, 44736, 29195, 1021168, 46029,
#'                  104552, 21997, 15262, 19908, 61741)) |>
#'   which_monoisotopes()
which_monoisotopes <- function(vec_iso) .Call(wrap__which_monoisotopes, vec_iso)

#' Returns the isotopic grouping for a given mass spectrum take in only mz, rt,
#' and abundance as vectors then build the Vec-IsotopeFeature- this gets around
#' the use of data.frame in R
#' @param vec_mz
#' A vector of numerical floats representing the mz component
#' (both vectors must be sorted on this value)
#' @param vec_int
#' A vector of numerical floats representing the ion intensity component
#' @export
#' @examples
#' group_isotopes(c(287.171, 288.119, 288.174, 290.161, 291.137,
#'                  291.164, 292.177, 293.124, 296.135, 298.139),
#'                c(218487, 44736, 29195, 1021168, 46029,
#'                  104552, 21997, 15262, 19908, 61741))
group_isotopes <- function(vec_mz, vec_int) .Call(wrap__group_isotopes, vec_mz, vec_int)

#' Helper function that provides the isotopic assignment given the output
#' from group_isotopes()
#' @param vec_iso
#' A vector of numerical floats representing the isotopic groups
#' @export
#' @examples
#' group_isotopes(c(287.171, 288.119, 288.174, 290.161, 291.137,
#'                  291.164, 292.177, 293.124, 296.135, 298.139),
#'                c(218487, 44736, 29195, 1021168, 46029,
#'                  104552, 21997, 15262, 19908, 61741)) |>
#'   label_isotopes()
label_isotopes <- function(vec_iso) .Call(wrap__label_isotopes, vec_iso)


# nolint end
