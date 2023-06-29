#' GGplot2 object of a ms2 spectrum
#'
#' @description
#' `plot_spectrum()` get the mass of a poly amino acid
#'
#' @param spectrum as data frame cols = `[mz, i]`
#' @param sequence as character string
#' @param charge as integer
#'
#' @return a string
#' @export
#'
plot_spectrum <- function(
    spectrum = NULL,
    sequence = NULL,
    charge = 1,
    label_size = 3,
    tolerance = 0.1
){

  if(!is.data.frame(spectrum)) { cli::cli_abort("spectrum is not a data table object") }
  if(!is.numeric(charge)) { cli::cli_abort("`charge` must be an integer") }
  charge <- intersect(charge, 1:50)[1]


  plot <- spectrum |>
    ggplot2::ggplot(ggplot2::aes(mz, i)) +
    ggplot2::geom_segment(ggplot2::aes(xend = mz, yend = 0)) +
    ggplot2::theme_classic()

  if(is.null(sequence)) { return(plot) }

  if(!is.character(sequence)) { cli::cli_abort("`sequence` must be a character string") }
  tbl_frag <- fragments(sequence, charge = 1:charge)
  tbl_assn <- assign_spectrum(spectrum, tbl_frag, tolerance = tolerance)

  plot <- plot +
    ggplot2::geom_text(data = tbl_assn,
                       ggplot2::aes(label = ion, color = type),
                       size = label_size,
                       hjust = 0, vjust = 0) +
    ggplot2::labs(title = sequence) +
    ggplot2::theme(legend.position = 'none') +
    ggplot2::scale_color_manual(values = c(
      'y' = 'blue',
      'b' = 'red',
      'precursor' = 'forestgreen'))

  return(plot)
}
