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
    peptide = NULL,
    label_size = 3,
    tolerance = 0.1,
    ...
){

  if(!is.data.frame(spectrum)) { cli::cli_abort("spectrum is not a data table object") }

  cn <- colnames(spectrum)
  w_mz <- which(grepl("^m", cn))
  w_int <- which(grepl("^i", cn))
  colnames(spectrum)[w_mz] <- 'mz'
  colnames(spectrum)[w_int] <- 'int'

  plot <- spectrum |>
    ggplot2::ggplot(ggplot2::aes(mz, int)) +
    ggplot2::geom_segment(ggplot2::aes(xend = mz, yend = 0)) +
    ggplot2::theme_classic()

  if(is.null(peptide)) { return(plot) }

  if(!is.character(peptide)) { cli::cli_abort("`sequence` must be a character string") }
  tbl_assn <- assign_spectrum(spectrum, peptide,
                              tolerance = tolerance,
                              ...)

  plot <- plot +
    ggplot2::geom_text(data = tbl_assn,
                       ggplot2::aes(label = ion, color = type),
                       size = label_size,
                       hjust = 0, vjust = 0) +
    ggplot2::labs(title = peptide) +
    ggplot2::theme(legend.position = 'none') +
    ggplot2::scale_color_manual(values = c(
      'y' = 'blue',
      'b' = 'red',
      'precursor' = 'forestgreen')) +
    ggplot2::scale_x_continuous(n.breaks = 7)

  return(plot)
}
