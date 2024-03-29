#' GGplot2 object of a ms2 spectrum.
#'
#' @description
#' `plot_spectrum()` generates a plot of a given mass spectrum, w/o matching
#' peptide(s) sequences.
#'
#' @param spectrum
#' as data frame cols = `[mz, i]`
#'
#' @param filter
#' A boolean to determine if the spectrum should be filtered to remove residual
#' precursor peaks, isotopes peaks and any low-level noise.
#'
#' @param peptides
#' An optional vector of peptide sequences that will result in plotting each in
#' separate facets for comparison.
#'
#' @param label_size
#' A value used to increase or decrease the assigned labels.
#'
#' @param tolerance
#' The tolerance in Th to allow predicted fragments to match observed values.
#'
#' @param ...
#' Flow-through parameters to `assign_sepectrum()`
#'
#' @export
#'
#' @examples
#'  # using the supplied spectrum from the msreadr package
#'  library(msreadr)
#'  mzml <- path_to_example() |>
#'          read_spectra()
#'  mzml |>
#'    subset(spectrum_num == 1) |>
#'    plot_spectrum('HAVSEGTK')
#'
plot_spectrum <- function(
    spectrum = NULL,
    filter = FALSE,
    peptides = NULL,
    label_size = 3,
    tolerance = 0.1,
    ...
){

  # visible bindings
  mz <- int <- type <- ion <- NULL

  spectrum <- spectrum |> spectrum_extract(filter)

  if(!is.data.frame(spectrum)) { cli::cli_abort("spectrum is not a data table object") }

  cn <- colnames(spectrum)
  w_mz <- which(grepl("^m", cn))
  w_int <- which(grepl("^i", cn))
  colnames(spectrum)[w_mz] <- 'mz'
  colnames(spectrum)[w_int] <- 'int'

  plot <- spectrum |>
    ggplot2::ggplot(ggplot2::aes(mz, int)) +
    ggplot2::geom_segment(ggplot2::aes(xend = mz, yend = 0)) +
    ggplot2::geom_hline(yintercept = max(spectrum$int) * 1.1, color = NA) +
    ggplot2::theme_classic() +
    ggplot2::scale_y_continuous(n.breaks = 5)

  if(is.null(peptides)) { return(plot) }

  plots <- list()
  for(i in 1:length(peptides)){

    peptide <- peptides[i]
    if(!is.character(peptide)) { cli::cli_abort("`sequence` must be a character string") }

    tbl_assn <- spectrum_assign(spectrum, peptide,
                                tolerance = tolerance,
                                ...)
    tbl_poss <- peptide |> fragments(...)

    v_err <- 1 - abs(tbl_assn$error)/tolerance
    v_int <- tbl_assn$int / sum(spectrum$int)
    score_01 <- sum(v_err) * sum(v_int)

    hits <- nrow(tbl_assn |> dplyr::filter(type != 'precursor'))
    dotp <- round(hits / nrow(tbl_poss |> dplyr::filter(type != 'precursor')),2)

    plots[[i]] <- plot +
      ggplot2::geom_text(data = tbl_assn,
                         ggplot2::aes(label = ion, color = type),
                         size = label_size,
                         hjust = 0, vjust = -0.1) +
      ggplot2::labs(title = glue::glue("{peptide}")) +
      ggplot2::annotate('text', x = -Inf, y = Inf,
                        vjust = 1, hjust = -0.1, size = 3,
                        label = glue::glue("score:{round(score_01,2)} ndp:{dotp}  n:{hits}")) +
      ggplot2::theme(legend.position = 'none') +
      ggplot2::scale_color_manual(values = c(
        'y' = 'blue',
        'b' = 'red',
        'precursor' = 'forestgreen')) +
      ggplot2::scale_x_continuous(n.breaks = 10) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 10))

  }

  return(gridExtra::grid.arrange(grobs = plots))
}
