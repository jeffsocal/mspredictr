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
#' @param rm_precursor
#' A boolean to remove the precursor.
#'
#' @param rm_isotopes
#' A boolean to remove isotopes.
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
    peptides = NULL,
    rm_precursor = FALSE,
    rm_isotopes = FALSE,
    label_size = 3,
    tolerance = 0.1,
    ...
){

  # visible bindings
  mz <- int <- type <- ion <- NULL

  spectrum <- spectrum |> spectrum_extract(rm_precursor, rm_isotopes)

  if(!is.data.frame(spectrum)) { cli::cli_abort("spectrum is not a data table object") }

  cn <- colnames(spectrum)
  w_mz <- which(grepl("^m", cn))
  w_int <- which(grepl("^int", cn))
  colnames(spectrum)[w_mz] <- 'mz'
  colnames(spectrum)[w_int] <- 'intensity'

  plot <- spectrum |>
    ggplot2::ggplot(ggplot2::aes(mz, intensity)) +
    ggplot2::geom_segment(ggplot2::aes(xend = mz, yend = 0)) +
    ggplot2::geom_hline(yintercept = max(spectrum$intensity) * 1.1, color = NA) +
    ggplot2::theme_classic() +
    ggplot2::scale_y_continuous(n.breaks = 5)

  if(is.null(peptides)) { return(plot) }

  tbl_assn <- list()
  tbl_score <- list()
  for(i in 1:length(peptides)){

    peptide <- peptides[i]
    if(!is.character(peptide)) { cli::cli_abort("`sequence` must be a character string") }

    tbl_assn[[i]] <- spectrum |>
      spectrum_assign(peptide,
                      tolerance = tolerance,
                      ...) |>
      dplyr::mutate(peptide = peptide)

    tbl_poss <- peptide |> fragments(...)

    # v_err <- 1 - abs(tbl_assn[[i]]$error)/tolerance
    # v_int <- tbl_assn[[i]]$intensity / sum(spectrum$intensity)
    # v_int_a <- tbl_assn[[i]]$intensity / median(spectrum$intensity, na.rm = TRUE)
    #
    # score <- sum(v_err) * sum(v_int)
    # score_a <- sum(v_int_a * v_err)
    score_b <- sum(tbl_assn[[i]]$ion_score)
    hits <- nrow(tbl_assn[[i]] |> dplyr::filter(type != 'precursor'))
    dotp <- round(hits / nrow(tbl_poss |> dplyr::filter(type != 'precursor')),2)

    tbl_score[[i]] <- data.frame(
      score = score_b,
      hits = hits,
      dotp = dotp,
      peptide = peptide
    )
  }

  tbl_assn <- tbl_assn |> dplyr::bind_rows()
  tbl_score <- tbl_score |> dplyr::bind_rows()

  plot <- plot +
    ggplot2::geom_text(data = tbl_assn,
                       ggplot2::aes(label = ion, color = type),
                       size = label_size,
                       hjust = 0, vjust = -0.1) +
    ggplot2::geom_text(data = tbl_score,
                       x = -Inf, y = Inf,
                       vjust = 1, hjust = -0.1, size = 3,
                       ggplot2::aes(label = glue::glue("score:{round(score,2)} ndp:{dotp}  n:{hits}"))) +
    ggplot2::theme(legend.position = 'none') +
    ggplot2::scale_color_manual(values = c(
      'y' = 'blue',
      'b' = 'red',
      'precursor' = 'forestgreen')) +
    ggplot2::scale_x_continuous(n.breaks = 10) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10)) +
    ggplot2::facet_wrap(~peptide)

  return(plot)
}
