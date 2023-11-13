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
    peptides = NULL,
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
    ggplot2::geom_hline(yintercept = max(spectrum$int) * 1.1, color = NA) +
    ggplot2::theme_classic() +
    ggplot2::scale_y_continuous(n.breaks = 4)

  if(is.null(peptides)) { return(plot) }

  plots <- list()
  for(i in 1:length(peptides)){

    peptide <- peptides[i]
    if(!is.character(peptide)) { cli::cli_abort("`sequence` must be a character string") }

    tbl_assn <- assign_spectrum(spectrum, peptide,
                                tolerance = tolerance,
                                ...)
    tbl_poss <- peptide |> fragments(...)

    v_err <- 1 - abs(tbl_assn$error)/tolerance
    v_int <- tbl_assn$int / sum(spectrum$int)
    score_01 <- sum(v_err) * sum(v_int)
    score_02 <- sum(v_err * v_int)

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
                        label = glue::glue("s:{round(score_01,2)}  x:{round(score_02,4)}  dp:{dotp}  n:{hits}")) +
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
