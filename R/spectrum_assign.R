#' assign a table of fragment masses to a spectrum data object
#'
#' @description
#' `assign_spectrum()` converts an amino acid character to a numeric
#'
#' @param spectrum a spectrum data object
#' @param table_fragments a fragment data object
#' @param tolerance a numeric
#'
#' @return a tibble
#' @export
#'
#'
spectrum_assign <- function(
    spectrum = NULL,
    peptide = NULL,
    tolerance = 0.1,
    ...){


  if(is.null(spectrum)) { cli::cli_abort("`spectrum` must not be null")}
  if(is.null(peptide)) { cli::cli_abort("`peptide` must not be null")}

  # convert to a dataframe for tidyverse
  if(is.matrix(spectrum)) {
    spectrum <- spectrum |> as.data.frame()
  }

  table_fragments <- fragments(sequence = peptide, ...)
  # spectrum <- spc
  # table_fragments <- table_predicted %>% unite(ion, ion, z, sep="_")

  cn <- colnames(spectrum)
  w_mz <- which(grepl("^m", cn))
  w_int <- which(grepl("^i", cn))
  colnames(spectrum)[w_mz] <- 'mz'
  colnames(spectrum)[w_int] <- 'int'

  spectrum <- spectrum %>%
    dplyr::mutate(rid = dplyr::row_number() %>% as.character(),
           ion = rid)

  out <- pairwise_delta(spectrum, table_fragments, 'mz', 'ion') %>%
    dplyr::filter(abs(dif) <= tolerance) %>%
    dplyr::rename(rid = cluster_id,
                  ion = feature_id) %>%
    dplyr::inner_join(spectrum %>% dplyr::select(mz, int, rid), by='rid') %>%
    dplyr::inner_join(table_fragments %>% dplyr::select(-mz), by=c('ion')) %>%
    dplyr::select(
      mz, int, ion, type, z,
      error = dif, seq, pos, pair
    ) %>%
    dplyr::group_by(mz) %>%
    dplyr::slice_min(abs(error), n=1, with_ties = F) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(ion) %>%
    dplyr::slice_max(int, n=1, with_ties = F) %>%
    dplyr::ungroup()

  return(out)
}
