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
assign_spectrum <- function(
    spectrum = NULL,
    table_fragments = NULL,
    tolerance = 0.1){

  if(!is.numeric(tolerance)) { cli::cli_abort("`tolerance` must be a numeric")}

  # spectrum <- spc
  # table_fragments <- table_predicted %>% unite(ion, ion, z, sep="_")

  spectrum <- spectrum %>%
    mutate(rid = row_number() %>% as.character(),
           ion = rid)

  out <- pairwise_delta(spectrum, table_fragments, 'mz', 'ion') %>%
    filter(abs(dif) <= tolerance) %>%
    dplyr::rename(rid = cluster_id,
                  ion = feature_id) %>%
    inner_join(spectrum %>% select(mz, i, rid), by='rid') %>%
    inner_join(table_fragments %>% select(-mz), by=c('ion')) %>%
    select(
      mz, i, ion, type, z,
      error = dif, seq, pos
    ) %>%
    group_by(mz) %>%
    slice_min(abs(error), n=1, with_ties = F) %>%
    ungroup() %>%
    group_by(ion) %>%
    slice_max(i, n=1, with_ties = F) %>%
    ungroup()

  return(out)
}
