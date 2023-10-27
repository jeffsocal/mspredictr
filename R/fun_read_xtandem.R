#' helper function to read in platform specific results
#'
#' @param x location of file to parse
#'
#' @return a tibble
#'
read_xtandem <- function(
    x,
    cpus = 1
){

  # x <- "~/Local/data/results_hela_oe480/xtandem_decoy/Hela50ng_45min_DDA_1.ms2.t.xml"
  d <- x |> xml2::read_html()
  h <- d |> xml2::xml_find_all('.//group[@id]', flatten = FALSE)

  xml_dataframe <- function(x){x |> unlist() |> t() |> as.data.frame()}

  tbl_peptides <- h |>
    xml2::xml_find_all('.//peptide/domain') |>
    xml2::xml_attrs() |>
    lapply(xml_dataframe) |>
    dplyr::bind_rows() |>
    tidyr::separate(id, into = c('spectrum_num', 'a' , 'b')) |>
    dplyr::mutate(spectrum_num = spectrum_num |> as.numeric(),
                  psm_peptide = NA)

  tbl_proteins <- h |>
    xml2::xml_attrs() |>
    lapply(xml_dataframe) |>
    dplyr::bind_rows() |>
    dplyr::select(spectrum_num = id, psm_protein = label) |>
    dplyr::mutate(psm_protein = stringr::str_remove(psm_protein, "\\s.+")) |>
    dplyr::mutate(spectrum_num = spectrum_num |> as.numeric())


  tbl <- tbl_proteins |>
    dplyr::full_join(tbl_peptides, by = "spectrum_num") |>
    dplyr::rename(psm_sequence = seq) |>
    dplyr::mutate(psm_score = expect |> as.numeric() |> log10() * -1,
           psm_nmass = as.numeric(mh),
           psm_dp = as.numeric(y_ions) + as.numeric(b_ions) / stringr::str_length(psm_sequence),
           file = sub("\\..*", "", basename(x))) |>
    dplyr::group_by(spectrum_num) |>
    dplyr::arrange(desc(psm_score)) |>
    dplyr::mutate(psm_rank = dplyr::row_number()) |>
    dplyr::ungroup() |>
    dplyr::mutate(origin = 'xtandem') |>
    dplyr::select(dplyr::matches('file|spectrum_num|origin|decoy|psm'))

  return(tbl)
}
