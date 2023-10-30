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

  x <- "~/Local/data/project_jester/results/simulation_ecoli_10kPeptides/parse_test/peptides10000_simPTM_ecoli-BL21DE3-4156_1.t.xml"
  d <- x |> xml2::read_html()
  h <- d |> xml2::xml_find_all('.//group[@id]', flatten = FALSE)

  xml_dataframe <- function(x){x |> unlist() |> t() |> as.data.frame()}

  l_peptides <- h |>
    xml2::xml_find_all('.//peptide/domain')

  tbl_seq <- list()
  tbl_mod <- list()
  n <- 0
  for(l_pep in l_peptides){
    n <- n + 1

    tbl_seq[[n]] <- l_pep |>
      xml2::xml_attrs() |>
      unlist() |> t() |> as.data.frame()

    l_mod <- l_pep |>
      xml2::xml_children() |>
      xml2::xml_attrs() |>
      unlist()

    if(length(l_mod) > 0){
      tbl_mod[[n]] <- l_mod |> t() |> as.data.frame()
    }


    # lapply(xml_dataframe) |>
    # dplyr::bind_rows() |>
    # tidyr::separate(id, into = c('spectrum_num', 'a' , 'b')) |>
    # dplyr::mutate(spectrum_num = spectrum_num |> as.numeric(),
    #               psm_peptide = NA) |>
    # tibble::as_tibble()
  }

  |>


    tbl_modifications <- h |>
    xml2::xml_find_all('.//peptide/domain/aa') |>
    xml2::xml_attrs() |>
    lapply(xml_dataframe) |>
    dplyr::bind_rows() |>
    tibble::as_tibble()

  tbl_proteins <- h |>
    xml2::xml_attrs() |>
    lapply(xml_dataframe) |>
    dplyr::bind_rows() |>
    dplyr::select(spectrum_num = id, psm_protein = label) |>
    dplyr::mutate(psm_protein = stringr::str_remove(psm_protein, "\\s.+")) |>
    dplyr::mutate(spectrum_num = spectrum_num |> as.numeric()) |>
    tibble::as_tibble()


  tbl <- tbl_proteins |>
    dplyr::full_join(tbl_peptides, by = "spectrum_num") |>
    dplyr::rename(psm_sequence = seq) |>
    dplyr::mutate(psm_score = expect |> as.numeric() |> log10() * -1,
                  psm_mh = as.numeric(mh),
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
