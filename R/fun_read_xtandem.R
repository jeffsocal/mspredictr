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

  d <- x |> xml2::read_html()
  h <- d |> xml2::xml_find_all('.//group[@id]', flatten = FALSE)

  xml_dataframe <- function(x){x |> unlist() |> t() |> as.data.frame()}

  l_peptides <- h |> xml2::xml_find_all('.//peptide/domain')

  tbl_seq <- list()
  n <- 0
  for(l_pep in l_peptides){
    n <- n + 1

    tbl_seq[[n]] <- l_pep |>
      xml2::xml_attrs() |>
      unlist() |> t() |>
      as.data.frame() |>
      tidyr::separate(id, into = c('spectrum_num', 'a' , 'b')) |>
      dplyr::rename(psm_sequence = seq) |>
      dplyr::mutate(psm_peptide = psm_sequence)

    tmp_mod <- l_pep |>
      xml2::xml_children() |>
      xml2::xml_attrs()

    if(length(tmp_mod) > 0){
      tbl_seq[[n]]$psm_peptide <- tbl_seq[[n]]$psm_peptide |> xtandem_peptide(tmp_mod, tbl_seq[[n]]$start)
    }

  }

  tbl_peptides <- tbl_seq |>
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


#' helper function to read in platform specific results
#'
#' @param sequence string
#' @param modifications string
#'
#' @return string
#'
xtandem_peptide <- function(
    sequence = NULL,
    modifications = NULL,
    val_start = 1
){

  sequence <- sequence |> stringr::str_extract_all('[A-Z]') |> unlist()

  val_masses <- modifications |> lapply(function(x){x[3] |> as.numeric()}) |> unlist() |> as.numeric()
  val_locate <- modifications |> lapply(function(x){x[2]}) |> unlist() |> as.numeric() - as.numeric(val_start) + 1

  out <- ''
  for(i in 1:length(sequence)){
    if(i %in% val_locate){
      # xtandem can place multiple modifications on the same residue
      res <- paste0("[", sequence[i], num_trunc(sum(val_masses[which(i == val_locate)]),2), "]")
      sequence[i] <- res
    }
    out <- paste0(out, sequence[i])
  }
  return(out)
}
