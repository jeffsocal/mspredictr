#' assign a table of fragment masses to a spectrum data object
#'
#' @description
#' `assign_string()` converts an amino acid character to a numeric
#'
#' @param spectrum a spectrum data object
#' @param table_fragments a fragment data object
#' @param tolerance a numeric
#'
#' @return a list
#' @export
#'
#'
assign_string <- function(
    spectrum = NULL,
    table_fragments = NULL,
    tolerance = 0.1){


  if(!is.numeric(tolerance)) { cli::cli_abort("`tolerance` must be a numeric")}

  tbl_ass <- assignSpectrum(spectrum, table_fragments, tolerance)

  sequence <- table_fragments %>%
    filter(type == 'precursor') %>%
    select(seq) %>%
    unlist() %>% as.character()

  seql <- massladder(sequence)
  seqn <- length(seql)

  o <- c()
  for(s in names(seql)){
    if(grepl("[0-9]", s)){
      sn <- str_extract(s, "[a-zA-Z]|[\\-\\+]*[0-9]+\\.*[0-9]*")
      s <- str_to_lower(sn[[1]][1])
    }

    o <- c(o, s, ' ')
  }
  o <- o[1:(length(o)-1)]
  b <- y <- gsub("[a-zA-Z]", " ", o)

  tbl_ass <- tbl_ass %>%
    filter(!grepl("a|w|MH", ion))

  ys <- tbl_ass$pos[tbl_ass$type == 'y'] %>% unique() %>% sort() * 2
  ys <- rep(seqn*2, length(ys)) - ys
  bs <- tbl_ass$pos[tbl_ass$type == 'b'] %>% unique() %>% sort() * 2

  b[bs] <- "|"
  y[ys] <- "|"

  out <- list(
    sequence = paste(o, collapse=''),
    b_ion = paste(b, collapse=''),
    y_ion = paste(y, collapse='')
  )

  return(out)
}
