#' Get the mass of an amino acid
#'
#' @description
#' `pairwise_delta()` converts an amino acid character to a numeric
#'
#' @param x a data frame
#' @param y a data frame
#' @param col_values a character string
#' @param col_names a character string
#'
#' @return a tibble
#'
#'
pairwise_delta <- function(x = NULL,
                           y = NULL,
                           col_values = NULL,
                           col_names = NULL){

  x <- as_tibble(x)
  y <- as_tibble(y)

  x_vals <- x %>%
    dplyr::select(all_of(col_values)) %>%
    unlist()

  y_vals <- y %>%
    dplyr::select(all_of(col_values)) %>%
    unlist()

  x_names <- x %>%
    dplyr::select(all_of(col_names)) %>%
    unlist()

  y_names <- y %>%
    dplyr::select(all_of(col_names)) %>%
    unlist()

  x_vals_n <- x %>% nrow()
  y_vals_n <- y %>% nrow()

  if(x_vals_n != length(unique(x_names))) {
    cli::cli_abort("Error: non-unique names for x")
  }

  if(y_vals_n != length(unique(y_names))) {
    cli::cli_abort("Error: non-unique names for y")
  }

  xy_matrix <- array(x_vals, c(x_vals_n, y_vals_n))

  # matrix form
  xy_dif <- (t(xy_matrix) - y_vals) %>% as.numeric()

  x_row <- array(x_names, c(x_vals_n, y_vals_n)) %>% t() %>% as.character()
  y_row <- array(y_names, c(y_vals_n, x_vals_n)) %>% as.character()

  df <- tibble(
    ref=t(xy_matrix) %>% as.numeric(),
    dif=xy_dif,
    cluster_id=x_row,
    feature_id=y_row)

  return(df)
}
