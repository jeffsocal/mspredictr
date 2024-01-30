#' A helper function to determine the closest value pairs between two vectors
#' in one dimensional space.
#'
#' @description
#' `pairwise_delta()` determine the closest values
#'
#' @param x
#' A data frame of values.
#'
#' @param y
#' A data frame of values.
#'
#' @param col_values
#' A numeric column of values in both x and y to cluster.
#'
#' @param col_names
#' A character string column of unique values in both x and y to cluster.
#'
pairwise_delta <- function(x = NULL,
                           y = NULL,
                           col_values = NULL,
                           col_names = NULL){

  x <- tibble::as_tibble(x)
  y <- tibble::as_tibble(y)

  x_vals <- x %>%
    dplyr::select(dplyr::all_of(col_values)) %>%
    unlist()

  y_vals <- y %>%
    dplyr::select(dplyr::all_of(col_values)) %>%
    unlist()

  x_names <- x %>%
    dplyr::select(dplyr::all_of(col_names)) %>%
    unlist()

  y_names <- y %>%
    dplyr::select(dplyr::all_of(col_names)) %>%
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

  df <- tibble::tibble(
    ref=t(xy_matrix) %>% as.numeric(),
    dif=xy_dif,
    cluster_id=x_row,
    feature_id=y_row)

  return(df)
}
