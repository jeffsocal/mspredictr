#' Convert a peptide string to a named variable
#'
#' @description
#' `spectrum_noise()` Generates a named variable
#'
#' @param df a data frame
#' @param n_peaks total number of peaks, including noise, the spectrum should have
#' @param noise_mean the level at which to inject noise, left NULL the Nth quantile will be used
#' @param noise_quantile the Nth quantile used to assess noise
#' @param noise_sd the fraction of the noise_quantile to use as the noise variability
#'
#' @return a data frame
#'
spectrum_noise <- function(
    df = NULL,
    n_peaks = 300,
    noise_mean = NULL,
    noise_quantile = 0.25,
    noise_sd = .75
){

  if(is.null(df)){ cli::cli_abort("`df` is null!") }
  if(!is.data.frame(df)){ cli::cli_abort("`df` is not a data.frame!") }

  if(!is.numeric(n_peaks)){ cli::cli_abort("`n_peaks` is not a numeric!") }
  if(n_peaks < 1){ cli::cli_abort("`n_peaks` is not a positive whole number.") }

  if(is.null(noise_mean)){
    noise_mean <- 10^quantile(log10(df[,2]), noise_quantile)[1]
  }

  if(!is.numeric(noise_sd)){ cli::cli_abort("`noise_sd` is not a numeric!") }
  if(noise_sd < 0){ cli::cli_abort("`noise_sd` must be greater than 0.") }
  if(noise_sd > 1){ cli::cli_abort("`noise_sd` must be less than 1.") }

  n_df <- floor(n_peaks) - nrow(df)
  if(n_df > 0){
    mass_max <- max(df[,1]) * (1 + 0.0125)
    mass_min <- min(df[,1]) * (1 - 0.0125)
    df_noise <- data.frame(
      mz = sample((mass_min * 100):(mass_max * 100), n_df) / 100,
      i = round(abs(rnorm(n_df, mean = noise_mean, sd = round(noise_mean * noise_sd))))
    )
    colnames(df_noise) <- colnames(df)

    df <- list(df, df_noise) |> dplyr::bind_rows()

    df[,2] <- df[order(df[,1]),2]
    df[,1] <- df[order(df[,1]),1]

  }
  return(df)
}
