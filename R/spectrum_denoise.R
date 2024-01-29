#' Clean up a mass spectrum to expose only the most abundant peaks
#'
#' @description
#' `denoise_spectrum()` Removes the noise from a given mass spectrum
#'
#' @param spectrum a spectrum data object
#' @param precursor the precursor mass to remove
#' @param hedge the mz space to remove on either side of a tall peak
#' @param n the number of peaks to retain
#'
#' @export
#'
#'
spectrum_denoise <- function(
    spectrum = NULL,
    precursor = NULL,
    hedge = 5,
    n = 30
){

  if(is.null(spectrum)) { cli::cli_abort("`spectrum` must not be null")}
  if(!is.null(precursor) & !is.numeric(precursor)) { cli::cli_abort("`precursor` must be a numeric")}
  if(!is.null(hedge) & !is.numeric(hedge)) { cli::cli_abort("`hedge` must be a numeric")}
  if(!is.null(n) & !is.numeric(n)) { cli::cli_abort("`n` must be a numeric")}

  if(!is.null(precursor)) {
    mz <- spectrum[,1]
    w_trim <- which(abs(mz - (precursor + 1.3)) > 1.5)
    spectrum <- spectrum[w_trim,]
  }

  if(n >= nrow(spectrum)) { return(spectrum) }

  new_spec <- list()
  for(i in 1:n){
    int <- spectrum[,2]
    mz <- spectrum[,1]

    w_keep <- which(int == max(int))[1]
    # test for monoisotope
    w_mono <- which(abs(mz - (mz[w_keep] - 1.002)) < 0.025)
    if(length(w_mono) > 0){ w_keep <- w_mono[1] }
    w_trim <- which(abs(mz - mz[w_keep]) > hedge)

    if(length(w_trim) == 0) break;

    new_spec[[i]] <- spectrum[w_keep,]
    spectrum <- spectrum[w_trim,]
    if(length(dim(spectrum)) <= 1) break;
  }

  new_spec <- new_spec |> dplyr::bind_rows()

  return(new_spec[order(new_spec[,1]),])

}
