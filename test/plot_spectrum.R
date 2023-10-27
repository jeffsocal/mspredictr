#

filter_mono <- function(spectra,
                        precursor = 600,
                        charge = 1,
                        filter = FALSE){
  n <- round(precursor * charge / 114.35) * 3
  spectra <- spectra |> as.data.frame()
  colnames(spectra) <- c('mz', 'intensity')
  if(filter == TRUE){
    spectra <- spectra |> rmstandem::denoise_spectrum(precursor = precursor, n = n)
  }
  return(spectra)
}

xspec <- function(x, filter = FALSE){ filter_mono(x$peaks, x$pre_mz, x$pre_z, filter) }

spec <- "~/Local/data/project_jester/mzml/human_Q15149_180-peptides.mzML" |>
  jester::spectra()


spec[1,] |>
  xspec() |>
  rmstandem::plot_spectrum(peptides = 'QQQQMEQER', charge = 1:2, tolerance = 0.1)
