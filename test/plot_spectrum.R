#
load_all()
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

spec <- "~/Local/data/project_jester/mzml/acquired/PEL_QC_20230214_QE_Hela200ng_Aur2hr_01.ms2.mzML" |>
  jester::spectra()

conf <-

spec[987,] |>
  xspec(filter = FALSE) |>
  plot_spectrum(peptides = c('SNAEDTLR','SNAE[216.0746]LR','SNAEDT[269.1852]','SN[200.0797]DTLR'),
                charge = 1, tolerance = 0.05)


jester:::psm_search(

)




spec[6,] |>
  xspec() |>
  assign_spectrum('SNAE[216.0746]LR') |>
  dplyr::arrange(mz) |>
  as.data.frame()

spec[6,] |>
  xspec() |>
  assign_spectrum('SN[200.0797]DTLR') |>
  dplyr::arrange(mz) |>
  as.data.frame()
