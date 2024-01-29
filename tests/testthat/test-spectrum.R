test_that("import and match peaks in spectra", {

  expect_no_error(
    suppressMessages(
      mzml <- msreadr::path_to_example() |>
        msreadr::read_spectra()
    )
  )

  expect_no_error(
    suppressMessages(
      match <- mzml |>
        msreadr::subset(spectrum_num == 1) |>
        spectrum_extract() |>
        spectrum_assign('HAVSEGTK')
    )
  )

  expect_equal(round(match$mz, 4),
               c(414.7139,138.0662,209.1033,308.1716,395.2037, 524.2462,
                 262.6271,581.2679,291.1642,682.3152,341.6611,147.1128,
                 248.1605,305.1819,434.2243,521.2567,261.1565,131.0816,
                 620.3247,691.3616,231.1338))

  expect_equal(match$ion,
               c("MH++","b1+","b2+","b3+","b4+","b5+","b5++","b6+","b6++",
                 "b7+","b7++","y1+","y2+","y3+","y4+","y5+","y5++","y5++++",
                 "y6+","y7+","y7+++"))

})
