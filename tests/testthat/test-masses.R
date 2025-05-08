test_that("mass calculations", {

  expect_equal(mass_atomic('Ti'), 47.947947)
  expect_equal(mass_residue('A'), 71.03712)
  expect_equal(mass_charged(1234.5, 5), 247.90728)
  expect_equal(mass_neutral(1234.5, 5), 6167.4636)

  expect_equal(peptide_mass('SAM[P78.0456]LER'), 880.4463)
  expect_equal(peptide_length('SAM[P78.0456]LER'), 7)

  expect_equal(mass_ladder('SAM[P78.0456]LER'),
               c(87.03203,71.03712,131.04048,175.098365,
                 113.08406,129.04259,156.10110))

  expect_equal(mass_fragments('SAM[P78.0456]LER'),
               c(175.11894,88.03931,304.16153,159.07643,417.24559,290.11691,
                 592.34396,465.21527,723.38444,578.29933,794.42156,707.34192))

  expect_equal(peptide_xleucine('SAM[P78.0456]IER'), "SAM[P78.0456]LER")
})
