test_that("peptide strings", {

  expect_equal(str_sequence(c("SAM[P78.0456]LER", "SA[M15.99]PLEK")),
               c("SAMPLER","SAMPLEK"))


  expect_equal(str_peptide("SAMP[78.0456]LER"), "SAM[P78.0455]LER")

  expect_equal(peptide_length("SAM[P78.0455]LER"), 7)

  expect_equal(peptide_reverse("SAM[P78.0455]LER"), "EL[P78.0455]MASR")

})
