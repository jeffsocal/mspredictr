test_that("small value truncation", {
  expect_equal(num_trunc(0.6789234567, 3), 0.678)
})
