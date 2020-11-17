# Test misc functions

context("Miscellaneous functions used in package")
library(SSGLM)

rev_vec


test_that("Output matches expectations", {
  x = 0:10
  expect_equal(as.numeric(rev_vec(x)),10:0)

  # Test dpois via py_xy
  expect_equal(as.numeric(py_xy(-1,1,1,2)),dpois(1,1))
  expect_equal(as.numeric(py_xy(-1,1:3,1,2)),dpois(1:3,1))

})

