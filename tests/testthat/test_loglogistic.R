# Test function log_logistic function

context("log-logistic function")
library(SSGLM)

tmp_fun <- function(x, a, b){
  x = x/a
  x = x^-b
  x = (1+x)^-1
  return(x)
}

test_that("log-logistic matches correct function", {
  x = seq(-10,10,.01)
  y = tmp_fun(x, 5,10)

  expect_is(log_logistic(x,5,10), "matrix")           # output from c++ is a 1-column matrix
  expect_equal(as.numeric(log_logistic(x,5,10)), y)   # equal tmp-function above
  expect_true(is.na(log_logistic(NA,5,10)))           # expect NAs
  expect_true(is.na(log_logistic(1,NA,10)))           # expect NAs
  expect_true(is.na(log_logistic(1,5,NA)))            # expect NAs
  expect_error(log_logistic(NA))                      # number of inputs must match
  expect_equal(as.numeric(log_logistic(x,5,0)), rep(.5,length(x)))   # probabilites are 0.5 regardless of x, when beta is 0

})

