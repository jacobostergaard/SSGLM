# Test function for probabilities

context("transition and output probabilities")
library(SSGLM)



test_that("transition probabilites match correct log-logistic parameters", {
  x = seq(-10,10,1)
  expect_equal(as.numeric(px_x(x,c(1,0,1,0))),rep(0.5, length(x)))

  x = seq(1,100,1)
  expect_equal(px_x(x,c(1,2,3,4)),log_logistic(x,1,2))

  x = -seq(1,100,1)
  expect_equal(px_x(x,c(1,2,3,4)),log_logistic(x,3,4))

})


test_that("observation probabilites match correct parameters and modes", {
  x = seq(-10,10,1)
  y = 0:4
  expect_equal(py_xy(x,y,1,2),cbind(matrix(dpois(y,1),nr=5,nc=11),matrix(dpois(y,2),nr=5,nc=10)))

})

