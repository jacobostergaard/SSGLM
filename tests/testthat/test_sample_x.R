# Test particle sampling

context("particle sampling")
library(SSGLM)



test_that("particles are consistently sampled", {

  tmp_fun <- function(x,pars){
    randvec = c(0.1137,0.6223,0.6093,0.6234,0.8609,0.6403,
                0.0095,0.2326,0.6661,0.5143,0.6936,0.5450,
                0.2827,0.9234,0.2923,0.8373,0.2862,0.2668,
                0.1867,0.2322,0.3166,0.3027,0.1590,0.0400,
                0.2188,0.8106,0.5257,0.9147,0.8313,0.0458)

    if(length(randvec)<length(x)){
      print(paste("Input must be max", length(randvec)))
      return(numeric(0))
    } else{
      x_prob = px_x(x,pars);
      idx = which(randvec[1:length(x_prob)] < x_prob)
      x[idx] = -sign(x[idx])
      idx = which(randvec[1:length(x_prob)] > x_prob)
      x[idx] = x[idx]+sign(x[idx])
      return(x)
    }
  }

  x = 50:70
  pars = c(50,0,50,0)

  set.seed(1234)
  expect_equal(tmp_fun(x, pars),as.numeric(sample_x(x, pars)))

  set.seed(1234)
  expect_equal(tmp_fun(-x, pars),as.numeric(sample_x(-x, pars)))


})

