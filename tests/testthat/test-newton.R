library(Temporal)
library(testthat)

test_that("Check Newton-Raphson maximization.", {
  
  # Normal density.
  f <- function(x){return(stats::dnorm(x, log = TRUE))}
  fit <- NewtonRaphson(init = 1, obj = f)
  expect_equal(fit, 0, tolerance = 1e-6)
  
  # Chi-squared density.
  f <- function(x){return(stats::dchisq(x, df = 3, log = TRUE))}
  fit <- NewtonRaphson(init = 0.2, obj = f)
  expect_equal(fit, 1, tolerance = 1e-6)
  
})