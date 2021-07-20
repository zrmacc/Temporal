library(Temporal)
library(testthat)

test_that("Check gamma.", {
  
  # Without censoring.
  withr::with_seed(
    seed = 102,
    {
      data <- GenData(n = 1e3, dist = "gamma", theta = c(1, 1))
      fit <- FitParaSurv(data, dist = "gamma") 
    }
  )
  # Check confidence intervals for shape and rate contain 1.
  expect_true(all(fit@Parameters$L <= 1, 1 <= fit@Parameters$U))
  
  # With censoring.
  withr::with_seed(
    seed = 103,
    {
      data <- GenData(n = 1e3, dist = "gamma", theta = c(1, 1), p = 0.2)
      fit <- FitParaSurv(data, dist = "gamma") 
    }
  )
  # Check confidence interval for rate contains 1.
  expect_true(all(fit@Parameters$L <= 1, 1 <= fit@Parameters$U))
  
})