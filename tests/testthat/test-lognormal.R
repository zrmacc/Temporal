library(Temporal)
library(testthat)

test_that("Check log normal.", {
  
  # Without censoring.
  withr::with_seed(
    seed = 104,
    {
      data <- GenData(n = 1e3, dist = "log-normal", theta = c(1, 1))
      fit <- FitParaSurv(data, dist = "log-normal") 
    }
  )
  # Check confidence intervals for location and scale contain 1.
  expect_true(all(fit@Parameters$L <= 1, 1 <= fit@Parameters$U))
  
  # With censoring.
  withr::with_seed(
    seed = 105,
    {
      data <- GenData(n = 1e3, dist = "log-normal", theta = c(1, 1), p = 0.2)
      fit <- FitParaSurv(data, dist = "log-normal") 
    }
  )
  # Check confidence intervals for location and scale contain 1.
  expect_true(all(fit@Parameters$L <= 1, 1 <= fit@Parameters$U))
  
})