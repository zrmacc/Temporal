library(Temporal)
library(testthat)

test_that("Check exponential.", {
  
  # Without censoring.
  withr::with_seed(
    seed = 100,
    {
      data <- GenData(n = 1e3, dist = "exp", theta = 1)
      fit <- FitParaSurv(data, dist = "exp") 
    }
  )
  # Check confidence interval for rate contains 1.
  expect_true(all(fit@Parameters$L <= 1, 1 <= fit@Parameters$U))
  
  # With censoring.
  withr::with_seed(
    seed = 101,
    {
      data <- GenData(n = 1e3, dist = "exp", theta = 1, p = 0.2)
      fit <- FitParaSurv(data, dist = "exp") 
    }
  )
  # Check confidence interval for rate contains 1.
  expect_true(all(fit@Parameters$L <= 1, 1 <= fit@Parameters$U))
  
})