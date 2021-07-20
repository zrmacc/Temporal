library(Temporal)
library(testthat)

test_that("Check generalized gamma.", {
  skip_on_cran()
  
  # Without censoring.
  withr::with_seed(
    seed = 102,
    {
      data <- GenData(n = 1e3, dist = "gen-gamma", theta = c(1, 1, 1))
      fit <- FitParaSurv(data, dist = "gen-gamma") 
    }
  )
  # Check confidence intervals for all parameters contain 1.
  expect_true(all(fit@Parameters$L <= 1, 1 <= fit@Parameters$U))
  
  # With censoring.
  withr::with_seed(
    seed = 103,
    {
      data <- GenData(n = 1e4, dist = "gen-gamma", theta = c(1, 1, 1), p = 0.2)
      fit <- FitParaSurv(data, dist = "gen-gamma") 
    }
  )
  # Check confidence intervals for all parameters contain 1.
  expect_true(all(fit@Parameters$L <= 1, 1 <= fit@Parameters$U))
  
})