library(Temporal)
library(testthat)

test_that("Exponential fit recovers rate.", {
  withr::with_seed(
    seed = 100,
    {
      data <- GenData(n = 2e3, dist = "exp", theta = c(2), p = 0.15)
      fit <- FitParaSurv(data, dist = "exp")
    }
  )
  # True rate 2 should lie within the confidence interval.
  expect_true(fit@Parameters$L[1] <= 2 && 2 <= fit@Parameters$U[1])
  expect_equal(fit@Distribution, "exp")
})

test_that("Exponential fit works without censoring.", {
  withr::with_seed(
    seed = 101,
    {
      data <- GenData(n = 500, dist = "exp", theta = c(1))
      fit <- FitParaSurv(data, dist = "exp")
    }
  )
  expect_true(all(fit@Parameters$L <= 1, 1 <= fit@Parameters$U))
})
