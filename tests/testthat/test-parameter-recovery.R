library(Temporal)
library(testthat)

# Parameter recovery with no censoring (p = 0).
# Point estimates should be close to generative parameters within tolerance.

test_that("Exponential: parameter recovery without censoring.", {
  withr::with_seed(seed = 200, {
    theta_true <- 2
    data <- GenData(n = 3e3, dist = "exp", theta = theta_true, p = 0)
    fit <- FitParaSurv(data, dist = "exp")
  })
  est <- as.numeric(fit@Parameters$Estimate[1])
  expect_equal(est, theta_true, tolerance = 0.08)
})

test_that("Gamma: parameter recovery without censoring.", {
  withr::with_seed(seed = 201, {
    theta_true <- c(2, 2)
    data <- GenData(n = 3e3, dist = "gamma", theta = theta_true, p = 0)
    fit <- FitParaSurv(data, dist = "gamma")
  })
  est <- as.numeric(fit@Parameters$Estimate)
  expect_equal(est[1], theta_true[1], tolerance = 0.1)
  expect_equal(est[2], theta_true[2], tolerance = 0.1)
})

test_that("Weibull: parameter recovery without censoring.", {
  withr::with_seed(seed = 202, {
    theta_true <- c(2, 2)
    data <- GenData(n = 3e3, dist = "weibull", theta = theta_true, p = 0)
    fit <- FitParaSurv(data, dist = "weibull")
  })
  est <- as.numeric(fit@Parameters$Estimate)
  expect_equal(est[1], theta_true[1], tolerance = 0.1)
  expect_equal(est[2], theta_true[2], tolerance = 0.1)
})

test_that("Log-normal: parameter recovery without censoring.", {
  withr::with_seed(seed = 203, {
    theta_true <- c(1, 1)
    data <- GenData(n = 3e3, dist = "log-normal", theta = theta_true, p = 0)
    fit <- FitParaSurv(data, dist = "log-normal")
  })
  est <- as.numeric(fit@Parameters$Estimate)
  expect_equal(est[1], theta_true[1], tolerance = 0.08)
  expect_equal(est[2], theta_true[2], tolerance = 0.08)
})

test_that("Generalized gamma: parameter recovery without censoring.", {
  skip_on_cran()
  withr::with_seed(seed = 204, {
    theta_true <- c(2, 2, 2)
    data <- GenData(n = 5e3, dist = "gen-gamma", theta = theta_true, p = 0)
    fit <- FitParaSurv(data, dist = "gen-gamma")
  })
  est <- as.numeric(fit@Parameters$Estimate)
  expect_equal(est[1], theta_true[1], tolerance = 0.15)
  expect_equal(est[2], theta_true[2], tolerance = 0.15)
  expect_equal(est[3], theta_true[3], tolerance = 0.15)
})
