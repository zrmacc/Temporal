library(Temporal)
library(testthat)

# vcov.fit returns the covariance matrix of the estimated parameters, i.e. the
# inverse of the observed information matrix.

test_that("vcov(fit) returns a matrix with correct dimensions.", {
  withr::with_seed(seed = 301, {
    data <- GenData(n = 200, dist = "exp", theta = 1, p = 0)
    fit <- FitParaSurv(data, dist = "exp")
  })
  v <- vcov(fit)
  expect_true(is.matrix(v))
  expect_equal(nrow(v), 1)
  expect_equal(ncol(v), 1)
  expect_equal(dimnames(v), list("Rate", "Rate"))
})

test_that("vcov(fit) diagonal equals squared SE from Parameters.", {
  withr::with_seed(seed = 302, {
    data <- GenData(n = 300, dist = "gamma", theta = c(2, 2), p = 0.1)
    fit <- FitParaSurv(data, dist = "gamma")
  })
  v <- vcov(fit)
  se_sq <- (fit@Parameters$SE)^2
  expect_equal(as.numeric(diag(v)), as.numeric(se_sq), tolerance = 1e-6)
})

test_that("vcov(fit) is inverse of observed information matrix.", {
  withr::with_seed(seed = 303, {
    data <- GenData(n = 250, dist = "weibull", theta = c(2, 0.5), p = 0.15)
    fit <- FitParaSurv(data, dist = "weibull")
  })
  v <- vcov(fit)
  info <- fit@Information
  id <- diag(nrow(info))
  expect_equal(c(v %*% info), c(id), tolerance = 1e-5)
  expect_equal(c(info %*% v), c(id), tolerance = 1e-5)
})

test_that("vcov(fit) is symmetric.", {
  withr::with_seed(seed = 304, {
    data <- GenData(n = 300, dist = "log-normal", theta = c(1, 1), p = 0.1)
    fit <- FitParaSurv(data, dist = "log-normal")
  })
  v <- vcov(fit)
  expect_equal(v, t(v), tolerance = 1e-10)
})

test_that("vcov(fit) has parameter names for multi-parameter distributions.", {
  withr::with_seed(seed = 305, {
    data <- GenData(n = 200, dist = "gamma", theta = c(2, 2), p = 0)
    fit <- FitParaSurv(data, dist = "gamma")
  })
  v <- vcov(fit)
  expect_equal(rownames(v), c("Shape", "Rate"))
  expect_equal(colnames(v), c("Shape", "Rate"))
})
