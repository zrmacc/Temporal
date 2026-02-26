library(Temporal)
library(testthat)

test_that("Invalid distribution is rejected.", {
  data <- data.frame(time = c(1, 2), status = c(1, 1))
  expect_error(
    FitParaSurv(data, dist = "invalid"),
    "Select distribution from among"
  )
})

test_that("Invalid arm is rejected.", {
  data <- GenData(n = 50, dist = "exp", theta = c(1), p = 0)
  data$arm <- c(rep(1, 25), rep(2, 25))
  expect_error(
    CompParaSurv(data, dist1 = "exp", dist0 = "exp"),
    "Arm should have 2 levels"
  )
})

test_that("Negative observation times are rejected.", {
  data <- data.frame(time = c(-1, 1), status = c(1, 1))
  expect_error(
    FitParaSurv(data, dist = "exp"),
    "Strictly positive observation times"
  )
})

test_that("FitParaSurv with tau returns RMST in fit object.", {
  withr::with_seed(seed = 110, {
    data <- GenData(n = 200, dist = "gamma", theta = c(2, 2), p = 0.1)
    fit <- FitParaSurv(data, dist = "gamma", tau = 0.5)
  })
  expect_true(nrow(fit@RMST) >= 1)
  expect_true("Tau" %in% names(fit@RMST))
  expect_true("Estimate" %in% names(fit@RMST))
})

test_that("Fit without tau has empty RMST slot and print/show work.", {
  withr::with_seed(seed = 111, {
    data <- GenData(n = 100, dist = "exp", theta = c(1), p = 0)
    fit <- FitParaSurv(data, dist = "exp")
  })
  expect_s3_class(fit@RMST, "data.frame")
  expect_equal(nrow(fit@RMST), 0)
  expect_error(show(fit), NA)
})

test_that("CompParaSurv without tau has empty RMST slot and print/show work.", {
  withr::with_seed(seed = 112, {
    df1 <- GenData(n = 80, dist = "weibull", theta = c(1, 1), p = 0)
    df1$arm <- 1
    df0 <- GenData(n = 80, dist = "weibull", theta = c(1, 1), p = 0)
    df0$arm <- 0
    data <- rbind(df1, df0)
    comp <- CompParaSurv(data, dist1 = "weibull", dist0 = "weibull")
  })
  expect_s3_class(comp@RMST, "data.frame")
  expect_equal(nrow(comp@RMST), 0)
  expect_error(show(comp), NA)
})
