library(Temporal)
library(testthat)

test_that("Check arm.", {
  
  # Valid.
  expect_invisible(CheckArm(c(0, 1)))
  
  # Invalid.
  expect_error(CheckArm(c(1, 2)))
  expect_error(CheckArm(c(1, 1)))
})


test_that("Check init.", {
  
  # Valid.
  expect_invisible(CheckInit(dist = "exp", init = NULL))
  expect_invisible(CheckInit(dist = "log-normal", init = list(loc = -1, scale = 1)))
  
  # Initialization not required.
  expect_warning(CheckInit(dist = "exp", init = list(rate = 1)))
  expect_warning(CheckInit(dist = "weibull", init = list(shape = 1, rate = 1)))
  
  # Initialization not supplied as list. 
  expect_error(CheckInit(dist = "gamma", init = c(shape = 1, rate = 1)))
  
  # Initialization invalid.
  expect_error(CheckInit(dist = "gamma", init = c(shape = 0, rate = 1)))
  expect_error(CheckInit(dist = "gen-gamma", init = c(shape = 0, rate = 1)))
  expect_error(CheckInit(dist = "log-normal", init = list(loc = 1, scale = -1)))
})


test_that("Check status.", {
  
  # Valid.
  expect_invisible(CheckStatus(c(0, 1)))
  expect_invisible(CheckStatus(c(1, 1)))
  
  # Invalid.
  expect_error(CheckStatus(c(1, 2)))
  expect_error(CheckStatus(c(0, 0)))
})


test_that("Check theta.", {
  
  # Valid.
  expect_invisible(CheckTheta("exp", c(1)))
  expect_invisible(CheckTheta("gamma", c(1, 1)))
  
  # Invalid.
  expect_error(CheckTheta("exp", c(0)))
  expect_error(CheckTheta("gamma", c(-1, 1)))
  expect_error(CheckTheta("gamma", list(shape = 1, rate = 1)))
})
  