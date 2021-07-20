library(Temporal)
library(testthat)

test_that("Check log-likelihood.", {
  
  data <- data.frame(time = 1, status = 1)
  
  # Exponential.
  ll <- SurvLogLik(data = data, dist = "exp", theta = 2)
  truth <- log(2) - 2
  expect_equal(ll, truth)
  
  # Gamma.
  ll <- SurvLogLik(data = data, dist = "gamma", theta = c(2, 2))
  truth <- stats::dgamma(x = 1, shape = 2, rate = 2, log = TRUE)
  expect_equal(ll, truth)
  
  # Generalized gamma.
  ll <- SurvLogLik(data = data, dist = "gen-gamma", theta = c(2, 2, 2))
  truth <- log(2) + log(2) - lgamma(2) + (2 * 2 - 1) * log(2) - (2)^2
  expect_equal(ll, truth)
    
  # Log-normal.
  ll <- SurvLogLik(data = data, dist = "log-normal", theta = c(1, 2))
  truth <- stats::dnorm(x = -0.5, log = TRUE) - log(2)
  expect_equal(ll, truth)
  
  # Weibull.
  ll <- SurvLogLik(data = data, dist = "weibull", theta = c(2, 2))
  truth <- log(2) + 2 * log(2) - (2)^2
  expect_equal(ll, truth)
  
})