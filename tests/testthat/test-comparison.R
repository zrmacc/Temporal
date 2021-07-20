library(Temporal)
library(testthat)

test_that("Check exponential.", {
  skip_on_cran()
  n <- 1000
  
  # Comparison, same medians.
  withr::with_seed(
    seed = 102,
    {
      df1 <- GenData(n = n, dist = "log-normal", theta = c(0, 2), p = 0.2)
      df1$arm <- 1
      df0 <- GenData(n = n, dist = "weibull", theta = c(2, sqrt(log(2))), p = 0.2)
      df0$arm <- 0
      data <- rbind(df1, df0)
      comp <- CompParaSurv(data, dist1 = "log-normal", dist0 = "weibull")
    }
  )
  
  # Expect p-values for difference and ratio of medians are not significant.
  expect_true(all(comp@Location$P[3:4] > 0.05))
  
  
  # Comparison, same means.
  withr::with_seed(
    seed = 103,
    {
      df1 <- GenData(n = n, dist = "gamma", theta = c(2, 2), p = 0.2)
      df1$arm <- 1
      df0 <- GenData(n = n, dist = "weibull", theta = c(2, sqrt(pi) / 2), p = 0.2)
      df0$arm <- 0
      data <- rbind(df1, df0)
      comp <- CompParaSurv(data, dist1 = "gamma", dist0 = "weibull")
    }
  )
  
  # Expect p-values for difference and ratio of means are not significant.
  expect_true(all(comp@Location$P[1:2] > 0.05))
  
  
  # Comparison, same distribution, different censoring levels.
  withr::with_seed(
    seed = 104,
    {
      df1 <- GenData(n = n, dist = "weibull", theta = c(2, 2), p = 0.2)
      df1$arm <- 1
      df0 <- GenData(n = n, dist = "weibull", theta = c(2, 2), p = 0.4)
      df0$arm <- 0
      data <- rbind(df1, df0)
      comp <- CompParaSurv(data, dist1 = "weibull")
    }
  )
  
  # Expect p-values for difference and ratio of means/medians are not significant.
  expect_true(all(comp@Location$P > 0.05))
  
  
})