library(Temporal)
library(testthat)

test_that("GenData produces correct structure for each distribution.", {
  withr::with_seed(seed = 200, {
    for (dist in c("exp", "gamma", "weibull", "log-normal", "gen-gamma")) {
      n <- 100
      theta <- switch(dist,
        "exp" = c(1),
        "gamma" = c(2, 2),
        "weibull" = c(2, 2),
        "log-normal" = c(0, 1),
        "gen-gamma" = c(1, 1, 1)
      )
      data <- GenData(n = n, dist = dist, theta = theta, p = 0.1)
      expect_equal(nrow(data), n)
      expect_true(all(c("time", "status") %in% names(data)))
      expect_true(all(data$time > 0))
      expect_true(all(data$status %in% c(0, 1)))
    }
  })
})

test_that("GenData with p = 0 gives all observed events.", {
  withr::with_seed(seed = 201, {
    data <- GenData(n = 100, dist = "weibull", theta = c(1, 1), p = 0)
    expect_equal(sum(data$status), 100)
  })
})
