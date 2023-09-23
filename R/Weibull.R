# Purpose: Estimation of Weibull distribution
# Updated: 2021-07-19

# -----------------------------------------------------------------------------
# Weibull distribution.
# -----------------------------------------------------------------------------

#' Weibull Rate MLE
#' 
#' Profile MLE of the Weibull rate as a function of the shape.
#'
#' @param data Data.frame.
#' @param shape Shape parameter.
#' @return Numeric rate.
WeiRate <- function(data, shape) {
  
  # Unpack.
  time <- data$time
  nobs <- sum(data$status)
  
  # Calculation.
  ta <- time^shape
  rate <- exp(-log(sum(ta) / nobs) / shape)
  return(rate)
}


#' Weibull Profile Score for Shape
#' 
#' Profile score equation for the Weibull shape parameter.
#' 
#' @param data Data.frame.
#' @param shape Shape parameter.
#' @return Numeric score.
WeiScore <- function(data, shape) {
  
  # Case of invalid shape. 
  if (shape <= 0) {
    return(Inf)
  }
  
  # Definitions.
  time <- data$time
  logtime <- log(time)
  nobs <- sum(data$status)
  tobs <- data$time[data$status == 1]
  logtobs <- log(tobs)
  
  # Calculation.
  ta <- time^shape
  score <- nobs / shape - nobs * sum(ta * logtime) / sum(ta) + sum(logtobs)
  return(score)
}


#' Weibull Initialization.
#' 
#' @param data Data.frame.
#' @param init Initialization list.
#' @return Numeric initial value for shape.
WeiInit <- function(data, init) {
  a0 <- init$shape
  if (!is.null(a0)) {
    tobs <- data$time[data$status == 1]
    logtobs <- log(tobs)
    q0 <- as.numeric(stats::quantile(x = tobs, probs = c(1 - exp(-1))))
    l0 <- 1 / q0
    a0 <- digamma(1) / (log(l0) + mean(logtobs))
  } 
  a0 <- max(a0, 1e-3)
  return(a0)
}


#' Weibull Information Matrix.
#' 
#' Information matrix for the Weibull shape and rate parameters.
#' 
#' @param data Data.frame.
#' @param shape Shape parameter, alpha.
#' @param rate Rate parameter, lambda.
#' @return Numeric information matrix.
WeiInfo <- function(
  data,
  shape,
  rate
) {
  
  # Unpack.
  time <- data$time
  logtime <- log(time)
  nobs <- sum(data$status)
  tobs <- data$time[data$status == 1]
  logtobs <- log(tobs)
  
  # Calculation.
  ta <- data$time^shape
  s0 <- sum(ta)
  s1 <- sum(ta * logtime)
  s2 <- sum(ta * (logtime)^2)
  
  # Information for shape.
  info_shape <- (nobs / (shape^2)) + 
    (rate^shape) * (log(rate))^2 * s0 + 
    2 * (rate^shape) * log(rate) * s1 + 
    (rate^shape) * s2
  
  # Information for rate.
  info_rate <- (nobs * shape) / (rate^2) + 
    shape * (shape - 1) * (rate^(shape - 2)) * s0
  
  # Cross information.
  cross_info <- -(nobs / rate) + 
    (rate^(shape - 1)) * s0 + 
    shape * (rate^(shape - 1)) * log(rate) * s0 + 
    shape * (rate^(shape - 1)) * s1
  
  # Output.
  info <- matrix(c(info_shape, cross_info, cross_info, info_rate), nrow = 2)
  dimnames(info) <- list(c("shape", "rate"), c("shape", "rate"))
  return(info)
}


#' Weibull Distribution Parameter Estimation
#'
#' Estimates parameters for Weibull event times subject to non-informative
#' right censoring. The Weibull distribution is parameterized in terms
#' of the shape \eqn{\alpha} and rate \eqn{\lambda}:
#' \deqn{f(t) = \alpha\lambda^{\alpha}t^{\alpha-1}e^{-(\lambda t)^{\alpha}}, t>0}
#'
#' @param data Data.frame.
#' @param init List containing the initial value for the shape, \eqn{\alpha}.
#' @param sig Significance level, for CIs.
#' @param status_name Name of the status indicator, 1 if observed, 0 if censored.
#' @param tau Optional truncation times for calculating RMSTs.
#' @param time_name Name of column containing the time to event.
#' @return An object of class \code{fit} containing the following:
#' \describe{
#'  \item{Parameters}{The estimated shape \eqn{\alpha} and rate \eqn{\lambda}.}
#'  \item{Information}{The observed information matrix.}
#'  \item{Outcome}{The fitted mean, median, and variance.}
#'  \item{RMST}{The estimated RMSTs, if tau was specified.}
#' }
#' @examples
#' # Generate Weibull data with 20% censoring.
#' data <- GenData(n = 1e3, dist = "weibull", theta = c(2, 2), p = 0.2)
#' 
#' # Estimate parameters.
#' fit <- FitParaSurv(data, dist = "weibull")
FitWeibull <- function(
  data, 
  init = list(),
  sig = 0.05, 
  status_name = "status", 
  tau = NULL, 
  time_name = "time"
) {
  
  # Data formatting.
  data <- data %>%
    dplyr::rename(
      status = {{ status_name }},
      time = {{ time_name }}
    )

  # Optimize.
  shape0 <- WeiInit(data, init)
  shape <- stats::uniroot(
    f = function(a) {WeiScore(data, a)}, 
    lower = 0, 
    upper = 2 * shape0, 
    extendInt = "downX")$root
  rate <- WeiRate(data, shape)

  # Observed information.
  info <- WeiInfo(data, shape, rate)
  inv_info <- solve(info)
  
  # Check information matrix for positive definiteness.
  if (any(diag(inv_info) < 0)) {
    stop("Information matrix not positive definite. Try another initialization")
  }

  # Parameters.
  params <- data.frame(
    Aspect = c("Shape", "Rate"), 
    Estimate = c(shape, rate), 
    SE = sqrt(diag(inv_info)),
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  # Outcome mean.
  mu <- (1 / rate) * gamma(1 + 1 / shape)
  grad <-  -(1 / rate) * gamma(1 + 1 / shape) * 
    c(digamma(1 + 1 / shape) / shape^2, 1 / rate)
  se_mu <- sqrt(QF(grad, inv_info))
  
  # Numerical verification.
  # g <- function(x) {(1 / x[2]) * gamma(1 + 1 / x[1])}
  # numDeriv::grad(g, x = c(shape, rate))

  # Outcome median.
  me <- (1 / rate) * (log(2))^(1 / shape)
  grad <- -(1 / rate) * (log(2))^(1 / shape) * c(
    log(log(2)) / shape^2, 
    (1 / rate)
  )
  se_me <- sqrt(QF(grad, inv_info))
  
  # Numerical verification.
  # g <- function(x){(1 / x[2]) * (log(2))^(1 / x[1])}
  # numDeriv::grad(g, x = c(shape, rate))

  # Outcome variance.
  v <- (1 / rate^2) * (gamma(1 + 2 / shape) - gamma(1 + 1 / shape)^2)
  grad <- -(2 / rate^2) * c(
    (1 / shape^2) * (gamma(1 + 2 / shape) * digamma(1 + 2 / shape) - 
                       (gamma(1 + 1 / shape)^2) * digamma(1 + 1 / shape)),
    (1 / rate) * (gamma(1 + 2 / shape) - gamma(1 + 1 / shape)^2)
  )
  se_v <- sqrt(QF(grad, inv_info))
  
  # Numerical verification.
  # g <- function(x){(1 / x[2]^2) * (gamma(1 + 2 / x[1]) - gamma(1 + 1 / x[1])^2)}
  # numDeriv::grad(g,x = c(shape, rate))

  # Outcome characteristics.
  outcome <- data.frame(
    Aspect = c("Mean", "Median", "Variance"),
    Estimate = c(mu, me, v), 
    SE = c(se_mu, se_me, se_v),
    stringsAsFactors = FALSE
  )

  # Confidence intervals.
  Estimate <- NULL
  SE <- NULL
  z <- stats::qnorm(1 - sig / 2)
  
  params <- params %>%
    dplyr::mutate(
      L = Estimate - z * SE,
      U = Estimate + z * SE
    )
  
  outcome <- outcome %>%
    dplyr::mutate(
      L = Estimate - z * SE,
      U = Estimate + z * SE
    )

  # Fitted survival function.
  surv <- function(t) {return(exp(-(rate * t)^shape))}

  # Format results.
  out <- methods::new(
    Class = "fit", 
    Distribution = "weibull", 
    Parameters = params, 
    Information = info, 
    Outcome = outcome, 
    S = surv
  )
  
  # RMST.
  if (is.numeric(tau)) {
    rmst <- ParaRMST(fit = out, sig = sig, tau = tau)
    out@RMST <- rmst
  }
  
  # Output.
  return(out)
}
