# Purpose: Log likelihood evaluation.
# Updated: 2021-07-12

#' Log Likelihood
#'
#' Evaluates the log-likelihood for a parametric survival distribution.
#'
#' The parameter vector theta should contain the following elements, in order,
#' depending on the distribution:
#' \describe{
#'  \item{Exponential}{Rate \eqn{\lambda}.}
#'  \item{Gamma}{Shape \eqn{\alpha}, rate \eqn{\lambda}.}
#'  \item{Generalized Gamma}{Shape 1 \eqn{\alpha}, shape 2 \eqn{\beta}, rate \eqn{\lambda}.}
#'  \item{Log-Normal}{Location \eqn{\mu}, scale \eqn{\sigma}.}
#'  \item{Weibull}{Shape \eqn{\alpha}, rate \eqn{\lambda}.}
#' }
#'
#' @param data Data.frame
#' @param dist Distribution, from among: "exp","gamma","gen-gamma","log-normal","weibull".
#' @param theta Parameters, which will vary according to the distribution.
#' @param log_scale Are strictly positive parameters on log-scale?
#' @param status_name Status indicator, coded as 1 if an event was observed, 0 if censored.
#' @param time_name Name of column containing the time to event.
#' @return Scalar value of the log likelihood.
#' @export
#' @examples
#' # Generate gamma event time data with 10% censoring.
#' data <- GenData(n = 1e3, dist = "gamma", theta = c(2, 2), p = 0.1)
#' 
#' # Evaluate log likelihood.
#' ll <- SurvLogLik(data, dist = "gamma", theta = c(2, 2))
#' 
#' # Generate Weibull event time data with 20% censoring.
#' data <- GenData(n = 1e3, dist = "weibull", theta = c(2, 2), p = 0.2)
#' 
#' # Evaluate log likelihood.
#' ll <- SurvLogLik(data, dist = "weibull", theta = c(2, 2))
SurvLogLik <- function(
  data, 
  dist, 
  theta,
  log_scale = FALSE, 
  status_name = "status",
  time_name = "time"
) {
  
  # Data formatting.
  status <- NULL
  time <- NULL
  data <- data %>%
    dplyr::rename(
      status = {{ status_name }},
      time = {{ time_name }}
    )
  
  time <- data$time
  status <- data$status

  # Input check.
  CheckDist(dist)
  if (!log_scale) {CheckTheta(dist, theta)}

  # Total observations.
  n <- length(data$time)
  
  # Observed events.
  nobs <- sum(data$status)
  
  # Event times.
  tobs <- time[status == 1]
  tcen <- time[status == 0]
  
  # Is censoring present?
  is_cens <- (length(tcen) > 0)

  # Exponential.
  if (dist == "exp") {
    rate <- theta[1]
    if (log_scale) {rate <- exp(rate)}
    ll <- nobs * log(rate) - rate * sum(time)
    return(ll)
  } 
  
  # Gamma.
  if (dist == "gamma") {
    shape <- theta[1]
    rate <- theta[2]
    if (log_scale) {
      shape <- exp(shape)
      rate <- exp(rate)
    }
    ll <- nobs * shape * log(rate) + shape * sum(log(tobs)) - 
      rate * sum(tobs) - n * lgamma(shape)
    if (is_cens) {ll <- ll + sum(log(expint::gammainc(shape, rate * tcen)))}
    return(ll)
  } 
  
  # Generalized gamma.
  if (dist == "gen-gamma") {
    alpha <- theta[1]
    beta <- theta[2]
    rate <- theta[3]
    if (log_scale) {
      alpha <- exp(alpha)
      beta <- exp(beta)
      rate <- exp(rate)
    }
    ll <- nobs * (log(beta) + alpha * beta * log(rate)) + 
      alpha * beta * sum(log(tobs)) - 
      (rate^beta) * sum(tobs^beta) - n * lgamma(alpha)
    if (is_cens) {
      ll <- ll + sum(log(expint::gammainc(alpha, (rate * tcen)^beta)))
    }
    return(ll)
  }
  
  # Log-normal.
  if (dist == "log-normal") {
    loc <- theta[1]
    scale <- theta[2]
    if (log_scale) {scale <- exp(scale)}
    zobs <- (log(tobs) - loc) / scale
    zcen <- (log(tcen) - loc) / scale
    ll <- -nobs * log(scale) - 
      sum(log(tobs)) +
      sum(stats::dnorm(x = zobs, log = TRUE)) + 
      sum(stats::pnorm(q = zcen, lower.tail = FALSE, log.p = TRUE))
    return(ll)
  }
  
  # Weibull.
  if (dist == "weibull") {
    shape <- theta[1]
    rate <- theta[2]
    if (log_scale) {
      shape <- exp(shape)
      rate <- exp(rate)
    }
    ll <- nobs * log(shape) + nobs * shape * log(rate) + 
      (shape - 1) * sum(log(tobs)) - (rate^shape) * sum(time^(shape))
    return(ll)
  }
}
