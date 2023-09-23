# Purpose: Estimation of the exponential distribution.
# Updated: 2021-07-12

# -----------------------------------------------------------------------------

#' Exponential Distribution Parameter Estimation
#'
#' Estimates parameters for exponential event times subject to non-informative
#' right censoring. The exponential distribution is parameterized in terms
#' of the rate \eqn{\lambda}: \deqn{f(t) = \lambda e^{-\lambda t}, t>0}
#'
#' @param data Data.frame.
#' @param sig Significance level, for CIs.
#' @param status_name Name of the status indicator, 1 if observed, 0 if censored.
#' @param tau Optional truncation times for calculating RMSTs.
#' @param time_name Name of column containing the time to event.
#' @return An object of class \code{fit} containing the following:
#' \describe{
#'  \item{Parameters}{The estimated model parameters.}
#'  \item{Information}{The observed information matrix.}
#'  \item{Outcome}{The fitted mean, median, and variance of the time to event distribution.}
#'  \item{RMST}{The estimated RMSTs, if tau was specified.}
#' }
#' @examples
#' # Generate exponential event time data with 20% censoring.
#' data <- GenData(n = 1e3, dist = "exp", theta = c(2), p = 0.2)
#' 
#' # Estimate parameters.
#' fit <- FitParaSurv(data, dist = "exp")
FitExp <- function(
  data, 
  sig = 0.05, 
  status_name = "status", 
  tau = NULL, 
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
  status <- data$status
  time <- data$time
  
  # Events.
  n <- length(time)
  nobs <- sum(status)
  
  # MLE of lambda.
  rate <- nobs / sum(time)
  
  # Information.
  info <- nobs / (rate^2)
  inv_info <- 1 / info

  # Parameter frame.
  params <- data.frame(
    Aspect = c("Rate"), 
    Estimate = c(rate), 
    SE = sqrt(inv_info),
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  # Outcome mean.
  mu <- 1 / rate
  grad <- -1 / (rate^2)
  se_mu <- sqrt(grad * inv_info * grad)

  # Outcome median.
  me <- (1 / rate) * log(2)
  grad <- -1 / (rate^2) * log(2)
  se_me <- sqrt(grad * inv_info * grad)

  # Outcome variance.
  v <- 1 / (rate^2)
  grad <- -2 / (rate^3)
  se_v <- sqrt(grad * inv_info * grad)

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
  surv <- function(t) {return(exp(-rate * t))}

  # Format results.
  info <- matrix(info)
  dimnames(info) <- list("rate", "rate")
  out <- methods::new(
    Class = "fit", 
    Distribution = "exp", 
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
