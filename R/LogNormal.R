# Purpose: Estimation of log-normal distribution
# Updated: 2021-07-13

# -----------------------------------------------------------------------------
# Log-normal distribution.
# -----------------------------------------------------------------------------

#' Log-Normal Distribution Parameter Estimation
#'
#' Estimates parameters for log-normal event times subject to non-informative
#' right censoring. The log-normal distribution is parameterized in terms
#' of the location \eqn{\mu} and scale \eqn{\sigma}:
#' \deqn{f(t) = \phi\left(\frac{\ln t-\mu}{\sigma}\right)\frac{1}{t\sigma}, t>0}
#'
#' @param data Data.frame.
#' @param eps Tolerance for Newton-Raphson iterations.
#' @param init List with initial values for the location (`loc`) \eqn{\mu} and
#'   `scale` \eqn{\sigma}.
#' @param maxit Maximum number of NR iterations.
#' @param report Report fitting progress?
#' @param sig Significance level, for CIs.
#' @param status_name Name of the status indicator, 1 if observed, 0 if censored.
#' @param tau Optional truncation times for calculating RMSTs.
#' @param time_name Name of column containing the time to event.
#' 
#' @return An object of class \code{fit} containing the following:
#' \describe{
#'  \item{Parameters}{The estimated location \eqn{\mu} and scale \eqn{\sigma}.}
#'  \item{Information}{The observed information matrix.}
#'  \item{Outcome}{The fitted mean, median, and variance.}
#'  \item{RMST}{The estimated RMSTs, if tau was specified.}
#' }
#'
#' @examples
#' # Generate log-normal data with 20% censoring.
#' data <- GenData(n = 1e3, dist = "log-normal", theta = c(0, 2), p = 0.2)
#' 
#' # Estimate parameters.
#' fit <- FitParaSurv(data, dist = "log-normal")

FitLogNormal <- function(
  data,
  eps = 1e-6, 
  init = list(), 
  maxit = 10, 
  report = FALSE,
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
  
  # Fitting procedure is simplified if complete data are available. Each branch
  # should return: estimated parameter vector and the information matrix.
  
  # Presence of censoring.
  is_cens <- any(data$status == 0)
  
  if (!is_cens) {

    theta <- FitLogNormComplete(data)
    info <- LogNormInfo(data, loc = theta["loc"], scale = theta["scale"])
    
  } else {

    # Log likelihood.
    ll <- function(theta) {
      out <- SurvLogLik(data, dist = "log-normal", theta = theta, log_scale = TRUE)
      return(as.numeric(out))
    }

    # Initialize.
    loc0 <- init$loc
    scale0 <- init$scale
    if (!all(is.numeric(loc0), is.numeric(scale0))) {
      theta0 <- FitLogNormComplete(data %>% dplyr::filter(status == 1))
      theta0["scale"] <- log(theta0["scale"])
    } else {
      theta0 <- c(loc0, log(scale0))
    }
    
    # Estimation.
    theta <- NewtonRaphson(
      init = theta0,
      obj = ll,
      eps = eps,
      maxit = maxit,
      report = report
    )
    theta[2] <- exp(theta[2])
    
    # Observed information.
    info <- -1 * numDeriv::hessian(
      func = function(t) {
        SurvLogLik(data, dist = "log-normal", theta = t)
      },
      x = theta
    )
    
  }

  # Inverse information.
  inv_info <- solve(info)
  
  # Check information matrix for positive definiteness.
  if (any(diag(inv_info) < 0)) {
    stop("Information matrix not positive definite. Try another initialization")
  }

  # Parameter frame.
  params <- data.frame(
    Aspect = c("Location", "Scale"), 
    Estimate = theta, 
    SE = sqrt(diag(inv_info)),
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  # Extract location and scale.
  loc <- theta[1]
  scale <- theta[2]
  
  # Outcome mean.
  mu <- exp(loc + scale^2 / 2)
  grad <-mu * c(1, scale)
  se_mu <- sqrt(QF(grad, inv_info))

  # Outcome median.
  me <- exp(loc)
  grad <- c(me, 0)
  se_me <- sqrt(QF(grad, inv_info))

  # Outcome variance.
  v <- (exp(scale^2) - 1) * exp(2 * loc + scale^2)
  grad <- 2 * exp(2 * loc + scale^2) * c(
      exp(scale^2) - 1, 
      scale * (2 * exp(scale^2) - 1)
    )
  se_v <- sqrt(QF(grad, inv_info))

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

  # Fitted survival function
  surv <- function(t) {
    return(stats::pnorm(q = log(t), mean = loc, sd = scale, lower.tail = FALSE))
  }

  # Format results.
  out <- methods::new(
    Class = "fit", 
    Distribution = "log-normal", 
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


# -----------------------------------------------------------------------------
# Case of no censoring.
# -----------------------------------------------------------------------------

#' Log-Normal Score Equation
#' 
#' Score equation for log-normal event times without censoring.
#' 
#' @param data Data.frame.
#' @param loc Location parameter.
#' @param scale Scale parameter.
#' @return Numeric score.

LogNormScore <- function(data, loc, scale) {
  tobs <- data$time
  nobs <- length(tobs)
  zobs <- (log(tobs) - loc) / scale
  score_loc <- (1 / scale) * sum(zobs)
  score_scale <- -1 * (nobs / scale) + (1 / scale) * sum(zobs^2)
  out <- c(score_loc, score_scale)
  return(out)
}


#' Log-Normal Observed Information
#' 
#' Observed information for log-normal event times without censoring.
#' 
#' @param data Data.frame.
#' @param loc Location parameter.
#' @param scale Scale parameter.
#' @param log_scale Is the scale parameter logged?
#' @return Numeric score.

LogNormInfo <- function(data, loc, scale, log_scale = FALSE) {
  n <- nrow(data)
  tobs <- data$time
  zobs <- (log(tobs) - loc) / scale
  
  # Information for location.
  info_loc <- n / (scale^2)
 
  # Information for scale.
  info_scale <- -n / (scale^2) + 3 / (scale^2) * sum(zobs^2)
  
  # Cross info.
  cross_info <- 2 / (scale^2) * sum(zobs)
  
  # Observed info
  info <- matrix(c(info_loc, cross_info, cross_info, info_scale), nrow = 2)
  return(info)
}


#' Log-Normal Parameter Estimation without Censoring
#'
#' @param data Data.frame.
#' @return Numeric vector containing the estimate location and scale parameters.

FitLogNormComplete <- function(data) {
  tobs <- data$time
  theta <- c(
    loc = mean(log(tobs)), 
    scale = stats::sd(log(tobs))
  )
  return(theta)
}
