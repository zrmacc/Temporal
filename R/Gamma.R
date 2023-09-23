# Purpose: Parameter estimation for gamma distribution.
# Updated: 2021-07-19

# -----------------------------------------------------------------------------
# Gamma distribution.
# -----------------------------------------------------------------------------

#' Gamma Distribution Parameter Estimation
#'
#' Estimates parameters for gamma event times subject to non-informative
#' right censoring. The gamma distribution is parameterized in terms
#' of the shape \eqn{\alpha} and rate \eqn{\lambda}:
#' \deqn{f(t) = \frac{\lambda}{\Gamma(\alpha)}(\lambda t)^{\alpha-1}e^{-\lambda t}, t>0}
#'
#' @param data Data.frame.
#' @param eps Tolerance for Newton-Raphson iterations.
#' @param init List with initial values for the `shape` \eqn{\alpha} and
#'   `rate` \eqn{\lambda}.
#' @param maxit Maximum number of NR iterations.
#' @param report Report fitting progress?
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
#' # Generate Gamma data with 20% censoring.
#' data <- GenData(n = 1e3, dist = "gamma", theta = c(2, 2), p = 0.2)
#' 
#' # Estimate parameters.
#' fit <- FitParaSurv(data, dist = "gamma")
FitGamma <- function(
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
  
  # Presence of censoring.
  is_cens <- any(data$status == 0)

  # Fitting procedure is simplified if complete data are available. Each branch
  # should return: estimated parameter vector and the information matrix.

  if (!is_cens) {
    
    theta <- FitGammaComplete(data, eps = eps)
    info <- GammaInfo(data, theta[1], theta[2])
  
  } else {

    # Log likelihood.
    ll <- function(theta) {
      out <- SurvLogLik(data, dist = "gamma", theta = theta, log_scale = TRUE)
      return(as.numeric(out))
    }
    
    # Initialize.
    shape0 <- init$shape
    rate0 <- init$rate
    if (!all(is.numeric(rate0), is.numeric(shape0))) {
      theta0 <- FitGammaComplete(data %>% dplyr::filter(status == 1))
    } else {
      theta0 <- c(shape0, rate0)
    }
    theta0 <- log(theta0)
    
    # Estimation.
    theta <- NewtonRaphson(
      init = theta0,
      obj = ll,
      eps = eps,
      maxit = maxit,
      report = report
    )
    theta <- exp(theta)
    
    # Observed information.
    info <- -1 * numDeriv::hessian(
      func = function(t) {
        SurvLogLik(data, dist = "gamma", theta = t)
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
    Aspect = c("Shape", "Rate"), 
    Estimate = c(theta[1], theta[2]), 
    SE = sqrt(diag(inv_info)),
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  # Extract shape and rate.
  shape <- theta[1]
  rate <- theta[2]
 
  # Outcome mean.
  mu <- shape / rate
  grad <- c(1 / rate, -1 * shape / (rate^2))
  se_mu <- sqrt(QF(grad, inv_info))

  # Outcome median.
  g <- function(x) {
    return(stats::qgamma(p = 0.5, shape = x[1], rate = x[2]))
  }
  me <- g(theta)
  grad <- numDeriv::grad(func = g, x = theta)
  se_me <- sqrt(QF(grad, inv_info))

  # Outcome variance.
  v <- shape / (rate^2)
  grad <- c(1 / (rate^2), -2 * shape / (rate^3))
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

  # Fitted survival function.
  surv <- function(t) {return(expint::gammainc(shape, rate * t) / gamma(shape))}

  # Format results.
  out <- methods::new(
    Class = "fit", 
    Distribution = "gamma", 
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

#' Gamma Profile Score for Shape
#' 
#' Profile score equation for gamma event times without censoring.
#' 
#' @param data Data.frame.
#' @param shape Shape parameter.
#' @return Numeric score.
GammaScore <- function(data, shape) {
  
  # Case of invalid shape. 
  if (shape <= 0) {
    return(Inf)
  }
  
  # Formatting.
  n <- nrow(data)
  sum_time <- sum(data$time)
  sum_log_time <- sum(log(data$time))
  
  # Calculation.
  score <- -n * digamma(shape) - n * log(sum_time) + 
      n * log(n * shape) + sum_log_time
  return(score)
}


#' Gamma Observed Information
#' 
#' Observed information for gamme event times without censoring.
#' 
#' @param data Data.frame.
#' @param shape Shape parameter \eqn{\alpha}.
#' @param rate Rate parameter \eqn{\lambda}.
#' @return Numeric information matrix.
GammaInfo <- function(data, shape, rate) {
  n <- nrow(data)
  
  # Information for shape.
  info_shape <- n * trigamma(shape)
  
  # Information for rate.
  info_rate <- n * shape / (rate^2)
  
  # Cross information.
  cross_info <- -(n / rate)
  
  # Information.
  info <- matrix(c(info_shape, cross_info, cross_info, info_rate), nrow = 2)
  dimnames(info) <- list(c("shape", "rate"), c("shape", "rate"))
  
  # Output.
  return(info)
}


#' Gamma Parameter Estimation without Censoring
#'
#' Paramter estimation for gamma event times without censoring.
#'
#' @param data Data.frame.
#' @param eps Tolerance for Newton-Raphson iterations.
#' @return Numeric vector containing the estimated shape and rate parameters.
FitGammaComplete <- function(data, eps = 1e-6) {

  # Moment estimator for initialization.
  shape0 <- mean(data$time)^2 / stats::var(data$time)

  # Numerically zero profile score
  shape <- stats::uniroot(
    f = function(shape) {GammaScore(data, shape)}, 
    lower = 0, 
    upper = 2 * shape0, 
    extendInt = "downX", 
    tol = eps
  )$root
  rate <- nrow(data) * shape / sum(data$time)
  theta <- c(shape = shape, rate = rate)
  return(theta)
}
