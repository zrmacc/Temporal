# Purpose: Parameter estimation for generalized gamma distribution.
# Updated: 2021-07-13

# -----------------------------------------------------------------------------
# Generalized Gamma Distribution
# -----------------------------------------------------------------------------

#' Generalized Gamma Distribution Parameter Estimation
#'
#' Estimates parameters for generalized gamma event times subject to non-informative
#' right censoring. The gamma distribution is parameterized in terms
#' of the shape parameters \eqn{(\alpha,\beta)}, and the rate \eqn{\lambda}:
#' \deqn{f(t) = \frac{\beta\lambda}{\Gamma(\alpha)} (\lambda t)^{\alpha\beta-1}e^{-(\lambda t)^{\beta}}, t>0}
#'
#' @param data Data.frame.
#' @param beta_lower If dist="gen-gamma", lower limit on possible values for beta.
#' @param beta_upper If dist="gen-gamma", upper limit on possible values for beta.
#' @param eps Tolerance for Newton-Raphson iterations.
#' @param init List with initial values for the shape `alpha`, `beta` and rate
#'   `lambda` parameters.
#' @param maxit Maximum number of NR iterations.
#' @param report Report fitting progress?
#' @param sig Significance level, for CIs.
#' @param status_name Name of the status indicator, 1 if observed, 0 if
#'   censored.
#' @param tau Optional truncation times for calculating RMSTs.
#' @param time_name Name of column containing the time to event.
#' @return An object of class \code{fit} containing the following:
#' \describe{
#'  \item{Parameters}{The estimated shape \eqn{(\alpha,\beta)} and rate \eqn{\lambda} parameters.}
#'  \item{Information}{The observed information matrix.}
#'  \item{Outcome}{The fitted mean, median, and variance.}
#'  \item{RMST}{The estimated RMSTs, if tau was specified.}
#' }
#' 
#' @importFrom dplyr "%>%"
#'
#' @examples
#' set.seed(103)
#' # Generate generalized gamma data with 20% censoring.
#' data <- GenData(n = 1e4, dist = "gen-gamma", theta = c(2, 2, 2), p = 0.2)
#' 
#' # Estimate parameters.
#' fit <- FitParaSurv(data, dist = "gen-gamma", report = TRUE)

FitGenGamma <- function(
  data,
  beta_lower = 0.1,
  beta_upper = 10,
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
  # should return: estimated theta, and a function for calculating observed
  # information.

  if (!is_cens) {

    # MLEs of theta
    theta <- FitGenGammaComplete(data, beta_lower, beta_upper)
    info <- GenGammaObsInfo(data, theta[1], theta[2], theta[3])

  } else {

    # Log lkelihood
    ll <- function(theta) {
      out <- SurvLogLik(data, dist = "gen-gamma", theta = theta, log_scale = TRUE)
      return(as.numeric(out))
    }

    # Initialize.
    alpha0 <- init$alpha
    beta0 <- init$beta
    lambda0 <- init$lambda
    if (!all(is.numeric(alpha0), is.numeric(beta0), is.numeric(lambda0))) {
      theta0 <- FitGenGammaComplete(
        data %>% dplyr::filter(status == 1),
        beta_lower = beta_lower,
        beta_upper = beta_upper
      )
    } else {
      theta0 <- c(alpha0, beta0, lambda0)
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
        SurvLogLik(data, dist = "gen-gamma", theta = t)
      },
      x = theta
    )
    
  }
  
  # Inverse information.
  inv_info <- solve(info)

  if (any(diag(inv_info) < 0)) {
    stop("Information matrix not positive definite. Try another initialization")
  }

  # Parameter frame.
  params <- data.frame(
    Aspect = c("Alpha", "Beta", "Lambda"), 
    Estimate = c(theta[1], theta[2], theta[3]), 
    SE = sqrt(diag(inv_info)),
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  # Outcome mean.
  g <- function(theta) {
    alpha <- theta[1]
    beta <- theta[2]
    lambda <- theta[3]
    mu <- gamma(alpha + 1 / beta) / (lambda * gamma(alpha))
    return(as.numeric(mu))
  }
  mu <- g(theta)
  grad <- numDeriv::grad(func = g, x = theta)
  se_mu <- sqrt(QF(grad, inv_info))

  # Outcome median.
  g <- function(theta) {
    alpha <- theta[1]
    beta <- theta[2]
    lambda <- theta[3]
    me <- stats::qgamma(
      p = 0.5, 
      shape = alpha, 
      rate = (lambda^beta)^(1 / beta)
    )
    return(me)
  }
  me <- g(theta)
  grad <- numDeriv::grad(func = g, x = theta)
  se_me <- sqrt(QF(grad, inv_info))

  # Outcome variance.
  g <- function(theta) {
    alpha <- theta[1]
    beta <- theta[2]
    lambda <- theta[3]
    v <- 1 / (lambda^2 * gamma(alpha)) * 
      (gamma(alpha + 2 / beta) - gamma(alpha + 1 / beta)^2 / gamma(alpha))
    return(as.numeric(v))
  }
  v <- g(theta)
  grad <- numDeriv::grad(func = g, x = theta)
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
  alpha <- theta[1]
  beta <- theta[2]
  lambda <- theta[3]
  surv <- function(t) {
    return(expint::gammainc(alpha, (lambda * t)^beta) / gamma(alpha))
  }

  # Format results.
  out <- methods::new(
    Class = "fit", 
    Distribution = "gen-gamma", 
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
  
  return(out)
}

# -----------------------------------------------------------------------------
# Case of no censoring.
# -----------------------------------------------------------------------------

#' Generalized Gamma Rate MLE
#' 
#' Profile MLE of the generalized gamma rate given the shape parameters.
#' 
#' @param data Data.frame.
#' @param alpha First shape parameter.
#' @param beta Second shape parameter.
#' @return Numeric MLE of the rate \eqn{\lambda}.

GenGammaRate <- function(data, alpha, beta) {
  n <- nrow(data)
  time <- data$time
  out <- 1 / (n * alpha) * sum(time^beta)
  out <- 1 / (out^(1 / beta))
  return(out)
}


#' Generalized Gamma Shape MLE
#' 
#' Profile MLE of the first shape parameter \eqn{\alpha} of the generalized gamma
#' given the second shape parameter \eqn{\beta}.
#' 
#' @param data Data.frame.
#' @param beta Second shape parameter.
#' @return Numeric MLE of the rate \eqn{\alpha}.

GenGammaShape <- function(data, beta) {
  n <- nrow(data)
  time <- data$time
  out <- sum((time^beta) * log(time)) / sum(time^beta) - sum(log(time)) / n
  out <- 1 / (beta * out)
  return(out)
}


#' Generalized Gamma Profile Log Likelihood
#' 
#' Profile log likelihood of the generalized gamma distribution as a function
#' of the second shape parameter \eqn{\beta}.
#' 
#' @param data Data.frame.
#' @param beta Second shape parameter.
#' @return Numeric profile log likelihood.

GenGammaProfileLogLik <- function(data, beta) {
  alpha <- GenGammaShape(data, beta)
  lambda <- GenGammaRate(data, alpha, beta)
  out <- SurvLogLik(
    data, 
    theta = c(alpha, beta, lambda), 
    dist = "gen-gamma"
  )
  return(out)
}


#' Generalized Gamma Score Equation
#' 
#' Score equation for the generalized gamma log likelihood in the absence
#' of censoring.
#' 
#' @param data Data.frame.
#' @param alpha First shape parameter.
#' @param beta Second shape parameter.
#' @param lambda Rate parameter.
#' @return Numeric score vector.

GenGammaScore <- function(data, alpha, beta, lambda) {
  n <- nrow(data)
  tobs <- data$time
  
  # Score for alpha.
  score_alpha <- n * beta * log(lambda) + 
    beta * sum(log(tobs))- n * digamma(alpha)
  score_alpha <- as.numeric(score_alpha)
  
  # Score for beta.
  score_beta <- n / beta + 
    n * alpha * log(lambda) + 
    alpha * sum(log(tobs)) - 
    (lambda^beta) * log(lambda) * sum(tobs^beta) - 
    (lambda^beta) * sum((tobs^beta) * log(tobs))
  score_beta <- as.numeric(score_beta)
  
  # Score for lambda.
  score_lambda <- n * alpha * beta / lambda - 
    beta * (lambda^(beta - 1)) * sum(tobs^beta)
  score_lambda <- as.numeric(score_lambda)
  
  # Overall score vector.
  score <- c(
    alpha = score_alpha,
    beta = score_beta,
    lambda = score_lambda
  )
  return(score)
}


#' Generalized Gamma Observed Information
#' 
#' Observed information for the generalized gamma log likelihood in the absence
#' of censoring.
#' 
#' @param data Data.frame.
#' @param alpha First shape parameter.
#' @param beta Second shape parameter.
#' @param lambda Rate parameter.
#' @return Numeric observed information matrix.

GenGammaObsInfo <- function(data, alpha, beta, lambda) {
  n <- nrow(data)
  tobs <- data$time
  
  # Hessian for alpha.
  hess_alpha <- -n * trigamma(alpha)
  hess_alpha <- as.numeric(hess_alpha)
  
  # Hessian for beta.
  hess_beta <- -n / (beta^2) - 
    lambda^beta * log(lambda)^2 * sum(tobs^beta) - 
    2 * lambda^beta * log(lambda) * sum((tobs^beta) * log(tobs)) - 
    lambda^beta * sum((tobs^beta) * (log(tobs)^2))
  hess_beta <- as.numeric(hess_beta)
  
  # Hessian for lambda.
  hess_lambda <- -n * alpha * beta / (lambda^2) - 
    beta * (beta - 1) * lambda^(beta - 2) * sum(tobs^beta)
  hess_lambda <- as.numeric(hess_lambda)
  
  # Mixed partials.
  hess_alpha_beta <- n * log(lambda) + sum(log(tobs))
  hess_alpha_beta <- as.numeric(hess_alpha_beta)
  
  hess_alpha_lambda <- n * beta / lambda
  hess_alpha_lambda <- as.numeric(hess_alpha_lambda)
  
  hess_beta_lambda <- n * alpha / lambda - 
    beta * lambda^(beta - 1) * log(lambda) * sum(tobs^beta) - 
    lambda^(beta - 1) * sum(tobs^beta) - 
    beta * lambda^(beta - 1) * sum((tobs^beta) * log(tobs))
  hess_beta_lambda <- as.numeric(hess_beta_lambda)

  # Observed info.
  info <- -1 * rbind(
    c(hess_alpha, hess_alpha_beta, hess_alpha_lambda),
    c(hess_alpha_beta, hess_beta, hess_beta_lambda),
    c(hess_alpha_lambda, hess_beta_lambda, hess_lambda)
  )
  return(info)
}


#' Generalized Gamma Parameter Estimation without Censoring
#'
#' Paramter estimation for generalized gamma event times without censoring.
#' 
#' @param data Data.frame.
#' @param beta_lower Lower limit on possible values for beta.
#' @param beta_upper Upper limit on possible values for beta.
#' @return Numeric vector containing the estimated shape and rate parameters.

FitGenGammaComplete <- function(data, beta_lower = 0.1, beta_upper = 10) {

  # Profile MLE of beta.
  beta <- stats::optimize(
    f = function(b) {GenGammaProfileLogLik(data, b)}, 
    lower = beta_lower, 
    upper = beta_upper, 
    maximum = TRUE)$maximum

  # Remaining MLEs.
  alpha <- GenGammaShape(data, beta)
  lambda <- GenGammaRate(data, alpha, beta)
  theta <- c(alpha = alpha, beta = beta, lambda = lambda)
  return(theta)
}
