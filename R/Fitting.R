# Purpose: Function to fit parameter survival distributions.
# Updated: 2021-07-17

# -----------------------------------------------------------------------------
# Master Fitting Function
# -----------------------------------------------------------------------------

#' Fit Parametric Survival Distribution
#'
#' Estimates parametric survival distributions using event times subject to
#' non-informative right censoring. Available distributions include:
#' exponential, gamma, generalized gamma, log-normal, and Weibull.
#'
#' @param data Data.frame containing the time to event and status.
#' @param beta_lower If dist="gen-gamma", lower limit on possible values for beta.
#' @param beta_upper If dist="gen-gamma", upper limit on possible values for beta.
#' @param dist String, distribution to fit, selected from among: exp, gamma, gen-gamma
#'   log-normal, and weibull.
#' @param eps Tolerance for Newton-Raphson iterations.
#' @param init List of initial parameters. See individual distributions for the
#'   expected parameters. 
#' @param maxit Maximum number of NR iterations.
#' @param report Report fitting progress?
#' @param status_name Name of the status indicator, 1 if observed, 0 if censored.
#' @param sig Significance level, for CIs.
#' @param tau Optional truncation time for calculating RMSTs.
#' @param time_name Name of column containing the time to event.
#' @return An object of class \code{fit} containing the following:
#' \describe{
#'   \item{Parameters}{The estimated shape and rate parameters.}
#'   \item{Information}{The observed information matrix.}
#'   \item{Outcome}{The fitted mean, median, and variance.}
#'  \item{RMST}{The estimated RMSTs, if tau was specified.}
#' }
#' 
#' @importFrom dplyr "%>%"
#' @export
#'
#' @seealso
#' \itemize{
#'   \item{Between group comparison of survival experience \code{\link{CompParaSurv}}}
#'   \item{Exponential distribution \code{\link{FitExp}}}
#'   \item{Gamma distribution \code{\link{FitGamma}}}
#'   \item{Generalized gamma distribution \code{\link{FitGenGamma}}}
#'   \item{Log-normal distribution \code{\link{FitLogNormal}}}
#'   \item{Weibull distribution \code{\link{FitWeibull}}}
#' }
#'
#' @examples
#' # Generate Gamma data with 20% censoring.
#' data <- GenData(n = 1e3, dist = "gamma", theta = c(2, 2), p = 0.2)
#' # Fit gamma distribution.
#' fit <- FitParaSurv(data, dist = "gamma")
#' 
#' # Generate Weibull data with 10% censoring.
#' data <- GenData(n = 1e3, dist = "weibull", theta = c(2, 2), p = 0.1)
#' # Fit weibull distribution, calculate RMST at tau=0.5.
#' fit <- FitParaSurv(data, dist = "weibull", tau = 0.5)

FitParaSurv <- function(
  data,
  beta_lower = 0.1,
  beta_upper = 10,
  dist = "weibull", 
  eps = 1e-6, 
  init = NULL, 
  maxit = 10, 
  report = FALSE,
  sig = 0.05,
  status_name = "status", 
  tau = NULL, 
  time_name = "time"
) {

  # Formatting.
  df <- data %>%
    dplyr::rename(
      time = {{ time_name }},
      status = {{ status_name }}
    )
  n <- length(df$time)
  
  # Positivity.
  if (min(df$time) < 0) {
    stop("Strictly positive observation times are required.")
  }

  # Status.
  CheckStatus(df$status)

  # Distribution selection.
  CheckDist(dist)

  # Initialization.
  CheckInit(dist = dist, init = init)

  # Model fitting.
  if (dist == "exp") {
    
    # Exponential.
    fit <- FitExp(data = df, sig = sig, tau = tau)
    
  } else if (dist == "gamma") {
    
    # Gamma.
    fit <- FitGamma(
      data = df, eps = eps, init = init, maxit = maxit,
      report = report, sig = sig, tau = tau)
    
  } else if (dist == "gen-gamma") {
    
    # Generalized gamma.
    fit <- FitGenGamma(
      data = df, beta_lower = beta_lower, beta_upper = beta_upper,
      eps = eps, init = init, maxit = maxit,
      report = report, sig = sig, tau = tau)
    
  } else if (dist == "log-normal") {
    
    # Log normal.
    fit <- FitLogNormal(
      data = df, eps = eps, init = init, maxit = maxit,
      report = report, sig = sig, tau = tau)
    
  } else if (dist == "weibull") {
    
    # Weibull.
    fit <- FitWeibull(data = df, init = init, sig = sig, tau = tau)
  }
  
  return(fit)
}
