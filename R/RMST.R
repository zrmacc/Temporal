# Purpose: Calculate area under a fitted parametric survival distribution.
# Updated: 2023-09-23


#' Restricted Mean Survival Time
#'
#' Calculates the RMST as the area under a fitted parametric survival
#' distribution.
#'
#' @param fit Fitted parametric survival distribution.
#' @param tau Numeric vector of truncation times.
#' @param sig Significance level, for CIs.
#' @return Data.frame containing the estimated RMST at each truncation time.
#' @export
#' @examples
#' # Generate Weibull data with 20% censoring.
#' data <- GenData(n = 1e3, dist = "weibull", theta = c(2, 0.5), p = 0.2)
#' 
#' # Fit Weibull distribution.
#' fit <- FitParaSurv(data, dist = "weibull")
#' 
#' # Calculate RMSTs.
#' rmst <- ParaRMST(fit = fit, tau = c(0.5, 1.0, 1.5, 2.0))
#' 
#' # Generate gamma data with 10% censoring.
#' data <- GenData(n = 1e3, dist = "gamma", theta = c(2, 2), p = 0.10)
#' 
#' # Fit gamma distribution.
#' fit <- FitParaSurv(data, dist = "gamma")
#' 
#' # Calculate RMSTs.
#' rmst <- ParaRMST(fit = fit, tau = c(0.5, 1.0, 1.5, 2.0))
ParaRMST <- function(fit, tau, sig = 0.05) {

  # Input check.
  if (!methods::is(fit, "fit")) {
    stop("Requires a fitted parametric survival distribution.")
  }
  if (!is.numeric(tau) | (min(tau) < 0)) {
    stop("Requires positive, numeric truncation times.")
  }
  
  # Critical value.
  z <- stats::qnorm(1 - sig / 2)
  
  # Observed parameters.
  theta_obs <- fit@Parameters$Estimate
  
  # Inverse information
  inv_info <- solve(fit@Information)
  
  # RMST function.
  rmst <- function(theta, upper) {
    surv_func <- SurvFunc(dist = fit@Distribution, theta = theta)
    out <- stats::integrate(f = surv_func, lower = 0, upper = upper)$value
    return(out)
  }

  # Calculate RMST and its SE for each truncation time.
  aux <- function(tt) {
    est <- rmst(theta_obs, tt)
    delta <- numDeriv::grad(
      func = function(theta) {return(rmst(theta, tt))}, 
      x = theta_obs
    )
    se <- sqrt(QF(delta, inv_info))
    out <- data.frame(
      Tau = tt,
      Estimate = est,
      SE = se
    )
    return(out)
  }

  # Loop over truncation times.
  out <- lapply(tau, aux)
  out <- do.call(rbind, out)
  
  # Add confidence interval.
  Estimate <- NULL
  SE <- NULL
  
  out <- out %>%
    dplyr::mutate(
      L = Estimate - z * SE,
      U = Estimate + z * SE
    )
  
  # Output.
  return(out)
}
