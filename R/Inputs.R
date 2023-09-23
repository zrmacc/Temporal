# Purpose: Functions to check inputs and set defaults.
# Updated: 2021-07-17

# -----------------------------------------------------------------------------

#' Set Default Parameters
#'
#' Function to select default parameter values for each distribution.
#'
#' @param dist String, distribution name.
#' @return Numeric parameter ist.
DefaultParam <- function(dist) {
  if (dist == "exp") {
    theta <- c(rate = 1)
  } else if (dist == "gamma") {
    theta <- c(shape = 1, rate = 1)
  } else if (dist == "gen-gamma") {
    theta <- c(alpha = 1, beta = 1, lambda = 1)
  } else if (dist == "log-normal") {
    theta <- c(loc = 0, scale = 1)
  } else if (dist == "weibull") {
    theta <- c(shape = 1, rate = 1)
  }
  return(theta)
}


# -----------------------------------------------------------------------------

#' Check Arm
#'
#' Check whether treatment arm is properly formatted.
#'
#' @param arm 0/1, treatment arm.
#' @return None.
CheckArm <- function(arm) {
  arm_levels <- sort(unique(arm))
  if (!all.equal(arm_levels, c(0, 1))) {
    stop("Arm should have 2 levels, coded 0 for reference, 1 for treatment.")
  }
  return(invisible(TRUE))
}


# -----------------------------------------------------------------------------

#' Check Distribution
#'
#' Check whether the distribution selected is available.
#'
#' @param dist String, distribution name.
#' @return None.
CheckDist <- function(dist) {
  choices <- c("exp", "gamma", "gen-gamma", "log-normal", "weibull")
  if (!(dist %in% choices)) {
    stop(
      "Select distribution from among the following choices:\n",
      paste(choices, collapse = ", ")
    )
  }
  return(invisible(TRUE))
}


# -----------------------------------------------------------------------------

#' Check Initialization
#'
#' Check whether the initialization is valid.
#'
#' @param dist String, distribution name.
#' @param init List of named parameters.
#' @return None.
CheckInit <- function(dist, init) {
  
  # Only check initialization if not null.
  if (is.null(init)) {return(invisible(TRUE))}
  
  # Exponential.
  if (dist == "exp") {
    warning("Exponential does not require initialization.")
  }
  
  # Check for named list.
  if (!is.list(init)){
    stop("A named list is required for initialization.")
  }
  
  params <- sort(names(init))
  
  # Gamma.
  if (dist == "gamma") {
    if (length(init) != 2) {
      stop("Gamma initialization requires 2 parameters.")
    }
    if (!all.equal(params, c("rate", "shape"))) {
      stop("Gamma initialization requires the following named parameters: rate, shape.")
    }
    if (!all(init > 0)) {
      stop("Gamma shape and rate must be strictly positive.")
    }
  }
  
  # Generalized Gamma.
  if (dist == "gen-gamma") {
    if (length(init) != 3) {
      stop("Generalized gamma initialization requires 2 parameters.")
    }
    if (!all.equal(params, c("alpha", "beta", "lambda"))) {
      stop("Generalized gamma initialization requires the following 
           named parameters: alpha, beta, lambda.")
    }
    if (!all(init > 0)) {
      stop("Generalized gamma parameters must be strictly positive.")
    }
  }
  
  # Log-normal.
  if (dist == "log-normal") {
    if (length(init) != 2) {
      stop("Log-normal initialization requires 2 parameters.")
    }
    if (!all.equal(params, c("loc", "scale"))) {
      stop("Log-normal initialization requires the following 
           named parameters: loc, scale.")
    }
    if (!init$scale > 0) {
      stop("Log-normal scale must be strictly positive.")  
    }
  }
  
  # Weibull.
  if (dist == "weibull") {
    if (length(init) > 1) {
      warning("Only the shape parameter is required for initializing Weibull.")
    }
    if (!(init$shape > 0)) {
      stop("Weibull shape must be strictly positive.")  
    }
  } 
  
  return(invisible(TRUE))
}


# -----------------------------------------------------------------------------

#' Status Check
#'
#' Function to ensure the status indicator is properly formatted
#'
#' @param status 0/1 status indicator.
#' @return None.
CheckStatus <- function(status) {
  status_levels <- sort(unique(status))
  n_status_levels <- length(status_levels)
  
  if (n_status_levels == 1) {
    if (status_levels == 1) {
      return(invisible(TRUE))
    } else {
      stop("If all events are observed, code status as 1 for all records.")
    }
  }
  
  if (n_status_levels == 2) {
    if (all(status_levels == c(0, 1))) {
      return(invisible(TRUE))
    } else {
      stop("Status should have two levels, coded as 0 for censored, 1 for observed.")
    }
  } else {
    stop("Status should only have two levels.")
  }
}

# -----------------------------------------------------------------------------

#' Check Theta
#'
#' Function to check the appropriate number of parameters are supplied for the
#' selected distribution. Used by \code{\link{GenData}}.
#'
#' @param dist String, distribution.
#' @param theta Numeric, parameter vector.
#' @return None.
CheckTheta <- function(dist, theta) {
  len_theta <- length(theta)
  
  # Check numeric.
  if (!is.numeric(theta)) {
    stop("Supply theta as a numeric vector.")
  }
  
  # Exponential.
  if (dist == "exp") {
    if (len_theta != 1) {
      stop("Exponential requires a single rate parameter.")
    }
    if (!(theta[1] > 0)) {
      stop("Exponential rate must be strictly positive.")
    }
  }
  
  # Gamma.
  if (dist == "gamma") {
    if (len_theta != 2) {
      stop("Gamma requires a shape and a rate parameter.")
    }
    if (!all(theta > 0)) {
      stop("Gamma shape and a rate must be strictly positive.")
    }
  }
  
  # Generalized gamma.
  if (dist == "gen-gamma") {
    if (len_theta != 3) {
      stop("Generalized gamma requires two shape parameters and one rate parameter.")
    }
    if (!all(theta > 0)) {
      stop("Gamma shape and a rate must be strictly positive.")
    }
  }
  
  # Log-normal.
  if (dist == "log-normal") {
    if (len_theta != 2) {
      stop("Log-normal requires a location and a scale parameter.")
    }
    if (!(theta[2] > 0)) {
      stop("Log-normal scale paramter must be strictly positive.")
    }
  }
  
  # Weibull.
  if (dist == "weibull") {
    if (len_theta != 2) {
      stop("Weibull requires a shape and a rate parameter.")
    }
    if (!all(theta > 0)) {
      stop("Weibull shape and a rate must be strictly positive.")
    }
  }
  
  return(invisible(TRUE))
}
