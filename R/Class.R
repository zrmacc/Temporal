# Purpose: Define class returned by the parametric fitting and contrast functions.
# Updated: 2021-06-05

# -----------------------------------------------------------------------------
# Auxiliary functions.
# -----------------------------------------------------------------------------

#' Round Data Frames
#' 
#' @param df Data.frame.
#' @param digits Integer.
#' @return Data.frame.
#' @importFrom dplyr "%>%"

RoundDF <- function(df, digits = 3) {
  out <- df %>% dplyr::mutate_if(is.numeric, ~ round(., digits = digits))
  return(out)
}


#' Distributions
#' 
#' @param dist Argument passed to FitParaSurv
#' @return String.

DistProperName <- function(dist) {
  lookup <- data.frame(
    arg = c("exp", "gamma", "gen-gamma", "log-normal", "weibull"),
    proper = c("Exponential", "Gamma", "Generalized Gamma", "Log-Normal", "Weibull")
  )
  out <- lookup$proper[lookup$arg == dist]
  return(out)
}


# -----------------------------------------------------------------------------

#' Fitted Survival Distribution
#'
#' Defines the object class returned by fitting functions.
#'
#' @slot Distribution Fitted distribution, string.
#' @slot Parameters Parameters, data.frame.
#' @slot Information Information components, matrix.
#' @slot Outcome Properties of the fitted distribution, data.frame.
#' @slot RMST Estimated restricted mean survival times, data.frame
#' @slot S Fitted survival function, function.
#' @name fit-class
#' @rdname fit-class
#' @exportClass fit

setClass(
  Class = "fit", 
  representation = representation(
    Distribution = "character",
    Parameters = "data.frame",
    Information = "matrix",
    Outcome = "data.frame",
    RMST = "data.frame",
    S = "function"
    )
  )

# -----------------------------------------------------------------------------

#' Print Method for Fitted Survival Distributions
#'
#' Print method for objects of class \code{fit}.
#'
#' @param x An object of class \code{fit}.
#' @param ... Unused.
#' @export

print.fit <- function(x, ...) {
  
  dist <- DistProperName(x@Distribution)
  params <- RoundDF(x@Parameters)
  outcome <- RoundDF(x@Outcome)
  if (length(x@RMST) > 0) {
    rmst <- RoundDF(x@RMST)
  }

  # Display.
  cat(paste0("Fitted ", dist, " Distribution."), "\n")
  cat("Estimated Parameters:\n")
  print(params)
  cat("\n")
  cat("Distribution Properties:\n")
  print(outcome)
  cat("\n")
  
  if (length(x@RMST) > 0) {
    cat("Restricted Mean Survival Times:\n")
    print(rmst)
    cat("\n")
  }
}


# -------------------------------------------------------------------------------------------------------------------------------------------------------

#' Show Method for Fitted Survival Distributions
#'
#' @param object An object of class \code{fit}.
#' @rdname fit-method
#' @importFrom methods show

setMethod(
  f = "show",
  signature = c(object = "fit"),
  definition = function(object) {
    print.fit(x = object)
  }
)


# -----------------------------------------------------------------------------
# Model Contrast
# -----------------------------------------------------------------------------

#' Contrast of Survival Distributions.
#'
#' Defines the object class returned by the comparison function.
#'
#' @slot Dist1 Distribution fit to the target group, string.
#' @slot Dist0 Distribution fit to the reference group, string.
#' @slot Model1 Fitted model for the target group, fit.
#' @slot Model0 Fitted model for the reference group, fit.
#' @slot Location Contrasts of means and medians, data.frame.
#' @slot RMST Contrasts of RMSTs, data.frame.
#' @name contrast-class
#' @rdname contrast-class
#' @exportClass contrast

setClass(
  Class = "contrast",
  representation = representation(
    Dist1 = "character",
    Dist0 = "character",
    Model1 = "fit",
    Model0 = "fit",
    Location = "data.frame",
    RMST = "data.frame"
  )
)

# -----------------------------------------------------------------------------

#' Print Method for a Contrast of Survival Distributions.
#'
#' Print method for an object of class \code{contrast}.
#'
#' @param x A \code{contrast} object.
#' @param ... Unused.
#' @export

print.contrast <- function(x, ...) {
  
  dist1 <- DistProperName(x@Dist1)
  dist0 <- DistProperName(x@Dist0)
  same_dist <- identical(dist1, dist0)
  
  group1 <- RoundDF(x@Model1@Outcome)
  group0 <- RoundDF(x@Model0@Outcome)
  
  loc <- RoundDF(x@Location)
  
  if (length(x@RMST > 0)) {
    rmst <- RoundDF(x@RMST)
  }

  # Display.
  if (same_dist) {
    cat(paste0("Contrast of Fitted ", dist1, " Distributions."), "\n\n")
  } else {
    cat(paste0("Contrast of Fitted ", dist1, " v. ", dist0), "\n\n")
  }
  
  cat("Fitted Characteristics for Group 1:\n")
  print(group1)
  cat("\n")
  
  cat("Fitted Characteristics for Group 0:\n")
  print(group0)
  cat("\n")
  
  cat("Location:\n")
  print(loc)
  cat("\n")
  
  if (length(x@RMST > 0)) {
    cat("RMST:\n")
    print(rmst)
    cat("\n")
  }
}

# -----------------------------------------------------------------------------

#' Show Method for a Contrast of Survival Distributions.
#'
#' @param object An object of class \code{contrast}.
#' @rdname contrast-method
#' @importFrom methods show

setMethod(
  f = "show",
  signature = c(object = "contrast"),
  definition = function(object) {
    print.contrast(x = object)
  }
)

