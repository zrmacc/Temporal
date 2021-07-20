#' Survival Functions
#'
#' Constructs the survival function for a parameter distribution.
#'
#' The parameter vector theta should contain the following elements, in order,
#' according to the distribution:
#' \describe{
#'  \item{Exponential}{Rate \eqn{\lambda}.}
#'  \item{Gamma}{Shape \eqn{\alpha}, rate \eqn{\lambda}.}
#'  \item{Generalized Gamma}{Shape 1 \eqn{\alpha}, shape 2 \eqn{\beta}, rate \eqn{\lambda}.}
#'  \item{Log-Normal}{Locaion \eqn{\mu}, scale \eqn{\sigma}.}
#'  \item{Weibull}{Shape \eqn{\alpha}, rate \eqn{\lambda}.}
#' }
#'
#' @param dist String, distribution name.
#' @param theta Numeric parameter vector.
#' @return Survival function.
#' @export
#' 
#' @examples 
#' # Survival function for the generalized gamma.
#' surv <- SurvFunc(dist = "gen-gamma", theta = c(2, 2, 2))
#' 
#' # Evaluation.
#' surv(1.0)

SurvFunc <- function(dist, theta) {

  # Input checks.
  CheckDist(dist)
  CheckTheta(dist, theta)

  # Define survival function.
  if (dist == "exp") {
    rate <- theta[1]
    surv <- function(t) {
      return(exp(-rate * t))
    }
  } else if (dist == "gamma") {
    shape <- theta[1]
    rate <- theta[2]
    surv <- function(t) {
      return(expint::gammainc(shape, rate * t) / gamma(shape))
    }
  } else if (dist == "gen-gamma") {
    alpha <- theta[1]
    beta <- theta[2]
    rate <- theta[3]
    surv <- function(t) {
      return(expint::gammainc(alpha, (rate * t)^beta) / gamma(alpha))
    }
  } else if (dist == "log-normal") {
    loc <- theta[1]
    scale <- theta[2]
    surv <- function(t) {
      return(stats::pnorm(q = log(t), mean = loc, sd = scale, lower.tail = FALSE))
    }
  } else if (dist == "weibull") {
    shape <- theta[1]
    rate <- theta[2]
    surv <- function(t) {
      return(exp(-(rate * t)^shape))
    }
  }

  # Output.
  return(surv)
}
