# Purpose: Derivatives of the log normal cdf
# Updated: 180829

#' Logarithmic Derivatives of the Normal Survival Function
#'
#' Evaluates logarithmic derivatives of the normal survival, defined as:
#' \deqn{\ln\Phi(s)=\ln\int_{s}^{\infty}\frac{e^{-(u^2)/2}}{\sqrt{2\pi}}du}
#'
#' @param s Value of the lower limit \eqn{s} at which to evaluate.
#' @param Order Order of the derivative, select from among 1 and 2.
#'
#' @return Scalar value of the derivative.
#'
#' @importFrom stats dnorm pnorm

dLogNormSurv = function(s,Order=1){
  if(!(Order %in% c(1,2))){stop("Select order from among: c(1,2).")};
  # First order
  if(Order==1){
    # Numerator
    num = -1*dnorm(s);
    # Denomintor
    denom = pnorm(s,lower.tail=F);
  } else {
    # Second order
    # Components
    a = dnorm(s);
    b = pnorm(s,lower.tail=F);
    # Numerator
    num = a*(s*b-a);
    # Denomintor
    denom = (b^2);
  }
  # Output
  Out = (num/denom);
  return(Out);
}
