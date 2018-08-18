# Purpose: Derivatives of the log normal cdf
# Updated: 180816

#' First Logarithmic Derivative
#' 
#' @param t Evaluation point
#' @importFrom stats dnorm pnorm

D1LogNormSurv = function(t){
  # Numerator
  num = -1*dnorm(t);
  # Denomintor
  denom = pnorm(t,lower.tail=F);
  # Output
  Out = (num/denom);
  return(Out);
}

#' Second Logarithmic Derivative
#' 
#' @param t Evaluation point
#' @importFrom stats dnorm pnorm

D2LogNormSurv = function(t){
  # Components
  a = dnorm(t);
  b = pnorm(t,lower.tail=F);
  # Numerator
  num = a*(t*b-a);
  # Denomintor
  denom = (b^2);
  # Output
  Out = (num/denom);
  return(Out);
}