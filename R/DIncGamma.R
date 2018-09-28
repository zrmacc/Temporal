# Purpose: Functions related to derivatives of the incomplete gamma
# Updated: 180828

# Verified numerically.

########################
# Derivatives of Incomplete Gamma
########################

#' Derivatives of the Incomplete Gamma Function
#'
#' Evaluates derivatives of the upper incomplete gamma function, defined as:
#' \deqn{\Gamma(\alpha,s) = \int_{s}^{\infty}u^{\alpha-1}e^{-u}du}
#'
#' @param a Value of the shape \eqn{\alpha} at which to evaluate.
#' @param s Value of the lower limit \eqn{s} at which to evaluate.
#' @param Dirn Direction in which to differentiate. Select from among "a", "s",
#'   and "as".
#' @param Order Order of the derivative, if the direction is either "a", or "s".
#'   Select from among 1 and 2.
#'
#' @return Scalar value of the partial derivative in the direction of interest.
#'
#' @importFrom stats integrate
#' @export
#'
#' @examples
#' # First partial in shape at (a=1,s=1)
#' dIncGamma(a=1,s=1,Dirn="a",Order=1);
#' # Second partial in lower limit at (a=1,s=1)
#' dIncGamma(a=1,s=1,Dirn="s",Order=2);
#' # Mixed partial at (a=1,s=1);
#' dIncGamma(a=1,s=1,Dirn="as");

dIncGamma = function(a,s,Dirn="a",Order=1){
  # Input check
  if(a<=0){stop("a>0 is required.")};
  if(s<=0){stop("s>0 is required.")};
  Dirns = c("a","s","as");
  if(!(Dirn %in% Dirns)){stop("Select dirn from among: c('a','s','as').")};
  if(!(Order %in% c(1,2))){stop("Select order from among: c(1,2).")};

  # Derivatives in a
  if(Dirn=="a"){
    # First order
    if(Order==1){
      # Integrand
      g = function(u){u^(a-1)*log(u)*exp(-u)};
    } else {
      # Integrand
      g = function(u){u^(a-1)*(log(u))^2*exp(-u)};
    }
    # Final derivative
    d = integrate(f=g,lower=s,upper=Inf);
    d = d$value;
  } else if(Dirn=="s"){
    # Derivatives in s
    # First order
    if(Order==1){
      d = -s^(a-1)*exp(-s);
    } else {
      d = s^(a-2)*exp(-s)*(s-(a-1));
    }
    # Final derivative
  } else if(Dirn=="as"){
   # Mixed partial
    d = -s^(a-1)*log(s)*exp(-s);
  }
  return(d);
}

########################
# Derivatives of Log Incomplete Gamma
########################

#' Derivatives of the Log Incomplete Gamma Function
#'
#' Evaluates derivatives of the log of the upper log incomplete gamma function,
#' defined as: \deqn{\ln\Gamma(\alpha,s) = \int_{s}^{\infty}u^{\alpha-1}e^{-u}du}
#'
#' @param a Value of the shape \eqn{\alpha} at which to evaluate.
#' @param s Value of lower limit \eqn{s} at which to evaluate.
#' @param Dirn Direction in which to differentiate. Select from among "a", "s",
#'   and "as".
#' @param Order Order of the derivative, if the direction is either "a", or "s".
#'   Select from among 1 and 2.
#'
#' @return Scalar value of the partial derivative in the direction of interest.
#'
#' @importFrom expint gammainc
#' @importFrom stats integrate
#' @export
#'
#' @examples
#' # First partial in shape at (a=1,s=1)
#' dLogIncGamma(a=1,s=1,Dirn="a",Order=1);
#' # Second partial in lower limit at (a=1,s=1)
#' dLogIncGamma(a=1,s=1,Dirn="s",Order=2);
#' # Mixed partial at (a=1,s=1);
#' dLogIncGamma(a=1,s=1,Dirn="as");

dLogIncGamma = function(a,s,Dirn="a",Order=1){
  # Input check
  if(a<=0){stop("a>0 is required.")};
  if(s<=0){stop("s>0 is required.")};
  Dirns = c("a","s","as");
  if(!(Dirn %in% Dirns)){stop("Select dirn from among: c('a','s','as').")};
  if(!(Order %in% c(1,2))){stop("Select order from among: c(1,2).")};

  # Incomplete gamma evaluation
  G = gammainc(a=a,x=s);

  # Derivatives in a
  if(Dirn=="a"){
    # First order
    if(Order==1){
      d = dIncGamma(a=a,s=s,Dirn="a",Order=1)/G;
    } else {
      # Second order
      a2 = dIncGamma(a=a,s=s,Dirn="a",Order=2);
      a1 = dIncGamma(a=a,s=s,Dirn="a",Order=1);
      d = (a2*G-a1*a1)/(G*G);
    }
  } else if(Dirn=="s"){
    # Derivatives in s
    # First order
    if(Order==1){
      d = dIncGamma(a=a,s=s,Dirn="s",Order=1)/G;
    } else {
      # Second order
      s2 = dIncGamma(a=a,s=s,Dirn="s",Order=2);
      s1 = dIncGamma(a=a,s=s,Dirn="s",Order=1);
      d = (s2*G-s1*s1)/(G*G);
    }
    # Final derivative
  } else {
    # Mixed partial
    m = dIncGamma(a=a,s=s,Dirn="as");
    a1 = dIncGamma(a=a,s=s,Dirn="a",Order=1);
    s1 = dIncGamma(a=a,s=s,Dirn="s",Order=1);
    d = (m*G-a1*s1)/(G*G);
  }
  return(d);
}
