# Purpose: Functions related to derivatives of the incomplete gamma
# Updated: 180815

########################
# Derivatives in Alpha
########################

#' First Derivative in Alpha
#'
#' Evaluates the first derivative of the incomplete gamma function
#' with respect to the shape parameter.
#'
#' @param a Shape
#' @param x Lower limit of integration
#' @param h Step size
#' @param method Select "fd" for finite difference, "ni" for numeric
#'   integration.
#' @importFrom expint gammainc
#' @importFrom stats integrate

D1aIncGamma = function(a,x,h=1e-4,method="fd"){
  # Input check
  if(a-2*h<0){stop("Shape-2*h should exceed zero.")};
  if(x<0){stop("Lower limit of integration should be positive.")};
  if(!(method %in% c("fd","ni"))){stop("Select 'fd' for finite difference, 'ni' for numerical integration.")};
  # Finite difference
  if(method=="fd"){
    # Evaluation points
    u = a + c(-2*h,-1*h,h,2*h);
    # Evaluations
    v = gammainc(a=u,x=x);
    # Approximate derivative
    d = as.numeric(v %*% c(1,-8,8,-1))/(12*h);
    return(d)
  }
  if(method=="ni"){
    # Integrand
    g = function(u){u^(a-1)*log(u)*exp(-u)};
    d = integrate(f=g,lower=x,upper=Inf);
    return(d$value);
  }
}

#' Second Derivative in Alpha
#'
#' Evaluates the second derivative of the incomplete gamma function
#' with respect to the shape parameter.
#'
#' @param a Shape
#' @param x Lower limit of integration
#' @param h Step size
#' @param method Select "fd" for finite difference, "ni" for numeric
#'   integration.
#' @importFrom expint gammainc
#' @importFrom stats integrate

D2aIncGamma = function(a,x,h=1e-4,method="ni"){
  # Input check
  if(a-2*h<0){stop("Shape-2*h should exceed zero.")};
  if(x<0){stop("Lower limit of integration should be positive.")};
  if(!(method %in% c("fd","ni"))){stop("Select 'fd' for finite difference, 'ni' for numerical integration.")};
  # Finite difference
  if(method=="fd"){
    # Evaluation points
    u = a + c(-2*h,-1*h,h,2*h);
    d1 = function(a){D1aIncGamma(a=a,x=x,h=h)};
    d1 = Vectorize(d1);
    # Evaluations
    v = d1(u);
    # Approximate derivative
    d2 = as.numeric(v %*% c(1,-8,8,-1))/(12*h);
    return(d2)
  }
  # Numerical integration
  if(method=="ni"){
    # Integrand
    g = function(u){u^(a-1)*(log(u))^2*exp(-u)};
    d2 = integrate(f=g,lower=x,upper=Inf);
    return(d2$value);
  }
}

#' First Logarithmic Derivative in Alpha
#'
#' Evaluates the first derivative of the log incomplete gamma function
#' with respect to the shape parameter.
#'
#' @param a Shape
#' @param x Lower limit of integration
#' @param h Step size
#' @param method Select "fd" for finite difference, "ni" for numeric
#'   integration.
#' @importFrom expint gammainc

D1aLogIncGamma = function(a,x,h=1e-4,method="fd"){
  # Input check
  if(a-2*h<0){stop("Shape-2*h should exceed zero.")};
  if(x<0){stop("Lower limit of integration should be positive.")};
  if(!(method %in% c("fd","ni"))){stop("Select 'fd' for finite difference, 'ni' for numerical integration.")};
  # Numerator
  num = D1aIncGamma(a=a,x=x,h=h,method=method);
  # Denominator
  denom = gammainc(a=a,x=x);
  # Output
  Out = (num/denom);
  return(Out);
}

#' Second Logarithmic Derivative in Alpha
#'
#' Evaluates the second derivative of the log incomplete gamma function
#' with respect to the shape parameter.
#'
#' @param a Shape
#' @param x Lower limit of integration
#' @param h Step size
#' @param method Select "fd" for finite difference, "ni" for numeric
#'   integration.
#' @importFrom expint gammainc

D2aLogIncGamma = function(a,x,h=1e-4,method="fd"){
  # Input check
  if(a-2*h<0){stop("Shape-2*h should exceed zero.")};
  if(x<0){stop("Lower limit of integration should be positive.")};
  if(!(method %in% c("fd","ni"))){stop("Select 'fd' for finite difference, 'ni' for numerical integration.")};
  # Incomplete gamma
  e = gammainc(a=a,x=x);
  # Numerator
  num = D2aIncGamma(a=a,x=x,h=h,method=method)*e-(D1aIncGamma(a=a,x=x,h=h,method=method))^2;
  # Denominator
  denom = e^2;
  # Output
  Out = (num/denom);
  return(Out);
}

########################
# Derivatives in Lambda
########################

#' First Derivative in Lambda
#'
#' Evaluates the first derivative of the incomplete gamma function
#' with respect to the rate parameter.
#'
#' @param a Shape.
#' @param l Rate.
#' @param u Time.

D1lIncGamma = function(a,l,u){
  # Output
  Out = -1*(l^(a-1))*(u^a)*exp(-l*u);
  return(Out);
}

#' Second Derivative in Lambda
#'
#' Evaluates the second derivative of the incomplete gamma function
#' with respect to the rate parameter.
#'
#' @param a Shape.
#' @param l Rate.
#' @param u Time.

D2lIncGamma = function(a,l,u){
  # Output
  Out = -1*(a-1)*(l^(a-2))*(u^a)*exp(-l*u)+(l^(a-1))*(u^(a+1))*exp(-l*u);
  return(Out);
}

#' First Logarithmic Derivative in Lambda
#'
#' Evaluates the first derivative of the log incomplete gamma function
#' with respect to the rate parameter.
#'
#' @param a Shape.
#' @param l Rate.
#' @param u Time.
#' @importFrom expint gammainc

D1lLogIncGamma = function(a,l,u){
  # Numerator
  num = -1*(l^(a-1))*(u^a)*exp(-l*u);
  # Denominor
  denom = gammainc(a,l*u);
  # Output
  Out = (num/denom);
  return(Out);
}

#' Second Logarithmic Derivative in Lambda
#'
#' Evaluates the second derivative of the log incomplete gamma function
#' with respect to the rate parameter.
#'
#' @param a Shape.
#' @param l Rate.
#' @param u Time.
#' @importFrom expint gammainc

D2lLogIncGamma = function(a,l,u){
  # Incomplete gamma
  e = gammainc(a,l*u);
  # Numerator
  num = (D2lIncGamma(a=a,l=l,u=u)*e)-(D1lIncGamma(a=a,l=l,u=u))^2;
  # Denominor
  denom = e^2;
  # Output
  Out = (num/denom);
  return(Out);
}

########################
# Mixed Derivatives
########################

#' Mixed Derivative
#'
#' Evaluates the mixed partial of the incomplete gamma
#' function with respect to the shape and rate parameters.
#'
#' @param a Shape.
#' @param l Rate.
#' @param u Time.

D2mixIncGamma = function(a,l,u){
  # Output
  Out = -log(l*u)*(l^(a-1))*(u^a)*exp(-l*u);
  return(Out);
}

#' Mixed Logarithmic Derivative
#'
#' Evaluates the mixed partial of the log incomplete gamma
#' function with respect to the shape and rate parameters.
#'
#' @param a Shape.
#' @param l Rate.
#' @param u Time.
#' @param h Step size.
#' @importFrom expint gammainc

D2mixLogIncGamma = function(a,l,u,h=1e-4){
  # Incomplete gamma
  e = gammainc(a,l*u);
  # Numerator
  num = (D2mixIncGamma(a=a,l=l,u=u)*e-D1aIncGamma(a=a,x=l*u,h=h)*D1lIncGamma(a=a,l=l,u=u));
  # Denominor
  denom = e^2;
  # Output
  Out = (num/denom);
  return(Out);
}
