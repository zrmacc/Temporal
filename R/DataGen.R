# Purpose: Censored Data Generation
# Updated: 180828

########################
# Gamma
########################

#' Simulation from the Gamma Distribution
#'
#' Generates gamma event times with shape parameter \eqn{\alpha} and rate
#' parameter \eqn{\lambda}. See \code{\link{fit.Gamma}} for the parameterization. If
#' a censoring proportion \eqn{p} is provided, the event times are subject to
#' non-informative random right censoring.
#'
#' @param n Sample size.
#' @param a Shape.
#' @param l Rate.
#' @param p Expected censoring proportion.
#'
#' @return A data.frame including the observation times and status indicators.
#'
#' @importFrom stats rgamma
#' @importFrom plyr aaply
#' @export
#'
#' @examples
#' # Gamma event times with shape 2 and rate 2
#' # Expected censoring proportion of 20%
#' D = rGamma(n=1e3,a=2,l=2,p=0.2);

rGamma = function(n,a=1,l=1,p=0){
  # Input checks
  if(a<=0){stop("Positive shape parameter is required.")};
  if(l<=0){stop("Positive rate parameter is required.")};
  if(min(p)<0|max(p)>=1){stop("Expected censoring proportion should fall in [0,1).")};
  # Draw gammas
  time = rgamma(n=n,shape=a,rate=l);
  # Return time if no censoring
  if(p==0){
    return(data.frame("time"=time,"status"=rep(1,n)));
  } else {
    # Censoring rate
    q = (1-p)^(1/a);
    b = l*((1-q)/q);
    # Censoring times
    cen = rgamma(n=n,shape=1,rate=b);
    # Observations
    U = cbind(time,cen);
    Umin = aaply(.data=U,.margins=1,.fun=min);
    # Status
    aux = function(x){1*(x[1]<=x[2])};
    d = aaply(.data=U,.margins=1,.fun=aux);
    # Output
    Out = data.frame("time"=Umin,"status"=d);
    return(Out);
  }
}

########################
# Generalized Gamma
########################

#' Simulation from the Generalized Gamma Distribution
#'
#' Generates generalized gamma event times with shape parameters
#' \eqn{(\alpha,\beta)}, and rate parameter \eqn{\lambda}. See
#' \code{\link{fit.GenGamma}} for the parameterization. If a censoring
#' proportion \eqn{p} is provided, the event times are subject to
#' non-informative random right censoring.
#'
#' @param n Sample size.
#' @param a First shape parameter, \eqn{\alpha}.
#' @param b Second shape parameter, \eqn{\beta}. For the standard gamma
#'   distribution, set \eqn{\beta=1}.
#' @param l Rate.
#' @param p Expected censoring proportion.
#'
#' @return A data.frame including the observation times and status indicators.
#'
#' @importFrom plyr aaply
#' @export
#'
#' @examples
#' # Generalized gamma event times with shapes (2,2) and rate 2
#' # Expected censoring proportion of 20%
#' D = rGenGamma(n=1e3,a=2,b=2,l=2,p=0.2);

rGenGamma = function(n,a=1,b=1,l=1,p=0){
  # Input checks
  if((a<=0)|(b<=0)){stop("Positive shape parameters are required.")};
  if(l<=0){stop("Positive rate parameter is required.")};
  if(min(p)<0|max(p)>=1){stop("Expected censoring proportion should fall in [0,1).")};
  # Draw gammas
  time = rgamma(n=n,shape=a,rate=(l^b));
  # Transform to general gammas
  time = (time)^(1/b);

  # Return time if no censoring
  if(p==0){
    return(data.frame("time"=time,"status"=rep(1,n)));
  } else {
    # Censoring rate
    q = (1-p)^(1/a);
    lc = l*((1-q)/q)^(1/b);
    # Censoring times
    cen = rWeibull(n=n,a=b,l=lc)$time;
    # Observations
    U = cbind(time,cen);
    Umin = aaply(.data=U,.margins=1,.fun=min);
    # Status
    aux = function(x){1*(x[1]<=x[2])};
    d = aaply(.data=U,.margins=1,.fun=aux);
    # Output
    Out = data.frame("time"=Umin,"status"=d);
    return(Out);
  }
}

########################
# Log-Logistic
########################

#' Quantile Function for the Log-Logistic Distribution
#'
#' Quantile function for the log-logistic distribution. See \code{\link{fit.LogLogistic}}
#' for the parameterization.
#'
#' @param p Probability.
#' @param a Shape.
#' @param l Rate.
#'
#' @return Scalar quantile.
#'
#' @export
#'
#' @examples
#' # Median of standard log-logistic distribution
#' qLogLogistic(p=0.5);

qLogLogistic = function(p,a=1,l=1){
  # Input checks
  if(a<0){stop("Positive shape parameter is required.")};
  if(l<0){stop("Positive rate parameter is required.")};
  if(min(p)<=0|max(p)>=1){stop("Probability should fall in (0,1).")};
  # Transform
  Out = (1/l)*(p/(1-p))^(1/a);
  return(Out)
}

#' Simulation from the Log-Logistic Distribution
#'
#' Generates log-logistic event times with shape parameter \eqn{\alpha} and rate
#' parameter \eqn{\lambda}. See \code{\link{fit.LogLogistic}} for the
#' parameterization. If a censoring proportion \eqn{p} is provided, the event
#' times are subject to non-informative random right censoring.
#'
#' @param n Sample size.
#' @param a Shape.
#' @param l Rate.
#' @param p Expected censoring proportion.
#'
#' @return A data.frame including the observation times and status indicators.
#'
#' @importFrom stats rnorm qnorm
#' @importFrom plyr aaply
#' @export
#'
#' @examples
#' # Log-logistic event times with shape 4 and rate 1
#' # Expected censoring proportion of 20%
#' D = rLogLogistic(n=1e3,a=4,l=1,p=0.2);

rLogLogistic = function(n,a=4,l=1,p=0){
  # Input checks
  if(a<0){stop("Positive shape parameter is required.")};
  if(l<0){stop("Positive rate parameter is required.")};
  if(min(p)<0|max(p)>=1){stop("Expected censoring proportion should fall in [0,1).")};
  # Draw uniforms
  u = runif(n=n);
  # Event times
  time = qLogLogistic(p=u,a=a,l=l);
  # Return time if no censoring
  if(p==0){
    return(data.frame("time"=time,"status"=rep(1,n)));
  } else {
    # Probability that T arrives before C
    g = function(t){
      if(t==0){
        return(0)
      } else if (t==1){
        return(0.5)
      } else {
        return(t*((t-1)-log(t))/(t-1)^2);
      }
    }
    # If censoring is less than 0.5:
    if(p>0.5){
      q = 1-p;
    } else {
      q = p;
    }
    h = function(x){g(x)-q};
    R = uniroot(f=h,lower=0,upper=1)$root;
    if(p>0.5){
      theta = R;
    } else {
      theta = 1/R;
    }
    # Final censoring rate
    lc = l/(theta^(1/a));
    # Draw uniforms
    v = runif(n=n);
    # Censoring times
    cen = qLogLogistic(p=v,a=a,l=lc);
    # Observations
    U = cbind(time,cen);
    Umin = aaply(.data=U,.margins=1,.fun=min);
    # Status
    aux = function(x){1*(x[1]<=x[2])};
    d = aaply(.data=U,.margins=1,.fun=aux);
    # Output
    Out = data.frame("time"=Umin,"status"=d);
    return(Out);
  }
}

########################
# Log-Normal
########################

#' Simulation from the Log-Normal Distribution
#'
#' Generates log-normal event times with location parameter \eqn{\mu} and scale
#' parameter \eqn{\sigma}. See \code{\link{fit.LogNormal}} for the
#' parameterization. If a censoring proportion \eqn{p} is provided, the event
#' times are subject to non-informative random right censoring.
#'
#' @param n Sample size.
#' @param m Location.
#' @param s Scale.
#' @param p Expected censoring proportion.
#'
#' @return A data.frame including the observation times and status indicators.
#'
#' @importFrom stats rnorm qnorm
#' @importFrom plyr aaply
#' @export
#'
#' @examples
#' # Log-normal event times with location 0 and scale 1
#' # Expected censoring proportion of 20%
#' D = rLogNormal(n=1e3,m=0,s=1,p=0.2);

rLogNormal = function(n,m=0,s=1,p=0){
  # Input checks
  if(s<=0){stop("Positive scale parameter is required.")};
  if(min(p)<0|max(p)>=1){stop("Expected censoring proportion should fall in [0,1).")};
  # Draw gammas
  time = exp(rnorm(n=n,mean=m,sd=s));
  # Return time if no censoring
  if(p==0){
    return(data.frame("time"=time,"status"=rep(1,n)));
  } else {
    # Censoring mean
    mc = m + sqrt(2)*s*qnorm(1-p);
    # Censoring times
    cen = exp(rnorm(n=n,mean=mc,sd=s));
    # Observations
    U = cbind(time,cen);
    Umin = aaply(.data=U,.margins=1,.fun=min);
    # Status
    aux = function(x){1*(x[1]<=x[2])};
    d = aaply(.data=U,.margins=1,.fun=aux);
    # Output
    Out = data.frame("time"=Umin,"status"=d);
    return(Out);
  }
}

########################
# Weibull
########################

#' Quantile Function for the Weibull Distribution
#'
#' Quantile function for the Weibull distribution. See \code{\link{fit.Weibull}}
#' for the parameterization.
#'
#' @param p Probability.
#' @param a Shape.
#' @param l Rate.
#'
#' @return Scalar quantile.

qWeibull = function(p,a=1,l=1){
  # Input checks
  if(a<0){stop("Positive shape parameter is required.")};
  if(l<0){stop("Positive rate parameter is required.")};
  if(min(p)<=0|max(p)>=1){stop("Probability should fall in (0,1).")};
  # Transform
  return((1/l)*(-log(p))^(1/a));
}

#' Simulation from the Weibull Distribution
#'
#' Generates Weibull event times with shape parameter \eqn{\alpha} and rate
#' parameter \eqn{\lambda}. See \code{\link{fit.Weibull}} for the parameterization. If
#' a censoring proportion \eqn{p} is provided, the deviates are subject to
#' non-informative random right censoring.
#'
#' @param n Sample size.
#' @param a Shape.
#' @param l Rate.
#' @param p Expected censoring proportion.
#'
#' @return A data.frame including the observation times and status indicators.
#'
#' @importFrom stats runif
#' @importFrom plyr aaply
#' @export
#'
#' @examples
#' # Weibull event times with shape 2 and rate 2
#' # Expected censoring proportion of 20%
#' D = rWeibull(n=1e3,a=2,l=2,p=0.2);

rWeibull = function(n,a=1,l=1,p=0){
  # Input checks
  if(a<0){stop("Positive shape parameter is required.")};
  if(l<0){stop("Positive rate parameter is required.")};
  if(min(p)<0|max(p)>=1){stop("Expected censoring proportion should fall in [0,1).")};
  # Draw uniforms
  u = runif(n=n);
  # Transform
  time = qWeibull(p=u,a=a,l=l);
  # Return time if no censoring
  if(p==0){
    return(data.frame("time"=time,"status"=rep(1,n)));
  } else {
    # Censoring rate
    lc = (p/(1-p))^(1/a)*l;
    # Censoring times
    v = runif(n=n);
    cen = qWeibull(p=v,a=a,l=lc);
    # Observations
    U = cbind(time,cen);
    Umin = aaply(.data=U,.margins=1,.fun=min);
    # Status
    aux = function(x){1*(x[1]<=x[2])};
    d = aaply(.data=U,.margins=1,.fun=aux);
    # Output
    Out = data.frame("time"=Umin,"status"=d);
    return(Out);
  }
}
