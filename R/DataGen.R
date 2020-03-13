# Purpose: Censored Data Generation
# Updated: 20/03/07

########################
# Master Function
########################

#' Data Generation with Censoring
#' 
#' Generates data from survival distributions as parameterized in this package,
#' with optional non-informative random right censoring. 
#' 
#' The parameter vector theta should contain the following elements, in order,
#' depending on the distribution:
#' \describe{
#'  \item{Exponential}{Rate \eqn{\lambda}.}
#'  \item{Gamma}{Shape \eqn{\alpha}, rate \eqn{\lambda}.}
#'  \item{Generalized Gamma}{Shape 1 \eqn{\alpha}, shape 2 \eqn{\beta}, rate \eqn{\lambda}.}
#'  \item{Log-Normal}{Locaion \eqn{\mu}, scale \eqn{\sigma}.}
#'  \item{Weibull}{Shape \eqn{\alpha}, rate \eqn{\lambda}.}
#' }
#' 
#' @param n Integer sample size.
#' @param dist String, distribution name selected from among:
#'   "exp","gamma","gen-gamma","log-normal","weibull".
#' @param theta Numeric parameter vector. Elements will vary according to the distribution. 
#' @param p Expected censoring proportion.
#'
#' @return data.frame including the observation times and status indicators.
#'
#' @export
#' 
#' @examples
#' # Gamma event times with shape 2 and rate 2
#' # Expected censoring proportion of 20%
#' data = genData(n=1e3,dist="gamma",theta=c(2,2),p=0.20);
#' 
#' # Generalized gamma event times with shapes (2,3) and rate 1
#' # Expected censoring proportion of 15%
#' data = genData(n=1e3,dist="gen-gamma",theta=c(2,3,1),p=0.15);
#' 
#' # Log-normal event times with location 0 and rate 1
#' # Expected censoring proportion of 10%
#' data = genData(n=1e3,dist="log-normal",theta=c(0,1),p=0.10);
#' 
#' # Weibull event times with shape 2 and rate 2
#' # Expected censoring proportion of 5%
#' data = genData(n=1e3,dist="weibull",theta=c(2,2),p=0.05);

genData = function(n,dist="exp",theta=NULL,p=0){
  
  # Defaults
  if(is.null(theta)){
    theta = defaultParam(dist);
  }
  
  # Input check
  pass = checkDist(dist);
  if(!pass){
    stop("Distribution check failed.");
  }
  
  pass = checkTheta(dist,theta);
  if(!pass){
    stop("Theta incorrectly specified.");
  }
  
  # Data generation
  data = NULL;
  if(dist=="exp"){
    data = rWeibull(n=n,a=1,l=theta[1],p=p);
  } else if(dist=="gamma"){
    data = rGamma(n=n,a=theta[1],l=theta[2],p=p);
  } else if(dist=="gen-gamma"){
    data = rGenGamma(n=n,a=theta[1],b=theta[2],l=theta[3],p=p);
  } else if(dist=="log-normal"){
    data = rLogNormal(n=n,m=theta[1],s=theta[2],p=p);
  } else if(dist=="weibull"){
    data = rWeibull(n=n,a=theta[1],l=theta[2],p=p);
  }
  
  # Output
  return(data);
}


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
#' @return Data.frame including the observation times and status indicators.
#'
#' @importFrom stats rgamma
#' @importFrom plyr aaply

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
#' @return DSata.frame including the observation times and status indicators.
#'
#' @importFrom plyr aaply

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
#' @return Data.frame including the observation times and status indicators.
#'
#' @importFrom stats rnorm qnorm
#' @importFrom plyr aaply

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
#' @return Data.frame including the observation times and status indicators.
#'
#' @importFrom stats runif
#' @importFrom plyr aaply

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
