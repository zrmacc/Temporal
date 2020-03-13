# Purpose: Log likelihood evaluation function
# Updated: 20/03/01

#' Log Likelihood
#'
#' Evaluates the log-likelihood for a parametric survival distribution.
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
#' @param time Numeric observation times.
#' @param status Status indicator, coded as 1 if an event was observed, 0 if censored.
#' @param dist Distribution, from among: "exp","gamma","gen-gamma","log-normal","weibull".
#' @param theta Parameters, which will vary according to the distribution.
#' @param log.scale Are strictly positive parameters on log-scale?
#'
#' @return Scalar value of the log likelihood.
#' @export
#' 
#' @importFrom expint gammainc
#' @importFrom stats pnorm
#'
#' @examples
#' # Generate Weibull data with 20% censoring.
#' data = genData(n=1e3,dist="weibull",theta=c(2,2),p=0.2);
#' # Log likelihood
#' ll = survLogLik(time=data$time,status=data$status,dist="weibull",theta=c(2,2));
#' 
#' # Generate Gamam data with 10% censoring. 
#' data = genData(n=1e3,dist="gamma",theta=c(2,2),p=0.1);
#' # Log likelihood
#' ll = survLogLik(time=data$time,status=data$status,dist="gamma",theta=c(2,2));

survLogLik = function(time,status=NULL,dist,theta,log.scale=FALSE){
  
  # Input check
  ## Distribution
  pass = checkDist(dist);
  if(!pass){
    stop("Distribution check failed.");
  }
  
  # Events
  n = length(time);
  # Status
  if(is.null(status)){
    status = rep(1,n);
  }
  # Observed events
  nobs = sum(status);
  # Observed event times
  tobs = time[status==1];
  tcen = time[status==0];
  # Presence of censoring
  flag = (length(tcen)>0);

  # Exponential distribution
  if(dist=="exp"){
    # Extract parameters
    l = theta[1];
    if(log.scale){l = exp(l);}
    # Log likelihood
    ll = nobs*log(l)-l*sum(time);
  } else if(dist=="gamma") {
    # Extract parameters
    a = theta[1];
    l = theta[2];
    if(log.scale){
      a = exp(a);
      l = exp(l);
    }
    # Log likelihood
    ll = nobs*a*log(l)+a*sum(log(tobs))-l*sum(tobs)-n*lgamma(a);
    # Add corrections for censoring
    if(flag){ll = ll+sum(log(gammainc(a,l*tcen)));}
  } else if (dist=="gen-gamma"){
    # extract parameters
    a = theta[1];
    b = theta[2];
    l = theta[3];
    if(log.scale){
      a = exp(a);
      b = exp(b);
      l = exp(l);
    }
    # Log likelihood
    ll = nobs*(log(b)+a*b*log(l))+a*b*sum(log(tobs))-(l^b)*sum(tobs^b)-n*lgamma(a);
    # Add corrections for censoring
    if(flag){ll = ll+sum(log(gammainc(a,(l*tcen)^b)));}
  } else if (dist=="log-logistic"){
    # Extract parameters
    a = theta[1];
    l = theta[2];
    if(log.scale){
      a = exp(a);
      l = exp(l);
    }
    # Log likelihood
    ll = nobs*log(a)+nobs*a*log(l)+a*sum(log(tobs))-sum((1+status)*log(1+(l*time)^a));
  } else if (dist=="log-normal"){
    # Extract parameters
    m = theta[1];
    s = theta[2];
    if(log.scale){s = exp(s);}
    # Z scores
    zobs = (log(tobs)-m)/s;
    zcen = (log(tcen)-m)/s;
    # Log likelihood
    ll = -nobs*log(s)-(1/2)*sum(zobs^2)+sum(pnorm(q=zcen,lower.tail=F,log.p=T));
  } else if(dist=="weibull"){
    # Extract parameters
    a = theta[1];
    l = theta[2];
    if(log.scale){
      a = exp(a);
      l = exp(l);
    }
    # Log likelihood
    ll = nobs*log(a)+nobs*a*log(l)+(a-1)*sum(log(tobs))-(l^a)*sum(time^(a));
  }
  # Output
  return(ll);
}
