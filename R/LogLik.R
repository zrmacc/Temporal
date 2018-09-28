# Purpose: Log likelihood evaluation function
# Updated: 180926

#' Log Likelihood
#'
#' Evaluation of log-likelihood for parametric survival distribution.
#'
#' @param time Observation times.
#' @param status Status indicator, coded as 1 if an event was observed, 0 if censored.
#' @param dist Distribution.
#' @param theta List of parameters, which will vary according to the distribution.
#' @param log.scale Are positive parameters on log-scale?
#'
#' @return Scalar value of the log likelihood.
#'
#' @importFrom expint gammainc

survLogLik = function(time,status,dist,theta,log.scale=T){
  # Events
  n = length(time);
  nobs = sum(status);
  # Observed event times
  tobs = time[status==1];
  tcen = time[status==0];
  # Presence of censoring
  flag = (length(tcen)>0);

  # Exponential distribution
  if(dist=="exp"){
    # Extract parameters
    l = theta$l;
    # Log likelihood
    ll = nobs*log(l)-l*sum(time);
  }
  # Gamma distribution
  if(dist=="gamma"){
    # Extract parameters
    if(log.scale){
      a = exp(theta$la);
      l = exp(theta$ll);
    } else {
      a = theta$a;
      l = theta$l;
    }
    # Log likelihood
    ll = nobs*a*log(l)+a*sum(log(tobs))-l*sum(tobs)-n*lgamma(a);
    # Add corrections for censoring
    if(flag){ll = ll+sum(log(gammainc(a,l*tcen)));}
  }
  # Generalized gamma
  if(dist=="gengamma"){
    # Extract parameters
    if(log.scale){
      a = exp(theta$la);
      b = exp(theta$lb);
      l = exp(theta$ll);
    } else {
      a = theta$a;
      b = theta$b;
      l = theta$l;
    }
    # Log likelihood
    ll = nobs*(log(b)+a*b*log(l))+a*b*sum(log(tobs))-(l^b)*sum(tobs^b)-n*lgamma(a);
    # Add corrections for censoring
    if(flag){ll = ll+sum(log(gammainc(a,(l*tcen)^b)));}
  }
  # Log logistic
  if(dist=="log-logistic"){
    # Extract parameters
    if(log.scale){
      a = exp(theta$la);
      l = exp(theta$ll);
    } else {
      a = theta$a;
      l = theta$l;
    }
    # Log likelihood
    ll = nobs*log(a)+nobs*a*log(l)+a*sum(log(tobs))-sum((1+status)*log(1+(l*time)^a));
  }
  # Log normal
  if(dist=="log-normal"){
    # Extract parameters
    m = theta$m;
    if(log.scale){
      s = exp(theta$ls);
    } else {
      s = theta$s;
    }
    # Z scores
    zobs = (log(tobs)-m)/s;
    zcen = (log(tcen)-m)/s;
    # Log likelihood
    ll = -nobs*log(s)-(1/2)*sum(zobs^2)+sum(pnorm(q=zcen,lower.tail=F,log.p=T));
  }
  # Weibull distribution
  if(dist=="weibull"){
    # Extract parameters
    a = theta$a;
    l = theta$l;
    # Log likelihood
    ll = nobs*log(a)+nobs*a*log(l)+(a-1)*sum(log(tobs))-(l^a)*sum(time^(a));
  }
  # Output
  return(ll);
}
