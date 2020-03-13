# Purpose: Estimation of Weibull distribution
# Updated: 20/03/07

# Notes
# Numerically verified: Score, Hessian, Mean, Median, Variance

########################
# Weibull Distribution
########################

#' Weibull Distribution Parameter Estimation
#'
#' Estimates parameters for Weibull event times subject to non-informative
#' right censoring. The Weibull distribution is parameterized in terms
#' of the shape \eqn{\alpha} and rate \eqn{\lambda}:
#' \deqn{f(t) = \alpha\lambda^{\alpha}t^{\alpha-1}e^{-(\lambda t)^{\alpha}}, t>0}
#'
#' @param time Numierc observation times.
#' @param status Status indicator, coded as 1 if an event was observed, 0 if censored.
#' @param sig Significance level, for CIs.
#' @param tau Optional truncation times for calculating RMSTs
#' @param init Numeric vector containing the initial value for \eqn{\alpha}. 
#'
#' @importFrom methods new
#' @importFrom stats quantile uniroot
#'
#' @return An object of class \code{fit} containing the following:
#' \describe{
#'  \item{Parameters}{The estimated shape \eqn{\alpha} and rate \eqn{\lambda}.}
#'  \item{Information}{The observed information matrix.}
#'  \item{Outcome}{The fitted mean, median, and variance.}
#'  \item{RMST}{The estimated RMSTs, if tau was specified.}
#' }
#'
#' @seealso
#' \itemize{
#'   \item{Fitting function for parametric survival distributions \code{\link{fitParaSurv}}}
#' }
#'
#' @examples
#' # Generate Weibull data with 20% censoring
#' data = genData(n=1e3,dist="weibull",theta=c(2,2),p=0.2);
#' # Estimate
#' fit = fitParaSurv(time=data$time,status=data$status,dist="weibull");

fit.Weibull = function(time,status,sig=0.05,tau=NULL,init=NULL){

  # Events
  n = length(time);
  nobs = sum(status);
  # Observed events
  tobs = time[status==1];
  # Log times
  logtime = log(time);
  logtobs = log(tobs);
  
  # MLE of rate
  rate = function(a){
    ta = time^a;
    l = (sum(ta)/nobs)^(-1/a);
    return(l);
  }
  
  # Profile score for shape
  score = function(a){
    if(a==0){
      out = Inf;
    }  else {
      # Scaled time
      ta = time^a;
      # Score
      out = nobs/a-nobs*sum(ta*logtime)/sum(ta)+sum(logtobs);
    }
    return(out);
  }
  
  ## Initialize
  a0 = init[1];
  if(!is.numeric(a0)){
    q0 = as.numeric(quantile(x=tobs,probs=c(1-exp(-1))));
    l0 = 1/q0;
    a0 = digamma(1)/(log(l0)+mean(logtobs));
  };

  # Optimize
  ## MLE for alpha
  a = uniroot(f=score,lower=0,upper=2*a0,extendInt="downX")$root;
  ## MLE for lambda
  l = rate(a);
  
  # Observed information
  obsInfo = function(a,l){
    ta = time^a;
    ## Sums
    S1 = sum(ta);
    S2 = sum(ta*logtime);
    S3 = sum(ta*(logtime)^2);
    ## Information components
    Jaa = (nobs/(a^2))+(l^a)*(log(l))^2*S1+2*(l^a)*log(l)*S2+(l^a)*S3;
    Jll = (nobs*a)/(l^2)+a*(a-1)*(l^(a-2))*S1;
    Jla = -(nobs/l)+(l^(a-1))*S1+a*(l^(a-1))*log(l)*S1+a*(l^(a-1))*S2;
    ## Information matrix
    J = matrix(c(Jaa,Jla,Jla,Jll),nrow=2);
    colnames(J) = rownames(J) = c("a","l");
    return(J);
  }
  
  J = obsInfo(a,l);
  Ji = matInv(J);
  
  if(any(diag(Ji)<0)){
    stop("Information matrix not positive definite. Try another initialization");
  }

  # Parameters
  P = data.frame(c("Shape","Rate"),c(a,l),sqrt(diag(Ji)));
  colnames(P) = c("Aspect","Estimate","SE");

  # Fitted Distribution
  ## Estimate mean
  mu = (1/l)*gamma(1+1/a);
  dg = -(1/l)*gamma(1+1/a)*c(digamma(1+1/a)/a^2,1/l);
  dg = matrix(dg,ncol=1);
  se.mu = sqrt(as.numeric(matQF(dg,Ji)));
  ## Numerical check
  # g = function(x){(1/x[2])*gamma(1+1/x[1])};
  # grad(func=g,x=c(a,l));

  ## Estimate median
  me = (1/l)*(log(2))^(1/a);
  dg = -(1/l)*(log(2))^(1/a)*c(log(log(2))/a^2,(1/l));
  dg = matrix(dg,ncol=1);
  se.me = sqrt(as.numeric(matQF(dg,Ji)));
  ## Numerical check
  # g = function(x){(1/x[2])*(log(2))^(1/x[1])};
  # grad(func=g,x=c(a,l));

  ## Estimate variance
  v = (1/l^2)*(gamma(1+2/a)-gamma(1+1/a)^2);
  dg = -(2/l^2)*c((1/a^2)*(gamma(1+2/a)*digamma(1+2/a)-(gamma(1+1/a)^2)*digamma(1+1/a)),
                   (1/l)*(gamma(1+2/a)-gamma(1+1/a)^2));
  dg = matrix(dg,ncol=1);
  se.v = sqrt(as.numeric(matQF(dg,Ji)));
  ## Numerical check
  # g = function(x){(1/x[2]^2)*(gamma(1+2/x[1])-gamma(1+1/x[1])^2)};
  # grad(func=g,x=c(a,l));

  # Outcome characteristics
  Y = data.frame(c("Mean","Median","Variance"),c(mu,me,v),c(se.mu,se.me,se.v));
  colnames(Y) = c("Aspect","Estimate","SE");

  # CIs
  z = qnorm(1-sig/2);
  P$L = P$Estimate-z*P$SE;
  P$U = P$Estimate+z*P$SE;
  Y$L = Y$Estimate-z*Y$SE;
  Y$U = Y$Estimate+z*Y$SE;

  # Fitted survival function
  S = function(t){exp(-(l*t)^a)};

  # Format Results
  Out = new(Class="fit",Distribution="weibull",Parameters=P,Information=J,Outcome=Y,S=S);
  
  ## Add RMST if requested. 
  if(is.numeric(tau)){
    rmst = paraRMST(fit=Out,sig=sig,tau=tau);
    Out@RMST = rmst;
  }
  
  return(Out);
}
