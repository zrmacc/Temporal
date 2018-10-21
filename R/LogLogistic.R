# Purpose: Estimation of log-logistic distribution
# Updated: 180928

# Notes
# Converted parameters to log-scale for optimization
# Numerically verified: Score, Hessian, Mean, Median, Variance

########################
# Log-Logistic Distribution
########################

#' Log-Logistic Distribution Parameter Estimation
#'
#' Estimates parameters for log-logistic event times subject to non-informative
#' right censoring. The log-logistic distribution is parameterized in terms of
#' the shape \eqn{\alpha} and rate \eqn{\lambda}:
#' \deqn{f(t) = \frac{\alpha\lambda(\lambda t)^{\alpha-1}}{[1+(\lambda t)^{\alpha}]^{2}}, t>0}
#'
#' For the log-logistic distribution, the mean is only defined if the shape
#' parameter \eqn{\alpha>1}, and the variance if the shape parameter
#' \eqn{\alpha>2}. Consequently, estimates of the fitted mean and variance are
#' only returned if the estimated shape parameter exceeds these thresholds. For
#' \eqn{\alpha\ll 1}, the fitting function may fail.
#'
#' @param time Observation times.
#' @param status Status indicator, coded as 1 if an event was observed, 0 if censored.
#' @param sig Significance level, for CIs.
#' @param init List of initial parameter values, including the log of the shape
#'   parameter "la" and the log of the rate parameter "ll".
#' @param eps Tolerance for Newton-Raphson iterations.
#' @param maxit Maximum number of NR iterations.
#' @param report Report fitting progress?
#'
#' @importFrom methods new
#' @importFrom stats median
#'
#' @return An object of class \code{fit} containing the following:
#' \describe{
#'  \item{Parameters}{The estimated location \eqn{\mu} and scale \eqn{\sigma}.}
#'  \item{Information}{The observed information matrix.}
#'  \item{Outcome}{The fitted mean, median, and variance.}
#' }
#'
#' @seealso
#' \itemize{
#'   \item{Fitting function for parametric survival distributions \code{\link{fitParaSurv}}}
#' }
#'
#' @examples
#' # Simulate
#' D = rLogLogistic(n=1e3,a=4,l=2);
#' # Estimate
#' M = fitParaSurv(time=D$time,status=D$status,dist="log-logistic");

fit.LogLogistic = function(time,status,sig=0.05,init=NULL,eps=1e-6,maxit=10,report=F){
  # Input check
  if(!is.null(init)){
    if(!all(names(init)%in%c("la","ll"))){stop("For log-logistic, init should contain log shape 'la' and log rate 'll'.")};
  }
  # Events
  n = length(time);
  nobs = sum(status);
  # Observed events
  tobs = time[status==1];
  tcen = time[status==0];
  # Reusable sums
  sum.log.tobs = sum(log(tobs));

  ## Score
  Score = function(theta,log.scale=T){
    # Extract parameters
    if(log.scale){
      a = exp(theta$la);
      l = exp(theta$ll);
    } else {
      a = theta$a;
      l = theta$l;
    }
    ## Calculate score
    # Scaled time
    b = (l*time)^a;
    # For alpha
    Ua = nobs/a+nobs*log(l)+sum.log.tobs-sum((1+status)*b*log(l*time)/(1+b));
    # For lambda
    Ul = nobs*a/l-sum((1+status)*a*(b/l)/(1+b));
    # Score
    Out = c(Ua,Ul);
    # Map to log scale
    if(log.scale){
      Out = Out*c(a,l);
    }
    # Output
    return(Out);
  }

  ## Observed information
  obsInfo = function(theta,log.scale=T){
    # Extract parameters
    if(log.scale){
      a = exp(theta$la);
      l = exp(theta$ll);
      v = c(a,l);
    } else {
      a = theta$a;
      l = theta$l;
    }
    ## Calculate information
    # Scaled time
    b = (l*time)^a;
    # Observed information for alpha
    Jaa = nobs/(a^2)+sum((1+status)*b*(log(l*time))^2/((1+b)^2));
    # Observed information for lambda
    Jll = nobs*a/(l^2)+sum((1+status)*a*(b/(l^2))*((a-1)-b)/((1+b)^2));
    # Cross information
    Jla = -1*nobs/l+sum((1+status)*(b/l)*(1+b+a*log(l*time))/((1+b)^2));
    # Observed info
    Out = matrix(c(Jaa,Jla,Jla,Jll),nrow=2);
    # Map to log scale
    if(log.scale){
      Out = Out*matOP(v,v)-diag(Score(theta=theta,log.scale=T));
    }
    # Output
    return(Out);
  }

  ## NR Update
  Update = function(theta){
    # Initial log likelihood
    q0 = survLogLik(time=time,status=status,dist="log-logistic",theta=theta);
    # Score
    U = Score(theta,log.scale=T);
    # Inverse observed information
    Ji = matInv(obsInfo(theta,log.scale=T));
    # Proposal
    Prop = c(theta$la,theta$ll)+MMP(Ji,U);
    Prop = as.numeric(Prop);
    # Proposed theta
    Out = list("la"=Prop[1],"ll"=Prop[2]);
    # Final log likelihood
    q1 = survLogLik(time=time,status=status,dist="log-logistic",theta=Out);
    # Increment
    Out$d = (q1-q0);
    return(Out);
  }

  ## Initialize
  theta0 = list();
  # Rate
  ll0 = init$ll;
  if(is.null(ll0)){theta0$ll = -log(median(tobs));} else {theta0$ll=ll0};
  # Shape
  la0 = init$la;
  if(is.null(la0)){
    q6 = as.numeric(quantile(tobs,probs=0.6));
    q4 = as.numeric(quantile(tobs,probs=0.4));
    theta0$la = log(log(1.5))-log(2)+log(1/(log(q6)+theta0$ll)-1/(log(q4)+theta0$ll));
  } else {
    theta0$la = la0;
  }

  ## Maximzation
  for(i in 1:maxit){
    # Update
    theta1 = Update(theta0);
    # Accept if increment is positive
    if(theta1$d>0){
      theta0 = theta1;
      if(report){cat("Objective increment: ",signif(theta1$d,digits=3),"\n")}
    }
    # Terminate if increment is below tolerance
    if(theta1$d<eps){
      rm(theta1);
      break;
    }
  };

  ## Fitting report
  if(report){
    if(i<maxit){
      cat(paste0(i-1," update(s) performed before tolerance limit.\n\n"));
    } else {
      cat(paste0(i," update(s) performed without reaching tolerance limit.\n\n"));
    }
  };

  ## Final parameters
  # Log scale
  la1 = theta0$la;
  ll1 = theta0$ll;
  log.point = c(la1,ll1);
  # Standard scale
  a1 = exp(la1);
  l1 = exp(ll1);
  point = c(a1,l1);

  ## Final information
  # Log scale
  J = obsInfo(theta0,log.scale=T);
  Ji = matInv(J);
  # Standard scale
  I = obsInfo(list("a"=a1,"l"=l1),log.scale=F);
  Ii = matInv(I);

  # Standard errors
  log.se = sqrt(diag(Ji));
  se = sqrt(diag(Ii));

  # Parameter frame
  P = data.frame(c("Shape","Rate"),point,se);
  colnames(P) = c("Aspect","Estimate","SE");

  ## Fitted Distribution
  # Estimate mean
  if(a1>1){
    mu = (pi/(a1*l1))/sin(pi/a1);
    dg = -mu*c((sin(pi/a1)-(pi/a1)*cos(pi/a1))/(a1*sin(pi/a1)),(1/l1));
    se.mu = sqrt(as.numeric(matQF(X=dg,A=Ii)));
  } else {
    mu = NA;
    se.mu = NA;
  }

  # Estimate median
  me = (1/l1);
  dg = c(0,-1/(l1^2));
  se.me = sqrt(as.numeric(matQF(X=dg,A=Ii)));

  # Estimate variance
  if(a1>2){
    v = (1/l1)^2*(2*(pi/a1)/sin(2*pi/a1)-(pi/a1)^2/(sin(pi/a1))^2);
    da1 = -2*pi/(a1^2)/(sin(2*pi/a1)^2)*(sin(2*pi/a1)-(2*pi/a1)*cos(2*pi/a1));
    da2 = -2*(pi^2)/(a1^3)/(sin(pi/a1)^4)*(sin(pi/a1)^2-pi/(2*a1)*sin(2*pi/a1));
    dg = c((1/(l1^2))*(da1-da2),-2/(l1)*v);
    se.v = sqrt(as.numeric(matQF(X=dg,A=Ii)));
  } else {
    v = NA;
    se.v = NA;
  }

  # Distribution characteristics
  Y = data.frame(c("Mean","Median","Variance"),c(mu,me,v),c(se.mu,se.me,se.v));
  colnames(Y) = c("Aspect","Estimate","SE");

  # CIs
  z = qnorm(1-sig/2);
  P$L = exp(log.point-z*log.se);
  P$U = exp(log.point+z*log.se);
  Y$L = Y$Estimate-z*Y$SE;
  Y$U = Y$Estimate+z*Y$SE;

  # Fitted survival function
  S = function(t){1/(1+(l1*t)^a1)};

  ## Format Results
  Out = new(Class="fit",Distribution="Log-Logistic",Parameters=P,Information=I,Outcome=Y,S=S);
  return(Out);
}
