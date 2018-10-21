# Purpose: Estimation of log-normal distribution
# Updated: 180927

# Notes
# Numerically verified: Score, Hessian, Mean, Median, Variance

########################
# Log-Normal Distribution
########################

#' Log-Normal Distribution Parameter Estimation
#'
#' Estimates parameters for log-normal event times subject to non-informative
#' right censoring. The log-normal distribution is parameterized in terms
#' of the location \eqn{\mu} and scale \eqn{\sigma}:
#' \deqn{f(t) = \phi\left(\frac{\ln t-\mu}{\sigma}\right)\frac{1}{t\sigma}, t>0}
#'
#' @param time Observation times.
#' @param status Status indicator, coded as 1 if an event was observed, 0 if
#'   censored.
#' @param sig Significance level, for CIs.
#' @param init List of initial parameter values, including the location "m", and
#'   the log of the scale parameter "ls".
#' @param eps Tolerance for Newton-Raphson iterations.
#' @param maxit Maximum number of NR iterations.
#' @param report Report fitting progress?
#'
#' @importFrom methods new
#' @importFrom stats var
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
#' D = rLogNormal(n=1e3,m=0,s=1);
#' # Estimate
#' M = fitParaSurv(time=D$time,status=D$status,dist="log-normal");

fit.LogNormal = function(time,status,sig=0.05,init=NULL,eps=1e-6,maxit=10,report=F){
  # Input check
  if(!is.null(init)){
    if(!all(names(init)%in%c("m","ls"))){stop("For log-normal, init should contain location 'm' and log scale 'ls'.")};
  }
  # Events
  n = length(time);
  nobs = sum(status);
  # Observed events
  tobs = time[status==1];
  tcen = time[status==0];

  ## Score
  Score = function(theta,log.scale=T){
    # Extract parameters
    m = theta$m;
    if(log.scale){
      s = exp(theta$ls);
    } else {
      s = theta$s;
    }
    ## Calculate score
    # Standardize
    zobs = (log(tobs)-m)/s;
    zcen = (log(tcen)-m)/s;
    # Score for mu
    Um = (1/s)*sum(zobs)-(1/s)*sum(dLogNormSurv(zcen));
    # Score for sigma
    Us = -1*(nobs/s)+(1/s)*sum(zobs^2)-(1/s)*sum(dLogNormSurv(zcen)*zcen);
    # Score
    Out = c(Um,Us);
    # Map to log scale
    if(log.scale){
      Out = Out*c(1,s);
    }
    # Output
    return(Out);
  }

  ## Function to calculate observed information
  obsInfo = function(theta,log.scale=T){
    # Extract parameters
    m = theta$m;
    if(log.scale){
      s = exp(theta$ls);
    } else {
      s = theta$s;
    }
    ## Calculate information
    # Standardize
    zobs = (log(tobs)-m)/s;
    zcen = (log(tcen)-m)/s;
    # Information for mu
    Jmm = nobs/(s^2)-1/(s^2)*sum(dLogNormSurv(zcen,Order=2));
    # Information for sigma
    Jss = -nobs/(s^2)+3/(s^2)*sum(zobs^2)-2/(s^2)*sum(dLogNormSurv(zcen)*zcen)-1/(s^2)*sum(dLogNormSurv(zcen,Order=2)*(zcen)^2);
    # Cross information
    Jms = 2/(s^2)*sum(zobs)-1/(s^2)*sum(dLogNormSurv(zcen))-1/(s^2)*sum(dLogNormSurv(zcen,Order=2)*zcen);
    # Observed info
    Out = matrix(c(Jmm,Jms,Jms,Jss),nrow=2);
    # Map to log scale
    if(log.scale){
      Out = Out*matOP(c(1,s),c(1,s))-diag(c(0,Score(theta,log.scale=T)[2]));
    }
    # Output
    return(Out);
  }

  ## NR Update
  Update = function(theta){
    # Initial log likelihood
    q0 = survLogLik(time=time,status=status,dist="log-normal",theta=theta);
    # Score
    U = Score(theta,log.scale=T);
    # Inverse observed information
    Ji = matInv(obsInfo(theta,log.scale=T));
    # Proposal
    Prop = c(theta$m,theta$ls)+MMP(Ji,U);
    Prop = as.numeric(Prop);
    # Proposed theta
    Out = list("m"=Prop[1],"ls"=Prop[2]);
    # Final log likelihood
    q1 = survLogLik(time=time,status=status,dist="log-normal",theta=Out);
    # Increment
    Out$d = (q1-q0);
    return(Out);
  }

  ## Initialize
  theta0 = list();
  # Location
  m0 = init$m;
  if(is.null(m0)){theta0$m = mean(log(tobs))} else {theta0$m = m0};
  # Scale
  ls0 = init$ls;
  if(is.null(ls0)){theta0$ls = (1/2)*log(var(log(tobs)))} else {theta0$ls = ls0};

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

  ## Final estimates
  m1 = theta0$m;
  ls1 = theta0$ls;
  s1 = exp(ls1);
  point = c(m1,s1);

  ## Final information
  I = obsInfo(theta=list("m"=m1,"s"=s1),log.scale=F);
  Ii = matInv(I);

  # Standard errors
  se = sqrt(diag(Ii));

  ## Parameter frame
  P = data.frame(c("Location","Scale"),point,se);
  colnames(P) = c("Aspect","Estimate","SE");

  ## Fitted Distribution
  # Estimate mean
  mu = exp(m1+s1^2/2);
  dg = mu*c(1,s1);
  se.mu = sqrt(as.numeric(matQF(X=dg,A=Ii)));

  # Estimate median
  me = exp(m1);
  dg = c(me,0);
  se.me = sqrt(as.numeric(matQF(X=dg,A=Ii)));

  # Estimate variance
  v = (exp(s1^2)-1)*exp(2*m1+s1^2);
  dg = 2*exp(2*m1+s1^2)*c(exp(s1^2)-1,s1*(2*exp(s1^2)-1));
  se.v = sqrt(as.numeric(matQF(X=dg,A=Ii)));

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
  S = function(t){pnorm(q=log(t),mean=m1,sd=s1,lower.tail=F)};

  ## Format Results
  Out = new(Class="fit",Distribution="Log-Normal",Parameters=P,Information=I,Outcome=Y,S=S);
  return(Out);
}
