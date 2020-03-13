# Purpose: Estimation of log-normal distribution
# Updated: 20/03/08

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
#' @param time Numeric observation times.
#' @param status Status indicator, coded as 1 if an event was observed, 0 if
#'   censored.
#' @param sig Significance level, for CIs.
#' @param tau Optional truncation times for calculating RMSTs. 
#' @param init Vecotor of initial values for the location \eqn{\mu} and scale \eqn{\sigma}. 
#' @param eps Tolerance for Newton-Raphson iterations.
#' @param maxit Maximum number of NR iterations.
#' @param report Report fitting progress?
#'
#' @importFrom methods new
#' @importFrom stats sd var
#'
#' @return An object of class \code{fit} containing the following:
#' \describe{
#'  \item{Parameters}{The estimated location \eqn{\mu} and scale \eqn{\sigma}.}
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
#' # Generate Log-Normal data with 20% censoring. 
#' data = genData(n=1e3,dist="log-normal",theta=c(0,2),p=0.2);
#' # Estimate
#' fit = fitParaSurv(time=data$time,status=data$status,dist="log-normal");

fit.LogNormal = function(time,status,sig=0.05,tau=NULL,init=NULL,eps=1e-6,maxit=10,report=F){
  # Events
  n = length(time);
  nobs = sum(status);
  # Observed events
  tobs = time[status==1];
  tcen = time[status==0];
  # Presence of censored events
  flag = (length(tcen)>0);
  
  # Fitting procedure is simplified if complete data are available. Each branch
  # should return: estimated theta, and a function for calculating observed
  # information.
  
  if(!flag){
    
    # MLEs
    theta = c(mean(log(tobs)),sd(log(tobs)));
    
    # Score
    ## Verified numerically. 
    score = function(theta,log.scale){
      # Extract parameters
      m = theta[1];
      s = theta[2];
      if(log.scale){
        s = exp(s);
      }
      # Calculate score
      zobs = (log(tobs)-m)/s;
      um = (1/s)*sum(zobs);
      us = -1*(nobs/s)+(1/s)*sum(zobs^2);
      out = c(um,us);
      # Map to log scale
      if(log.scale){
        out = out*c(1,s);
      }
      # Output
      return(out);
    }
    
    # Function to calculate observed information
    ## Verified numerically
    obsInfo = function(theta,log.scale){
      # Extract parameters
      m = theta[1];
      s = theta[2];
      if(log.scale){
        s = exp(s);
      }
      # Calculate information
      ## Standardize
      zobs = (log(tobs)-m)/s;
      Jmm = nobs/(s^2);
      Jss = -nobs/(s^2)+3/(s^2)*sum(zobs^2);
      Jms = 2/(s^2)*sum(zobs);
      # Observed info
      out = matrix(c(Jmm,Jms,Jms,Jss),nrow=2);
      # Map to log scale
      if(log.scale){
        v = matrix(c(1,s),ncol=1);
        u = score(theta,log.scale=TRUE);
        out = out*matOP(v,v)-diag(c(0,u[2]));
      }
      # Output
      return(out);
    }
    
  } else {
    
    # Log lkelihood
    ll = function(theta,log.scale){
      out = survLogLik(time=time,status=status,dist="log-normal",theta=theta,log.scale=log.scale);
      return(as.numeric(out));
    }
    
    # Numeric score
    score = function(theta,log.scale){
      aux = function(theta){ll(theta,log.scale=log.scale)};
      out = grad(func=aux,x=theta);
      return(out);
    }
    
    # Numeric observed information
    obsInfo = function(theta,log.scale){
      aux = function(theta){ll(theta,log.scale=log.scale)};
      out = -1*hessian(func=aux,x=theta);
      return(out);
    }
    
    # NR Update
    update = function(input,log.scale){
      theta = input$theta;
      # Initial log likelihood
      ll0 = ll(theta,log.scale=log.scale);
      # Score
      U = score(theta,log.scale=log.scale);
      # Inverse observed information
      Ji = matInv(obsInfo(theta,log.scale=log.scale));
      # Proposal
      prop = theta+as.numeric(MMP(Ji,matrix(U,ncol=1)));
      # Final log likelihood
      ll1 = ll(prop,log.scale=log.scale);
      # Increment
      out = list();
      out$theta = prop;
      out$d = ll1-ll0;
      return(out);
    }
    
    # Initialize
    m0 = init[1];
    s0 = init[2];
    theta0 = list();
    ## If not specified, apply ML to observed data
    if(!(is.numeric(m0)&is.numeric(s0))){
      theta0$theta = c(mean(log(tobs)),log(sd(log(tobs))));
    } else {
      theta0$theta = c(m0,log(s0));
    }
    
    ## Maximzation
    for(i in 1:maxit){
      # Update
      theta1 = update(theta0,log.scale=TRUE);
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
    }
    
    ## Fitting report
    if(report){
      if(i<maxit){
        cat(paste0(i-1," update(s) performed before tolerance limit.\n\n"));
      } else {
        cat(paste0(i," update(s) performed without reaching tolerance limit.\n\n"));
      }
    };
    
    ## Clear workspace
    theta = theta0$theta;
    theta[2] = exp(theta[2]);
    rm(theta0,score,update);
  }

  # Final information
  I = obsInfo(theta=theta,log.scale=FALSE);
  Ii = matInv(I);
  if(any(diag(Ii)<0)){
    stop("Information matrix not positive definite. Try another initialization");
  }

  ## Standard errors
  se = sqrt(diag(Ii));

  # Parameter frame
  P = data.frame(c("Location","Scale"),theta,se);
  colnames(P) = c("Aspect","Estimate","SE");

  # Fitted Distribution
  m = theta[1];
  s = theta[2];
  ## Estimate mean
  mu = exp(m+s^2/2);
  dg = mu*c(1,s);
  se.mu = sqrt(as.numeric(matQF(X=matrix(dg,ncol=1),A=Ii)));

  # Estimate median
  me = exp(m);
  dg = c(me,0);
  se.me = sqrt(as.numeric(matQF(X=matrix(dg,ncol=1),A=Ii)));

  # Estimate variance
  v = (exp(s^2)-1)*exp(2*m+s^2);
  dg = 2*exp(2*m+s^2)*c(exp(s^2)-1,s*(2*exp(s^2)-1));
  se.v = sqrt(as.numeric(matQF(X=matrix(dg,ncol=1),A=Ii)));

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
  S = function(t){pnorm(q=log(t),mean=m,sd=s,lower.tail=F)};

  ## Format Results
  Out = new(Class="fit",Distribution="log-normal",Parameters=P,Information=I,Outcome=Y,S=S);
  
  ## Add RMST if requested. 
  if(is.numeric(tau)){
    rmst = paraRMST(fit=Out,sig=sig,tau=tau);
    Out@RMST = rmst;
  }
  
  return(Out);
}
