# Purpose: Estimation of Gamma distribution
# Updated: 20/03/08

# Notes
# Switched to log-scale parameters for Newton Raphson
# Numerically verified

########################
# Gamma Distribution 
########################

#' Gamma Distribution Parameter Estimation
#'
#' Estimates parameters for gamma event times subject to non-informative
#' right censoring. The gamma distribution is parameterized in terms
#' of the shape \eqn{\alpha} and rate \eqn{\lambda}:
#' \deqn{f(t) = \frac{\lambda}{\Gamma(\alpha)}(\lambda t)^{\alpha-1}e^{-\lambda t}, t>0}
#'
#' @param time Numeric observation times.
#' @param status Status indicator, coded as 1 if an event was observed, 0 if censored.
#' @param sig Significance level, for CIs.
#' @param tau Optional truncation times for calculating RMSTs.
#' @param init Vector of initial values for shape \eqn{\alpha} and
#'   rate \eqn{\lambda} parameters, respectively.
#' @param eps Tolerance for Newton-Raphson iterations.
#' @param maxit Maximum number of NR iterations.
#' @param report Report fitting progress?
#'
#' @importFrom expint gammainc
#' @importFrom methods new
#' @importFrom numDeriv grad hessian
#' @importFrom stats qgamma var
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
#' # Generate Gamma data with 20% censoring. 
#' data = genData(n=1e3,dist="gamma",theta=c(2,2),p=0.2);
#' # Estimate
#' fit = fitParaSurv(time=data$time,status=data$status,dist="gamma");

fit.Gamma = function(time,status,sig=0.05,tau=NULL,init=NULL,eps=1e-6,maxit=10,report=FALSE){
  # Events
  n = length(time);
  nobs = sum(status);
  # Partition times
  tobs = time[status==1];
  tcen = time[status==0];
  # Presence of censored events
  flag = (length(tcen)>0);
  # Reusable sums
  sum.tobs = sum(tobs);
  sum.log.tobs = sum(log(tobs));
  
  # Fitting procedure is simplified if complete data are available. Each branch
  # should return: estimated theta, and a function for calculating observed
  # information.
  
  if(!flag){
    
    # MLEs of theta
    theta = fit.Gamma.Complete(time=time,eps=eps);
    
    # Score
    ## Verified numerically. 
    score = function(theta,log.scale){
      a = theta[1];
      l = theta[2];
      if(log.scale){
        a = exp(a);
        l = exp(l);
      }
      # Calculate score
      ua = n*log(l)+sum.log.tobs-n*digamma(a);
      ul = n*a/l-sum.tobs;
      out = c(ua,ul);
      # Map to log scale
      if(log.scale){
        out = out*c(a,l);
      }
      # Output
      return(out);
    }
    
    # Observed information
    ## Verified numerically 
    obsInfo = function(theta,log.scale){
      # Extract parameters
      a = theta[1];
      l = theta[2];
      if(log.scale){
        a = exp(a);
        l = exp(l);
      }
      # Calculate information
      Jaa = n*trigamma(a);
      Jll = nobs*a/(l^2);
      Jla = -(nobs/l);
      # Observed info
      Out = matrix(c(Jaa,Jla,Jla,Jll),nrow=2);
      # Map to log scale
      if(log.scale){
        v = matrix(c(a,l),ncol=1);
        Out = Out*matOP(v,v)-diag(score(theta=theta,log.scale=TRUE));
      }
      # Output
      return(Out);
    };
    
  } else {
    
    # Log lkelihood
    ll = function(theta,log.scale){
      out = survLogLik(time=time,status=status,dist="gamma",theta=theta,log.scale=log.scale);
      return(as.numeric(out));
    }
    
    # Numeric evaluation of score and hessian is faster than analytic due 
    # to the presence of incomplete gamma functions. 
    
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
    a0 = init[1];
    l0 = init[2];
    theta0 = list();
    ## If not specified, apply ML to observed data
    if(!(is.numeric(a0)&is.numeric(l0))){
      theta0$theta = log(fit.Gamma.Complete(time=tobs));
    } else {
      theta0$theta = log(c(a0,l0));
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
    theta = exp(theta0$theta);
    rm(theta0,update);
  }
  
  # Final parameters
  ## Log scale
  log.theta = log(theta);
  
  # Information 
  ## Log scale
  J = obsInfo(log.theta,log.scale=TRUE);
  Ji = matInv(J);
  ## Standard scale
  I = obsInfo(theta,log.scale=FALSE);
  Ii = matInv(I);
  
  if(any(diag(Ji)<0)|any(diag(Ii)<0)){
    stop("Information matrix not positive definite. Try another initialization");
  }

  ## Standard errors
  log.se = sqrt(diag(Ji));
  se = sqrt(diag(Ii));

  # Parameter frame
  P = data.frame(c("Shape","Rate"),theta,se);
  colnames(P) = c("Aspect","Estimate","SE");

  # Fitted Distribution
  a = theta[1];
  l = theta[2];
  ## Estimate mean
  mu = a/l;
  dg = c(1/l,-1*a/(l^2));
  se.mu = sqrt(as.numeric(matQF(X=matrix(dg,ncol=1),A=Ii)));

  # Estimate median
  g = function(x){qgamma(p=0.5,shape=x[1],rate=x[2])};
  me = g(theta);
  dg = grad(func=g,x=theta);
  se.me = sqrt(as.numeric(matQF(X=matrix(dg,ncol=1),A=Ii)));

  # Estimate variance
  v = a/(l^2);
  dg = c(1/(l^2),-2*a/(l^3));
  se.v = sqrt(as.numeric(matQF(X=matrix(dg,ncol=1),A=Ii)));

  # Distribution characteristics
  Y = data.frame(c("Mean","Median","Variance"),
                 c(mu,me,v),c(se.mu,se.me,se.v),
                 stringsAsFactors=FALSE);
  colnames(Y) = c("Aspect","Estimate","SE");
  

  # CIs
  z = qnorm(1-sig/2);
  P$L = exp(log.theta-z*log.se);
  P$U = exp(log.theta+z*log.se);
  Y$L = Y$Estimate-z*Y$SE;
  Y$U = Y$Estimate+z*Y$SE;

  # Fitted survival function
  S = function(t){expint::gammainc(a,l*t)/gamma(a)};

  ## Format Results
  Out = new(Class="fit",Distribution="gamma",Parameters=P,Information=I,Outcome=Y,S=S);
  
  ## Add RMST if requested. 
  if(is.numeric(tau)){
    rmst = paraRMST(fit=Out,sig=sig,tau=tau);
    Out@RMST = rmst;
  }
  
  return(Out);
}

########################
# Gamma Distribution without Censoring
########################

#' Gamma Distribution Parameter Estimation without Censoring
#'
#' Estimates parameters for gamma event times not subject to censoring. 
#'
#' @param time Observation times.
#' @param eps Tolerance for Newton-Raphson iterations.

fit.Gamma.Complete = function(time,eps=1e-6){
  
  # Events
  n = length(time);
  # Reusable sums
  sum.t = sum(time);
  sum.log.t = sum(log(time));
  
  # Profile score
  pscore = function(a){
    if(a==0){
      out = Inf;
    } else {
      out = -n*digamma(a)-n*log(sum.t)+n*log(n*a)+sum.log.t;
    }
    return(out);
  }
  
  # Moment estimators
  m = mean(time);
  v = var(time);
  a0 = m^2/v;
  
  # Numerically zero profile score
  a = uniroot(f=pscore,lower=0,upper=2*a0,extendInt="downX",tol=eps)$root;
  l = n*a/sum.t;
  theta = c("a"=a,"l"=l);
  return(theta);
}





