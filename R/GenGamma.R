# Purpose: Estimation of Gamma distribution
# Updated: 20/03/11

# Notes
# Switched to log-scale parameters for Newton Raphson
# Numerically verified

########################
# Generalized Gamma Distribution
########################

#' Generalized Gamma Distribution Parameter Estimation
#'
#' Estimates parameters for generalized gamma event times subject to non-informative
#' right censoring. The gamma distribution is parameterized in terms
#' of the shape parameters \eqn{(\alpha,\beta)}, and the rate \eqn{\lambda}:
#' \deqn{f(t) = \frac{\beta\lambda}{\Gamma(\alpha)} (\lambda t)^{\alpha\beta-1}e^{-(\lambda t)^{\beta}}, t>0}
#'
#' @param time Numeric observation times.
#' @param status Status indicator, coded as 1 if an event was observed, 0 if censored.
#' @param sig Significance level, for CIs.
#' @param tau Optional truncation times for calculating RMSTs.
#' @param init Vector of initial values for the shape \eqn{\alpha}, \eqn{\beta}, and
#'   rate \eqn{\lambda} parameters, respectively. 
#' @param bL Lower limit on possible values for beta.
#' @param bU Upper limit on possible values for beta. 
#' @param eps Tolerance for Newton-Raphson iterations.
#' @param maxit Maximum number of NR iterations.
#' @param report Report fitting progress?
#'
#' @importFrom expint gammainc
#' @importFrom methods new
#' @importFrom numDeriv grad
#' @importFrom stats qgamma var
#'
#' @return An object of class \code{fit} containing the following:
#' \describe{
#'  \item{Parameters}{The estimated shape \eqn{(\alpha,\beta)} and rate \eqn{\lambda} parameters.}
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
#' set.seed(103);
#' # Generalized Gamma data with 20% censoring.
#' data = genData(n=1e4,dist="gen-gamma",theta=c(2,2,2),p=0.2);
#' # Estimate
#' fit = fitParaSurv(time=data$time,status=data$status,dist="gen-gamma",report=TRUE);

fit.GenGamma = function(time,status,sig=0.05,tau=NULL,init,bL=0.1,bU=10,eps=1e-6,maxit=10,report=FALSE){
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
    theta = fit.GenGamma.Complete(time=tobs,bL=bL,bU=bU);
    
    # Score
    ## Verified numerically. 
    score = function(theta,log.scale){
      a = theta[1];
      b = theta[2];
      l = theta[3];
      if(log.scale){
        a = exp(a);
        b = exp(b);
        l = exp(l);
      }
      # Calculate score
      ua = n*b*log(l)+b*sum.log.tobs-n*digamma(a);
      ub = (n/b)+n*a*log(l)+a*sum.log.tobs-(l^b)*log(l)*sum(tobs^b)-(l^b)*sum((tobs^b)*log(tobs));
      ul = n*a*b/l-b*(l^(b-1))*sum(tobs^b);
      out = c(ua,ub,ul);
      # Map to log scale
      if(log.scale){
        out = out*c(a,b,l);
      }
      # Output
      return(out);
    }
    
    # Observed information
    ## Verified numerically 
    obsInfo = function(theta,log.scale){
      # Extract parameters
      a = theta[1];
      b = theta[2];
      l = theta[3];
      if(log.scale){
        a = exp(a);
        b = exp(b);
        l = exp(l);
      }
      # Reusable terms
      lb = l^b;
      lb1 = l^(b-1);
      lb2 = l^(b-2);
      ll = log(l);
      sum.t.b = sum(tobs^b);
      # Calculate hessian
      Haa = -n*trigamma(a);
      Hab = n*ll+sum.log.tobs;
      Hbb = -n/(b^2)-lb*ll^2*sum(tobs^b)-2*lb*ll*sum((tobs^b)*log(tobs))-lb*sum((tobs^b)*(log(tobs)^2));
      Hbl = n*a/l-b*lb1*ll*sum.t.b-lb1*sum.t.b-b*lb1*sum((tobs^b)*log(tobs));
      Hll = -n*a*b/(l^2)-b*(b-1)*lb2*sum.t.b;
      Hal = n*b/l;
      # Observed info
      Out = -1*rbind(c(Haa,Hab,Hal),
                     c(Hab,Hbb,Hbl),
                     c(Hal,Hab,Hll));
      # Map to log scale
      if(log.scale){
        v = matrix(c(a,b,l),ncol=1);
        Out = Out*matOP(v,v)-diag(score(theta=theta,log.scale=TRUE));
      }
      # Output
      return(Out);
    };
    
  } else {
    
    # Log lkelihood
    ll = function(theta,log.scale){
      out = survLogLik(time=time,status=status,dist="gen-gamma",theta=theta,log.scale=log.scale);
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
    b0 = init[2];
    l0 = init[3];
    theta0 = list();
    ## If not specified, apply ML to observed data
    if(!(is.numeric(a0)&is.numeric(b0)&is.numeric(l0))){
      theta0$theta = log(fit.GenGamma.Complete(time=tobs,bL=bL,bU=bU));
    } else {
      theta0$theta = log(c(a0,b0,l0));
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

  # Standard errors
  log.se = sqrt(diag(Ji));
  se = sqrt(diag(Ii));

  # Parameter frame
  P = data.frame(c("ShapeA","ShapeB","Rate"),theta,se);
  colnames(P) = c("Aspect","Estimate","SE");

  ## Outcome
  # Estimate mean
  g = function(x){
    a = exp(x[1]); b = exp(x[2]); l = exp(x[3]);
    return(gamma(a+1/b)/(l*gamma(a)));
  };
  mu = g(log.theta);
  dg = grad(func=g,x=log.theta);
  se.mu = sqrt(as.numeric(matQF(X=matrix(dg,ncol=1),A=Ji)));

  # Estimate median
  g = function(x){
    a = exp(x[1]); b = exp(x[2]); l = exp(x[3]);
    return(qgamma(p=0.5,shape=a,rate=l^b)^(1/b));
  };
  me = g(log.theta);
  dg = grad(func=g,x=log.theta);
  se.me = sqrt(as.numeric(matQF(X=matrix(dg,ncol=1),A=Ji)));

  # Estimate variance
  g = function(x){
    a = exp(x[1]); b = exp(x[2]); l = exp(x[3]);
    return(1/(l^2*gamma(a))*(gamma(a+2/b)-gamma(a+1/b)^2/gamma(a)));
  }
  v = g(log.theta);
  dg = grad(func=g,x=log.theta);
  se.v = sqrt(as.numeric(matQF(X=matrix(dg,ncol=1),A=Ji)));

  # Distribution characteristics
  Y = data.frame(c("Mean","Median","Variance"),c(mu,me,v),c(se.mu,se.me,se.v));
  colnames(Y) = c("Aspect","Estimate","SE");

  # CIs
  z = qnorm(1-sig/2);
  P$L = exp(log.theta-z*log.se);
  P$U = exp(log.theta+z*log.se);
  Y$L = Y$Estimate-z*Y$SE;
  Y$U = Y$Estimate+z*Y$SE;

  # Fitted survival function
  a = theta[1];
  b = theta[2];
  l = theta[3];
  S = function(t){expint::gammainc(a,(l*t)^b)/gamma(a)};

  # Format Results
  Out = new(Class="fit",Distribution="gen-gamma",Parameters=P,Information=I,Outcome=Y,S=S);
  
  ## Add RMST if requested. 
  if(is.numeric(tau)){
    rmst = paraRMST(fit=Out,sig=sig,tau=tau);
    Out@RMST = rmst;
  }
  
  return(Out);
}

########################
# Initialize Generalized Gamma
########################

#' Initialization for Generalized Gamma
#'
#' Initializes the parameters for the generalized gamma distribution via maximum
#' likelihood.
#'
#' @param time Observed event times.
#' @param bL Lower limit on possible values for beta.
#' @param bU Upper limit on possible values for beta.
#'
#' @return 3x1 numeric vector of estimated parameters, \eqn{\alpha},
#'   \eqn{\beta}, \eqn{\lambda}.
#'
#' @importFrom stats optimize

fit.GenGamma.Complete = function(time,bL,bU){
  
  # Events
  n = length(time);
  sum.log.t = sum(log(time));
  
  # Rate parameter
  rate = function(a,b){
    out = 1/(n*a)*sum(time^b);
    out = 1/(out^(1/b));
    return(out);
  }
  
  # Shape 1
  shape1 = function(b){
    out = sum((time^b)*log(time))/sum(time^b)-sum.log.t/n;
    out = 1/(b*out);
    return(out);
  }
  
  # Profile log likelihood
  ll = function(b){
    a = shape1(b);
    l = rate(a,b);
    out = survLogLik(time=time,theta=c(a,b,l),dist="gen-gamma");
    return(out);
  }
  
  b0 = optimize(f=ll,lower=bL,upper=bU,maximum=TRUE)$maximum;
  
  # Remaining MLEs
  a0 = shape1(b0);
  l0 = rate(a0,b0);
  theta = c("a"=a0,"b"=b0,"l"=l0);
  return(theta);
}
