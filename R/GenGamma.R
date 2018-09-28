# Purpose: Estimation of Gamma distribution
# Updated: 180928

# Notes
# Switched to log-scale parameters
# Numerically verified: Score, Hessian

########################
# Initialize Generalized Gamma
########################

#' Initialization for Generalized Gamma
#'
#' Initializes the parameters for the generalized gamma distribution using
#' moment estimators.
#'
#' @param tobs Observed event times.
#' @param L Lower bound on shape parameter \eqn{\alpha}.
#' @param U Upper bound on shape parameter \eqn{\alpha}.
#'
#' @return List containing the estimated shape and rate parameters on log scale.
#'
#' @importFrom stats optim

init.GenGamma = function(tobs,L=0.01,U=10){
  # Log time
  v = log(tobs);
  # Mean
  m1 = mean(v);
  # Higher central moments
  m2 = m3 = m4 = NULL;
  for(k in 2:4){
    assign(paste0("m",k),mean((v-m1)^k));
  }
  # Objective function for alpha
  Qa = function(a){
    # Derivatives of log gamma
    p0 = psigamma(a,deriv=0);
    p1 = psigamma(a,deriv=1);
    p2 = psigamma(a,deriv=2);
    p3 = psigamma(a,deriv=3);
    ## Third moment estimator
    # Observed value
    obs1 = 2*log(abs(m3))-3*log(m2);
    # Theoritcal value
    exp1 = 2*log(abs(p2))-3*log(p1);
    ## Fourth moment estimator
    # Observed value
    obs2 = log(m4)-2*log(m2);
    # Theoretical value
    exp2 = log(p3+3*(p1^2))-2*log(p1);
    ## GMM Objective
    Out = (obs1-exp1)^2+(obs2-exp2)^2;
    return(Out);
  }
  # Estimate for alpha
  a0 = optim(par=c(1),fn=Qa,method="Brent",lower=L,upper=U)$par;
  # Derivatives of log gamma
  p0 = psigamma(a0,deriv=0);
  p1 = psigamma(a0,deriv=1);
  p2 = psigamma(a0,deriv=2);
  p3 = psigamma(a0,deriv=3);
  # Estimates for beta
  b01 = (p1/m2)^(1/2);
  b02 = (p2/m3)^(1/3);
  b03 = ((p3+3*p1^2)/m4)^(1/4);
  # Combine estimates
  b0 = mean(c(b01,b02,b03));
  # Estimate log lambda
  n = length(tobs);
  ll0 = (1/b0)*(log(n*a0)-log(sum(tobs^b0)));
  # Output
  Out = list("la"=log(a0),"lb"=log(b0),"ll"=ll0);
  return(Out);
}

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
#' @param time Observation times.
#' @param status Status indicator, coded as 1 if an event was observed, 0 if censored.
#' @param sig Significance level, for CIs.
#' @param init List of initial parameter values, including the logs of the shape
#'   parameters "la" and "lb", and the log of the rate parameter "ll".
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
#' }
#'
#' @seealso
#' \itemize{
#'   \item{}{Fitting function for parametric survival distributions \code{\link{fitParaSurv}}}
#' }
#'
#' @examples
#' set.seed(100);
#' # Simulate
#' D = rGenGamma(n=1e4,a=2,b=2,l=2);
#' # Estimate
#' M = fitParaSurv(time=D$time,status=D$status,dist="gengamma");

fit.GenGamma = function(time,status,sig=0.05,init=NULL,eps=1e-6,maxit=10,report=F){
  # Input check
  if(!is.null(init)){
    if(!all(names(init)%in%c("la","lb","ll"))){stop("For gengamma, init should contain log shapes 'la' and 'lb', and log rate 'll'.")};
  }
  # Vectorize
  dLogIncGamma = Vectorize(dLogIncGamma,vectorize.args="s");
  # Events
  n = length(time);
  nobs = sum(status);
  # Observed events
  tobs = time[status==1];
  tcen = time[status==0];
  # Presence of censored events
  flag = (length(tcen)>0);
  # Reusable quantities
  log.tobs = log(tobs);
  sum.log.tobs = sum(log.tobs);

  ## Score
  Score = function(theta,log.scale=T){
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
    ## Calculate score
    # Reusable components
    tobs.b = tobs^b;
    # Score for alpha
    Ua = -1*n*digamma(a)+nobs*b*log(l)+b*sum.log.tobs;
    # Score for beta
    Ub = nobs/b+nobs*a*log(l)+a*sum.log.tobs-(l^b)*log(l)*sum(tobs.b)-(l^b)*sum((tobs.b)*log.tobs);
    # Score for lambda
    Ul = nobs*a*b/l-b*(l^(b-1))*sum(tobs.b);
    # Add corrections for censoring
    if(flag){
      ## Differentials of (lt)^b
      D1b = (l*tcen)^b*log(l*tcen);
      D1l = b*tcen*((l*tcen)^(b-1));
      # Updates for censoring
      Ua = Ua+sum(dLogIncGamma(a=a,s=(l*tcen)^b,Dirn="a",Order=1));
      Ub = Ub+sum(dLogIncGamma(a=a,s=(l*tcen)^b,Dirn="s")*D1b);
      Ul = Ul+sum(dLogIncGamma(a=a,s=(l*tcen)^b,Dirn="s")*D1l);
    }
    # Score
    Out = c(Ua,Ub,Ul);
    # Map to log scale
    if(log.scale){
      Out = Out*c(a,b,l);
    }
    return(Out);
  };

  ## Observed information
  obsInfo = function(theta,log.scale=T){
    # Extract parameters
    if(log.scale){
      a = exp(theta$la);
      b = exp(theta$lb);
      l = exp(theta$ll);
      v = c(a,b,l);
    } else {
      a = theta$a;
      b = theta$b;
      l = theta$l;
    }
    ## Calculate hessian
    # Reusable components
    tobs.b = tobs^b;
    ## Hessian components
    # For alpha
    Haa = -n*trigamma(a);
    # For beta
    Hbb = -nobs/(b^2)-(l^b)*(log(l)^2)*sum(tobs.b)-2*(l^b)*log(l)*sum((tobs.b)*log.tobs)-(l^b)*sum((tobs.b)*(log.tobs^2));
    # For lambda
    Hll = -nobs*a*b/(l^2)-b*(b-1)*(l^(b-2))*sum(tobs.b);
    # For alpha-beta
    Hab = nobs*log(l)+sum.log.tobs;
    # For alpha-lambda
    Hal = nobs*b/l;
    # For beta-lambda
    Hbl = nobs*a/l-(l^(b-1))*sum(tobs.b)*(1+b*log(l))-b*(l^(b-1))*sum((tobs.b)*log.tobs);
    # Add corrections for censoring
    if(flag){
      # Reusable components
      z = (l*tcen)^b;
      ## Differentials of (lt)^b;
      D1b = z*log(l*tcen);
      D2b = z*(log(l*tcen)^2);
      D1l = b*tcen*((l*tcen)^(b-1));
      D2l = b*(b-1)*(tcen^2)*((l*tcen)^(b-2));
      Dbl = b*tcen*((l*tcen)^(b-1))*log(l*tcen)+tcen*((l*tcen)^(b-1));
      ## Updates for censoring
      # For alpha
      Haa = Haa+sum(dLogIncGamma(a=a,s=z,Dirn="a",Order=2));
      # For beta
      Hbb = Hbb+sum(dLogIncGamma(a=a,s=z,Dirn="s",Order=2)*(D1b)^2);
      Hbb = Hbb+sum(dLogIncGamma(a=a,s=z,Dirn="s",Order=1)*D2b);
      # For lambda
      Hll = Hll+sum(dLogIncGamma(a=a,s=z,Dirn="s",Order=2)*(D1l)^2);
      Hll = Hll+sum(dLogIncGamma(a=a,s=z,Dirn="s",Order=1)*D2l);
      # For alpha-beta
      Hab = Hab+sum(dLogIncGamma(a=a,s=z,Dirn="as")*D1b);
      # For alpha-lambda
      Hal = Hal+sum(dLogIncGamma(a=a,s=z,Dirn="as")*D1l);
      # For beta-lambda
      Hbl = Hbl+sum(dLogIncGamma(a=a,s=z,Dirn="s",Order=2)*(D1b*D1l));
      Hbl = Hbl+sum(dLogIncGamma(a=a,s=z,Dirn="s",Order=1)*Dbl);
    }
    # Observed info
    Out = -1*matrix(c(Haa,Hab,Hal,Hab,Hbb,Hbl,Hal,Hbl,Hll),nrow=3,byrow=T);
    # Map to log scale
    if(log.scale){
      Out = Out*matOP(v,v)-diag(Score(theta=theta,log.scale=T));
    }
    return(Out);
  }

  ## NR Update
  Update = function(theta){
    # Initial log likelihood
    q0 = survLogLik(time=time,status=status,dist="gengamma",theta=theta,log.scale=T);
    # Score
    U = Score(theta,log.scale=T);
    # Inverse observed information
    Ji = matInv(obsInfo(theta,log.scale=T));
    # Proposal
    Prop = c(theta$la,theta$lb,theta$ll)+MMP(Ji,U);
    Prop = as.numeric(Prop);
    # Proposed theta
    Out = list("la"=Prop[1],"lb"=Prop[2],"ll"=Prop[3]);
    # Final log likelihood
    q1 = survLogLik(time=time,status=status,dist="gengamma",theta=Out,log.scale=T);
    # Increment
    Out$d = (q1-q0);
    return(Out);
  }

  ## Initialize
  theta0 = list();
  # Inputs
  la0 = init$la;
  lb0 = init$lb;
  ll0 = init$ll;
  miss = is.null(la0)|is.null(lb0)|is.null(ll0);
  # Attempt to initialize
  if(miss){theta0 = init.GenGamma(tobs=tobs)};
  # Overwrite supplied values
  if(!is.null(la0)){theta0$la = la0};
  if(!is.null(lb0)){theta0$lb = lb0};
  if(!is.null(ll0)){theta0$ll = ll0};

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
  lb1 = theta0$lb;
  ll1 = theta0$ll;
  log.point = c(la1,lb1,ll1);

  # Standard scale
  a1 = exp(la1);
  b1 = exp(lb1);
  l1 = exp(ll1);
  point = c(a1,b1,l1);

  ## Final Information
  # Log scale
  J = obsInfo(theta0);
  # Original scale
  theta1 = list("a"=a1,"b"=b1,"l"=l1);
  I = obsInfo(theta1,log.scale=F);

  ## Check positivity
  # Eigen-decomposition
  eJ = eigen(x=J,symmetric=T);
  eI = eigen(x=I,symmetric=T);
  # If not positive denfiite, try robust SEs
  if(min(eJ$values,eI$values)<0){
    warning("Observed information was not positive definite. Consider another parameter initialization.");
    if(report){
      cat("Eigenvalues of log-scale information:\n");
      cat(round(eJ$values,digits=2),"\n");
      cat("Eigenvalues of original-scale information:\n");
      cat(round(eI$values,digits=2),"\n\n");
    };
    # Robust SE, log scale
    u = Score(theta0,log.scale=T);
    B = matOP(u,u);
    Ji = matQF(X=matInv(J),B);
    # Robust SE, standard scale
    u = Score(theta1,log.scale=F);
    B = matOP(u,u);
    Ii = matQF(X=matInv(I),B);
  } else {
    Ji = matInv(J);
    Ii = matInv(I);
  }

  # Standard errors
  log.se = sqrt(diag(Ji));
  se = sqrt(diag(Ii));

  # Parameter frame
  P = data.frame(c("ShapeA","ShapeB","Rate"),point,se);
  colnames(P) = c("Aspect","Estimate","SE");

  ## Outcome
  # Estimate mean
  g = function(x){
    a1 = exp(x[1]); b1 = exp(x[2]); l1 = exp(x[3]);
    return(gamma(a1+1/b1)/(l1*gamma(a1)));
  };
  mu = g(log.point);
  dg = grad(func=g,x=log.point);
  se.mu = sqrt(as.numeric(matQF(X=dg,A=Ji)));

  # Estimate median
  g = function(x){
    a1 = exp(x[1]); b1 = exp(x[2]); l1 = exp(x[3]);
    return(qgamma(p=0.5,shape=a1,rate=l1^b1)^(1/b1));
  };
  me = g(log.point);
  dg = grad(func=g,x=log.point);
  se.me = sqrt(as.numeric(matQF(X=dg,A=Ji)));

  # Estimate variance
  g = function(x){
    a1 = exp(x[1]); b1 = exp(x[2]); l1 = exp(x[3]);
    return(1/(l1^2*gamma(a1))*(gamma(a1+2/b1)-gamma(a1+1/b1)^2/gamma(a1)));
  }
  v = g(log.point);
  dg = grad(func=g,x=log.point);
  se.v = sqrt(as.numeric(matQF(X=dg,A=Ji)));

  # Distribution characteristics
  Y = data.frame(c("Mean","Median","Variance"),c(mu,me,v),c(se.mu,se.me,se.v));
  colnames(Y) = c("Aspect","Estimate","SE");

  # CIs
  z = qnorm(1-sig/2);
  P$L = exp(log.point-z*log.se);
  P$U = exp(log.point+z*log.se);
  Y$L = Y$Estimate-z*Y$SE;
  Y$U = Y$Estimate+z*Y$SE;

  ## Format Results
  Out = new(Class="fit",Distribution="GenGamma",Parameters=P,Information=I,Outcome=Y);
  return(Out);
}
