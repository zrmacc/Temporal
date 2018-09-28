# Purpose: Estimation of Gamma distribution
# Updated: 180928

# Notes
# Converted parameters to log-scale for optimization
# Numerically verified: Score, Hessian, Mean, Median, Var

########################
# Gamma Distribution
########################

#' Gamma Distribution Parameter Estimation
#'
#' Estimates parameters for gamma event times subject to non-informative
#' right censoring. The gamma distribution is parameterized in terms
#' of the shape \eqn{\alpha} and rate \eqn{\lambda}:
#' \deqn{f(t) = \frac{\lambda{\Gamma(\alpha)} (\lambda t)^{\alpha-1}e^{-\lambda t}, t>0}}
#'
#' @param time Observation times.
#' @param status Status indicator, coded as 1 if an event was observed, 0 if
#'   censored.
#' @param sig Significance level, for CIs.
#' @param init List of initial parameter values, including the log of the shape
#'   parameter "la" and the log of the rate parameter "ll".
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
#'  \item{Parameters}{The estimated shape \eqn{\alpha} and rate \eqn{\lambda}.}
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
#' # Simulate
#' D = rGamma(n=1e3,a=2,l=2);
#' # Estimate
#' M = fitParaSurv(time=D$time,status=D$status,dist="gamma");

fit.Gamma = function(time,status,sig=0.05,init=NULL,eps=1e-6,maxit=10,report=F){
  # Input check
  if(!is.null(init)){
    if(!all(names(init)%in%c("la","ll"))){stop("For gamma, init should contain log shape 'la' and log rate 'll'.")};
  }
  # Vectorize
  dLogIncGamma = Vectorize(dLogIncGamma,vectorize.args="s");
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

  ## Function to calculate score
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
    # Score for alpha
    Ua = nobs*log(l)+sum.log.tobs-n*digamma(a);
    # Score for lambda
    Ul = nobs*a/l-sum.tobs;
    # Add corrections for censoring
    if(flag){
      Ua = Ua+sum(dLogIncGamma(a=a,s=l*tcen,Dirn="a",Order=1));
      Ul = Ul+sum(dLogIncGamma(a=a,s=l*tcen,Dirn="s",Order=1)*tcen);
    }
    # Score
    Out = c(Ua,Ul);
    # Map to log scale
    if(log.scale){
      Out = Out*c(a,l);
    }
    # Output
    return(Out);
  };

  ## Function to calculate observed information
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
    # Information for alpha
    Jaa = n*trigamma(a);
    # Information for lambda
    Jll = nobs*a/(l^2);
    # Cross information
    Jla = -(nobs/l);
    # Add corrections for censoring
    if(flag){
      Jaa = Jaa-sum(dLogIncGamma(a=a,s=l*tcen,Dirn="a",Order=2));
      Jll = Jll-sum(dLogIncGamma(a=a,s=l*tcen,Dirn="s",Order=2)*(tcen^2));
      Jla = Jla-sum(dLogIncGamma(a=a,s=l*tcen,Dirn="as")*tcen);
    }
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
    q0 = survLogLik(time=time,status=status,dist="gamma",theta=theta,log.scale=T);
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
    q1 = survLogLik(time=time,status=status,dist="gamma",theta=Out);
    # Increment
    Out$d = (q1-q0);
    return(Out);
  }

  ## Initialize
  theta0 = list();
  # Shape
  la0 = init$la;
  if(is.null(la0)){theta0$la = 2*log(mean(tobs))-log(var(tobs));} else {theta0$la=la0};
  # Rate
  ll0 = init$ll;
  if(is.null(ll0)){theta0$ll = log(mean(tobs))-log(var(tobs));} else {theta0$ll=ll0};

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
  mu = a1/l1;
  dg = c(1/l1,-1*a1/(l1^2));
  se.mu = sqrt(as.numeric(matQF(X=dg,A=Ii)));

  # Estimate median
  g = function(x){qgamma(p=0.5,shape=x[1],rate=x[2])};
  me = g(point);
  dg = grad(func=g,x=point);
  se.me = sqrt(as.numeric(matQF(X=dg,A=Ii)));

  # Estimate variance
  v = a1/(l1^2);
  dg = c(1/(l1^2),-2*a1/(l1^3));
  se.v = sqrt(as.numeric(matQF(X=dg,A=Ii)));

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
  Out = new(Class="fit",Distribution="Gamma",Parameters=P,Information=I,Outcome=Y);
  return(Out);
}
