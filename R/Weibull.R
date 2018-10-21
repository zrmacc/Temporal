# Purpose: Estimation of Weibull distribution
# Updated: 180926

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
#' @param time Observation times.
#' @param status Status indicator, coded as 1 if an event was observed, 0 if censored.
#' @param sig Significance level, for CIs.
#' @param init List containing an initial value for the shape parameter "a".
#'
#' @importFrom methods new
#' @importFrom stats quantile uniroot
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
#'   \item{Fitting function for parametric survival distributions \code{\link{fitParaSurv}}}
#' }
#'
#' @examples
#' # Simulate
#' D = rWeibull(n=1e3,a=2,l=2);
#' # Estimate
#' M = fitParaSurv(time=D$time,status=D$status,dist="weibull");

fit.Weibull = function(time,status,sig=0.05,init=NULL){
  # Input check
  if(!is.null(init)){
    if(!all(names(init)%in%c("a"))){stop("For weibull, init should contain shape 'a'.")};
  }

  # Events
  n = length(time);
  nobs = sum(status);
  # Observed events
  tobs = time[status==1];
  # Log times
  logtime = log(time);
  logtobs = log(tobs);

  # Profile score for shape
  Score = function(a){
    if(a==0){
      Out = Inf;
    }  else {
      # Scaled time
      ta = time^a;
      # Score
      Out = nobs/a-nobs*sum(ta*logtime)/sum(ta)+sum(logtobs);
    }
    return(Out);
  }

  ## Initialize
  a0 = init$a0;
  if(is.null(a0)){
    q0 = as.numeric(quantile(x=tobs,probs=c(1-exp(-1))));
    l0 = 1/q0;
    a0 = digamma(1)/(log(l0)+mean(logtobs));
  } else {
    a0 = a0;
  }
  # Search for upper limit
  U = a0;
  fU = Score(U);
  sU = sign(fU);
  if(sU==1){
    while(sU==1){
      U = 2*U;
      fU = Score(U);
      sU = sign(fU);
    }
  }

  ## Optimize
  # MLE for alpha
  a1 = uniroot(f=Score,lower=0,upper=U,f.upper=fU)$root;
  # MLE for lambda
  ta1 = time^a1;
  l1 = (sum(ta1)/nobs)^(-1/a1);

  ## Information
  S1 = sum(ta1);
  S2 = sum(ta1*logtime);
  S3 = sum(ta1*(logtime)^2);
  # Information components
  Jaa = (nobs/(a1^2))+(l1^a1)*(log(l1))^2*S1+2*(l1^a1)*log(l1)*S2+(l1^a1)*S3;
  Jll = (nobs*a1)/(l1^2)+a1*(a1-1)*(l1^(a1-2))*S1;
  Jla = -(nobs/l1)+(l1^(a1-1))*S1+a1*(l1^(a1-1))*log(l1)*S1+a1*(l1^(a1-1))*S2;
  # Information matrix
  J = matrix(c(Jaa,Jla,Jla,Jll),nrow=2);
  colnames(J) = rownames(J) = c("a","l");
  Ji = matInv(J);

  # Parameters
  P = data.frame(c("Shape","Rate"),c(a1,l1),sqrt(diag(Ji)));
  colnames(P) = c("Aspect","Estimate","SE");

  ## Fitted Distribution
  # Estimate mean
  mu = (1/l1)*gamma(1+1/a1);
  dg = -(1/l1)*gamma(1+1/a1)*c(digamma(1+1/a1)/a1^2,1/l1);
  se.mu = sqrt(as.numeric(matQF(dg,Ji)));

  # Estimate median
  me = (1/l1)*(log(2))^(1/a1);
  dg = -(1/l1)*(log(2))^(1/a1)*c(log(log(2))/a1^2,(1/l1));
  se.me = sqrt(as.numeric(matQF(dg,Ji)));

  # Estimate variance
  v = (1/l1^2)*(gamma(1+2/a1)-gamma(1+1/a1)^2);
  dg = -(2/l1^2)*c((1/a1^2)*(gamma(1+2/a1)*digamma(1+2/a1)-(gamma(1+1/a1)^2)*digamma(1+1/a1)),
                   (1/l1)*(gamma(1+2/a1)-gamma(1+1/a1)^2));
  se.v = sqrt(as.numeric(matQF(dg,Ji)));

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
  S = function(t){exp(-(l1*t)^a1)};

  ## Format Results
  Out = new(Class="fit",Distribution="Weibull",Parameters=P,Information=J,Outcome=Y,S=S);
  return(Out);
}
