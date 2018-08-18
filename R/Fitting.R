# Purpose: Parameter estimation 
# Updated: 180817

#' @useDynLib Temporal
#' @importFrom Rcpp sourceCpp
NULL

########################
# Master Fitting Function
########################

#' Fit Parametric Survival Distribution
#' 
#' Estimates parameters for parametric event times subject to non-informative 
#' right censoring. Available distributions include: exponential, gamma, 
#' log-logistic, log-normal, and Weibull.
#' 
#' @param time Observation times.
#' @param status Status indicator, coded as 1 if an event was observed, 0 if
#'   censored.
#' @param dist Distribution to fit, selected from among: exp, gamma,
#'   log-logistic, log-normal, and weibull.
#' @param alpha Significance level, for CIs. 
#' @param eps Tolerance for Newton-Raphson iterations.
#' @param maxit Maximum number of NR iterations.
#' @export
#' @return An object of class \code{fit} containing the following: \describe{ 
#'   \item{Parameters}{The estimated rate \eqn{\lambda}.} \item{Information}{The
#'   observed information.} \item{Outcome}{The fitted mean, median, and
#'   variance.} }
#' @examples 
#' # Generate cenored gamma data
#' D = rGamma(n=1e3,a=2,l=2,p=0.2);
#' # Fit gamma distribution
#' M = fitParaSurv(time=D$time,status=D$status,dist="gamma");
#' 
#' # Generate cenored weibull data
#' D = rWeibull(n=1e3,a=2,l=2,p=0.2);
#' # Fit weibull distribution
#' M = fitParaSurv(time=D$time,status=D$status,dist="weibull");

fitParaSurv = function(time,status,dist="weibull",alpha=0.05,eps=1e-6,maxit=10){
  ## Input checks
  n = length(time);
  # Positivity
  if(min(time)<0){stop("Strictly positive observation times required.")};
  # Status
  if(missing(status)){status=rep(1,n)};
  status.levels = sort(unique(status));
  if(length(status.levels)==1){
    if(status.levels!=1){status=rep(1,n)};
  }
  if(length(status.levels)==2){
    if(!all.equal(status.levels,c(0,1))){stop("Numeric 0,1 coding is expected for status.")};
  }
  if(length(status.levels)>2){stop("Only two levels are expected for status.")}
  # Distribution
  choices = c("exp","gamma","log-logistic","log-normal","weibull");
  if(!(dist %in% choices)){stop(c("Select distribution from among:\n",paste0(choices,collapse=" ")))};
  
  ## Model fitting
  if(dist=="exp"){
    M = fit.Exp(time=time,status=status,alpha=alpha);
  }
  if(dist=="gamma"){
    M = fit.Gamma(time=time,status=status,alpha=alpha,eps=eps,maxit=maxit);
  }
  if(dist=="log-logistic"){
    M = fit.LogLogistic(time=time,status=status,alpha=alpha,eps=eps,maxit=maxit);
  }
  if(dist=="log-normal"){
    M = fit.LogNormal(time=time,status=status,alpha=alpha,eps=eps,maxit=maxit);
  }
  if(dist=="weibull"){
    M = fit.Weibull(time=time,status=status,alpha=alpha);
  }
  # Output
  return(M);
}

########################
# Exponential Distribution
########################

#' Exponential Parameter Estimation
#'
#' Estimates parameters for exponential event times subject to non-informative
#' right censoring. The exponential distribution is parameterized in terms
#' of the rate \eqn{\lambda}: \deqn{f(t) = \lambda e^{-\lambda t}, t>0}
#'
#' @param time Observation times.
#' @param status Status indicator, coded as 1 if an event was observed, 0 if censored.
#' @param alpha Significance level, for CIs.
#' @importFrom methods new
#' @return An object of class \code{fit} containing the following:
#' \describe{
#'  \item{Parameters}{The estimated model parameters.}
#'  \item{Information}{The observed information matrix.}
#'  \item{Outcome}{The fitted mean, median, and variance of the time to event distribution.}
#' }

fit.Exp = function(time,status,alpha=0.05){
  # Events
  n = length(time);
  nobs = sum(status);
  # MLE of lambda
  l1 = nobs/sum(time);
  # Information
  J = nobs/(l1^2);
  Ji = 1/J;
  # Parameters
  P = data.frame(c("Rate"),c(l1),sqrt(Ji));
  colnames(P) = c("Aspect","Estimate","SE");
  ## Outcome
  # Estimate mean
  mu = 1/l1;
  dg = -1/(l1^2);
  se.mu = sqrt(dg*Ji*dg);
  # Estimate median
  me = (1/l1)*log(2);
  dg = -1/(l1^2)*log(2);
  se.me = sqrt(dg*Ji*dg);
  # Estimate variance
  v = 1/(l1^2);
  dg = -2/(l1^3);
  se.v = sqrt(dg*Ji*dg);
  # Outcome characteristics
  Y = data.frame(c("Mean","Median","Variance"),c(mu,me,v),c(se.mu,se.me,se.v));
  colnames(Y) = c("Aspect","Estimate","SE");
  # CIs
  z = qnorm(1-alpha/2);
  P$L = P$Estimate-z*P$SE;
  P$U = P$Estimate+z*P$SE;
  Y$L = Y$Estimate-z*Y$SE;
  Y$U = Y$Estimate+z*Y$SE;
  ## Format Results
  J = matrix(J);
  rownames(J) = colnames(J) = c("l");
  Out = new(Class="fit",Distribution="Exponential",Parameters=P,Information=J,Outcome=Y);
  return(Out);
}

########################
# Gamma Distribution
########################

#' Gamma Parameter Estimation
#'
#' Estimates parameters for gamma event times subject to non-informative
#' right censoring. The gamma distribution is parameterized in terms
#' of the shape \eqn{\alpha} and rate \eqn{\lambda}:
#' \deqn{f(t) = \frac{\lambda^{\alpha}}{\Gamma(\alpha)} t^{\alpha-1}e^{-\lambda t}, t>0}
#'
#' @param time Observation times.
#' @param status Status indicator, coded as 1 if an event was observed, 0 if censored.
#' @param alpha Significance level, for CIs.
#' @param eps Tolerance for Newton-Raphson iterations.
#' @param maxit Maximum number of NR iterations.
#' @importFrom expint gammainc
#' @importFrom methods new
#' @importFrom numDeriv grad
#' @importFrom stats qgamma var
#' @return An object of class \code{fit} containing the following:
#' \describe{
#'  \item{Parameters}{The estimated shape \eqn{\alpha} and rate \eqn{\lambda}.}
#'  \item{Information}{The observed information matrix.}
#'  \item{Outcome}{The fitted mean and variance.}
#' }

fit.Gamma = function(time,status,alpha=0.05,eps=1e-6,maxit=10){
  # Events
  n = length(time);
  nobs = sum(status);
  # Observed events
  tObs = time[status==1];
  tCen = time[status==0];
  # Flag presence of censored events
  flag = (length(tCen)>0);
  ## Function to calculate observed information
  D2aLogIncGamma = Vectorize(D2aLogIncGamma,vectorize.args="x");
  D2lLogIncGamma = Vectorize(D2lLogIncGamma,vectorize.args="u");
  D2mixLogIncGamma = Vectorize(D2mixLogIncGamma,vectorize.args="u");
  obsInfo = function(a,l){
    # Information for alpha
    Jaa = n*trigamma(a);
    # Information for lambda
    Jll = nobs*a/(l^2);
    # Cross information
    Jla = -(nobs/l);
    # Add corrections for censoring
    if(flag){
      Jaa = Jaa-sum(D2aLogIncGamma(a=a,x=l*tCen));
      Jll = Jll-sum(D2lLogIncGamma(a=a,l=l,u=tCen));
      Jla = Jla-sum(D2mixLogIncGamma(a=a,l=l,u=tCen));
    }
    # Output
    Out = matrix(c(Jaa,Jla,Jla,Jll),nrow=2);
    return(Out);
  }
  ## Function to calculate score
  Sa = sum(log(tObs));
  Sl = sum(tObs);
  D1aLogIncGamma = Vectorize(D1aLogIncGamma,vectorize.args="x");
  D1lLogIncGamma = Vectorize(D1lLogIncGamma,vectorize.args="u");
  Score = function(a,l){
    # Score for alpha
    Ua = nobs*log(l)+Sa-n*digamma(a)
    # Score for lambda
    Ul = nobs*a/l-Sl
    # Add corrections for censoring
    if(flag){
      Ua = Ua+sum(D1aLogIncGamma(a=a,x=l*tCen));;
      Ul = Ul+sum(D1lLogIncGamma(a=a,l=l,u=tCen));
    }
    # Output
    Out = c(Ua,Ul);
    return(Out);
  }
  ## Objective function
  Q = function(a,l){
    # Log likelihood
    Out = nobs*a*log(l)+a*Sa-l*Sl-n*log(gamma(a));
    # Add corrections for censoring
    if(flag){
      Out = Out+sum(gammainc(a,l*tCen));
    }
    return(Out);
  }
  ## NR Update
  Update = function(a,l){
    # Current score
    U = Score(a,l);
    # Inverse observed information
    Ji = fastInv(obsInfo(a,l));
    # Proposal
    Prop = c(a,l) + fastMMp(Ji,U);
    # Output
    Out = list("a"=Prop[1,1],"l"=Prop[2,1]);
    return(Out);
  }
  ## Initialize
  theta0 = list();
  theta0$a = mean(tObs)^2/var(tObs);
  theta0$l = mean(tObs)/var(tObs);
  theta0$ll = Q(a=theta0$a,l=theta0$l);
  ## Newton-Raphson
  for(i in 1:maxit){
    # Propose update
    theta1 = Update(a=theta0$a,l=theta0$l);
    # Proposed objective
    theta1$ll = Q(a=theta1$a,l=theta1$l);
    # Accept first update,
    # Otherwise, check for improvement
    if(i==1){
      theta0 = theta1;
    } else {
      delta = theta1$ll-theta0$ll;
      if(delta>eps){
        theta0 = theta1;
      } else {
        break;
      }
    }
  }; # End NR
  ## Report
  if(i<maxit){
    cat(paste0(i-1," update(s) performed before tolerance limit."),"\n");
  } else {
    cat(paste0(i," update(s) performed without reaching tolerance limit."));
  }
  # Maximum likelihood estimates
  a1 = theta0$a;
  l1 = theta0$l;
  ## Final Information
  J = obsInfo(a=a1,l=l1);
  Ji = fastInv(J);
  # Parameters
  P = data.frame(c("Shape","Rate"),c(a1,l1),sqrt(diag(Ji)));
  colnames(P) = c("Aspect","Estimate","SE");
  ## Outcome
  # Estimate mean
  mu = a1/l1;
  dg = c(1/l1,-1*a1/(l1^2));
  se.mu = sqrt(as.numeric(fastQF(X=dg,A=Ji)));
  # Estimate median
  me = qgamma(p=0.5,shape=a1,rate=l1);
  g = function(x){qgamma(p=0.5,shape=x[1],rate=x[2])};
  dg = grad(func=g,x=c(a1,l1));
  se.me = sqrt(as.numeric(fastQF(X=dg,A=Ji)));
  # Estimate variance
  v = a1/(l1^2);
  dg = c(1/(l1^2),-2*a1/(l1^3));
  se.v = sqrt(as.numeric(fastQF(X=dg,A=Ji)));
  # Outcome characteristics
  Y = data.frame(c("Mean","Median","Variance"),c(mu,me,v),c(se.mu,se.me,se.v));
  colnames(Y) = c("Aspect","Estimate","SE");
  # CIs
  z = qnorm(1-alpha/2);
  P$L = P$Estimate-z*P$SE;
  P$U = P$Estimate+z*P$SE;
  Y$L = Y$Estimate-z*Y$SE;
  Y$U = Y$Estimate+z*Y$SE;
  ## Format Results
  Out = new(Class="fit",Distribution="Gamma",Parameters=P,Information=J,Outcome=Y);
  return(Out);
}

########################
# Log-Logistic Distribution
########################

#' Log-Logistic Parameter Estimation
#'
#' Estimates parameters for log-logistic event times subject to non-informative
#' right censoring. The log-logistic distribution is parameterized in terms
#' of the shape \eqn{\alpha} and rate \eqn{\lambda}:
#' \deqn{f(t) = \frac{\alpha\lambda(\lambda t)^{\alpha-1}}{[1+(\lambda t)^{\alpha}]^{2}}, t>0}
#' 
#' For the log-logistic distribution, the mean is only defined if the shape
#' parameter \eqn{\alpha>1}, and the variance if the shape parameter \eqn{\alpha>2}. 
#' Consequently, estimates of the outcome mean and variance are only returned
#' if the estimated shape parameter exceeds these thresholds. For \eqn{\alpha\ll 1}, 
#' the fitting function may fail. 
#'
#' @param time Observation times.
#' @param status Status indicator, coded as 1 if an event was observed, 0 if censored.
#' @param alpha Significance level, for CIs.
#' @param eps Tolerance for Newton-Raphson iterations.
#' @param maxit Maximum number of NR iterations.
#' @importFrom methods new
#' @importFrom stats median
#' @return An object of class \code{fit} containing the following:
#' \describe{
#'  \item{Parameters}{The estimated location \eqn{\mu} and scale \eqn{\sigma}.}
#'  \item{Information}{The observed information matrix.}
#'  \item{Outcome}{The fitted mean, median, and variance.}
#' }

fit.LogLogistic = function(time,status,alpha=0.05,eps=1e-6,maxit=10){
  # Events
  n = length(time);
  nobs = sum(status);
  # Observed events
  tObs = time[status==1];
  tCen = time[status==0];
  ## Function to calculate observed information
  obsInfo = function(a,l){
    b = (l*time)^a;
    # Observed information for alpha
    Jaa = nobs/(a^2)+sum((1+status)*b*(log(l*time))^2/((1+b)^2));
    # Observed information for lambda
    Jll = nobs*a/(l^2)+sum((1+status)*a*(b/(l^2))*((a-1)-b)/((1+b)^2));
    # Cross information
    Jal = -1*nobs/l+sum((1+status)*(b/l)*(1+b+a*log(l*time))/((1+b)^2));
    # Output
    Out = matrix(c(Jaa,Jal,Jal,Jll),nrow=2);
    return(Out);
  }
  ## Function to calculate score
  Score = function(a,l){
    b = (l*time)^a;
    # For alpha
    Ua = nobs/a+nobs*log(l)+sum(log(tObs))-sum((1+status)*b*log(l*time)/(1+b));
    # For lambda
    Ul = nobs*a/l-sum((1+status)*a*(b/l)/(1+b));
    Out = c(Ua,Ul);
    return(Out);
  }
  ## Objective function
  Q = function(a,l){
    Out = nobs*log(a)+nobs*a*log(l)+a*sum(log(tObs))-sum((1+status)*log(1+(l*time)^a));
    return(Out);
  }
  ## NR Update
  Update = function(a,l){
    # Current score
    U = Score(a,l);
    # Inverse observed information
    Ji = fastInv(obsInfo(a,l));
    # Proposal
    Prop = c(a,l) + fastMMp(Ji,U);
    # Output
    Out = list("a"=Prop[1,1],"l"=Prop[2,1]);
    return(Out);
  }
  ## Initialize
  theta0 = list();
  theta0$l = 1/median(tObs);
  theta0$a = (log(1.5)/2)*as.numeric(1/log(quantile(tObs,probs=0.6)*theta0$l)-1/log(quantile(tObs,probs=0.4)*theta0$l));
  theta0$ll = Q(a=theta0$a,l=theta0$l);
  ## Newton-Raphson
  for(i in 1:maxit){
    # Propose update
    theta1 = Update(a=theta0$a,l=theta0$l);
    # Proposed objective
    theta1$ll = Q(a=theta1$a,l=theta1$l);
    # Accept first update,
    # Otherwise, check for improvement
    if(i==1){
      theta0 = theta1;
    } else {
      delta = theta1$ll-theta0$ll;
      if(delta>eps){
        theta0 = theta1;
      } else {
        break;
      }
    }
  }; # End NR
  ## Report
  if(i<maxit){
    cat(paste0(i-1," update(s) performed before tolerance limit."),"\n");
  } else {
    cat(paste0(i," update(s) performed without reaching tolerance limit."));
  }
  # Maximum likelihood estimates
  a1 = theta0$a;
  l1 = theta0$l;
  ## Final Information
  J = obsInfo(a=a1,l=l1);
  Ji = fastInv(J);
  # Parameters
  P = data.frame(c("Shape","Rate"),c(a1,l1),sqrt(diag(Ji)));
  colnames(P) = c("Aspect","Estimate","SE");
  ## Outcome
  # Estimate mean
  if(a1>1){
    mu = (pi/(a1*l1))/sin(pi/a1);
    dg = -mu*c((sin(pi/a1)-(pi/a1)*cos(pi/a1))/(a1*sin(pi/a1)),(1/l1));
    se.mu = sqrt(as.numeric(fastQF(X=dg,A=Ji)));
  } else {
    mu = NA;
    se.mu = NA;
  }
  # Estimate median
  me = (1/l1);
  dg = c(0,-1/(l1^2));
  se.me = sqrt(as.numeric(fastQF(X=dg,A=Ji)));
  # Estimate variance
  if(a1>2){
    v = (1/l1)^2*(2*(pi/a1)/sin(2*pi/a1)-(pi/a1)^2/(sin(pi/a1))^2);
    da1 = -2*pi/(a1^2)/(sin(2*pi/a1)^2)*(sin(2*pi/a1)-(2*pi/a1)*cos(2*pi/a1));
    da2 = -2*(pi^2)/(a1^3)/(sin(pi/a1)^4)*(sin(pi/a1)^2-pi/(2*a1)*sin(2*pi/a1));
    dg = c((1/(l1^2))*(da1-da2),-2/(l1)*v);
    se.v = sqrt(as.numeric(fastQF(X=dg,A=Ji)));
  } else {
    v = NA;
    se.v = NA;
  }
  # Outcome characteristics
  Y = data.frame(c("Mean","Median","Variance"),c(mu,me,v),c(se.mu,se.me,se.v));
  colnames(Y) = c("Aspect","Estimate","SE");
  # CIs
  z = qnorm(1-alpha/2);
  P$L = P$Estimate-z*P$SE;
  P$U = P$Estimate+z*P$SE;
  Y$L = Y$Estimate-z*Y$SE;
  Y$U = Y$Estimate+z*Y$SE;
  ## Format Results
  Out = new(Class="fit",Distribution="Log-Logistic",Parameters=P,Information=J,Outcome=Y);
  return(Out);
}

########################
# Log-Normal Distribution
########################

#' Log-Normal Parameter Estimation
#'
#' Estimates parameters for log-normal event times subject to non-informative
#' right censoring. The log-normal distribution is parameterized in terms
#' of the location \eqn{\mu} and scale \eqn{\sigma}:
#' \deqn{f(t) = \phi\left(\frac{\ln t-\mu}{\sigma}\right)\frac{1}{t\sigma}, t>0}
#' 
#' In the log-normal density, \eqn{\phi(\cdot)} denotes the standard normal pdf. 
#'
#' @param time Observation times.
#' @param status Status indicator, coded as 1 if an event was observed, 0 if censored.
#' @param alpha Significance level, for CIs.
#' @param eps Tolerance for Newton-Raphson iterations.
#' @param maxit Maximum number of NR iterations.
#' @importFrom methods new
#' @importFrom stats var
#' @return An object of class \code{fit} containing the following:
#' \describe{
#'  \item{Parameters}{The estimated location \eqn{\mu} and scale \eqn{\sigma}.}
#'  \item{Information}{The observed information matrix.}
#'  \item{Outcome}{The fitted mean, median, and variance.}
#' }

fit.LogNormal = function(time,status,alpha=0.05,eps=1e-6,maxit=10){
  # Events
  n = length(time);
  nobs = sum(status);
  # Observed events
  tObs = time[status==1];
  tCen = time[status==0];
  ## Function to calculate observed information
  obsInfo = function(m,s){
    # Standardize
    zObs = (log(tObs)-m)/s;
    zCen = (log(tCen)-m)/s;
    # Information for mu
    Jmm = nobs/(s^2)-1/(s^2)*sum(D2LogNormSurv(zCen));
    # Information for sigma
    Jss = -nobs/(s^2)+3/(s^2)*sum(zObs^2)-2/(s^2)*sum(D1LogNormSurv(zCen)*zCen)-1/(s^2)*sum(D2LogNormSurv(zCen)*(zCen)^2);
    # Cross information
    Jms = 2/(s^2)*sum(zObs)-1/(s^2)*sum(D1LogNormSurv(zCen))-1/(s^2)*sum(D2LogNormSurv(zCen)*zCen);
    # Output
    Out = matrix(c(Jmm,Jms,Jms,Jss),nrow=2);
    return(Out);
  }
  ## Function to calculate score
  Score = function(m,s){
    # Standardize
    zObs = (log(tObs)-m)/s;
    zCen = (log(tCen)-m)/s;
    # Score for mu
    Um = (1/s)*sum(zObs)-(1/s)*sum(D1LogNormSurv(zCen));
    # Score for sigma
    Us = -1*(nobs/s)+(1/s)*sum(zObs^2)-(1/s)*sum(D1LogNormSurv(zCen)*zCen);
    # Output
    Out = c(Um,Us);
    return(Out);
  }
  ## Objective function
  Q = function(m,s){
    zObs = (log(tObs)-m)/s;
    zCen = (log(tCen)-m)/s;
    Out = -nobs*log(s)-(1/2)*sum(zObs^2)+sum(pnorm(q=zCen,lower.tail=F,log.p=T));
    return(Out);
  }
  ## NR Update
  Update = function(m,s){
    # Current score
    U = Score(m,s);
    # Inverse observed information
    Ji = fastInv(obsInfo(m,s));
    # Proposal
    Prop = c(m,s) + fastMMp(Ji,U);
    # Output
    Out = list("m"=Prop[1,1],"s"=Prop[2,1]);
    return(Out);
  }
  ## Initialize
  theta0 = list();
  theta0$m = mean(log(tObs));
  theta0$s = sqrt(var(log(tObs)));
  theta0$ll = Q(m=theta0$m,s=theta0$s);
  ## Newton-Raphson
  for(i in 1:maxit){
    # Propose update
    theta1 = Update(m=theta0$m,s=theta0$s);
    # Proposed objective
    theta1$ll = Q(m=theta1$m,s=theta1$s);
    # Accept first update,
    # Otherwise, check for improvement
    if(i==1){
      theta0 = theta1;
    } else {
      delta = theta1$ll-theta0$ll;
      if(delta>eps){
        theta0 = theta1;
      } else {
        break;
      }
    }
  }; # End NR
  ## Report
  if(i<maxit){
    cat(paste0(i-1," update(s) performed before tolerance limit."),"\n");
  } else {
    cat(paste0(i," update(s) performed without reaching tolerance limit."));
  }
  # Maximum likelihood estimates
  m1 = theta0$m;
  s1 = theta0$s;
  ## Final Information
  J = obsInfo(m=m1,s=s1);
  Ji = fastInv(J);
  # Parameters
  P = data.frame(c("Location","Scale"),c(m1,s1),sqrt(diag(Ji)));
  colnames(P) = c("Aspect","Estimate","SE");
  ## Outcome
  # Estimate mean
  mu = exp(m1+s1^2/2);
  dg = mu*c(1,s1);
  se.mu = sqrt(as.numeric(fastQF(X=dg,A=Ji)));
  # Estimate median
  me = exp(m1);
  dg = c(mu,0);
  se.me = sqrt(as.numeric(fastQF(X=dg,A=Ji)));
  # Estimate variance
  v = (exp(s1^2)-1)*exp(2*m1+s1^2);
  dg = 2*exp(2*m1+s1^2)*c(exp(s1^2)-1,s1*(2*exp(s1^2)-1));
  se.v = sqrt(as.numeric(fastQF(X=dg,A=Ji)));
  # Outcome characteristics
  Y = data.frame(c("Mean","Median","Variance"),c(mu,me,v),c(se.mu,se.me,se.v));
  colnames(Y) = c("Aspect","Estimate","SE");
  # CIs
  z = qnorm(1-alpha/2);
  P$L = P$Estimate-z*P$SE;
  P$U = P$Estimate+z*P$SE;
  Y$L = Y$Estimate-z*Y$SE;
  Y$U = Y$Estimate+z*Y$SE;
  ## Format Results
  Out = new(Class="fit",Distribution="Log-Normal",Parameters=P,Information=J,Outcome=Y);
  return(Out);
}

########################
# Weibull Distribution
########################

#' Weibull Parameter Estimation
#'
#' Estimates parameters for Weibull event times subject to non-informative
#' right censoring. The Weibull distribution is parameterized in terms
#' of the shape \eqn{\alpha} and rate \eqn{\lambda}:
#' \deqn{f(t) = \alpha\lambda^{\alpha}t^{\alpha-1}e^{-(\lambda t)^{\alpha}}, t>0}
#'
#' @param time Observation times.
#' @param status Status indicator, coded as 1 if an event was observed, 0 if censored.
#' @param alpha Significance level, for CIs.
#' @importFrom methods new
#' @importFrom stats quantile uniroot
#' @return An object of class \code{fit} containing the following:
#' \describe{
#'  \item{Parameters}{The estimated shape \eqn{\alpha} and rate \eqn{\lambda}.}
#'  \item{Information}{The observed information matrix.}
#'  \item{Outcome}{The fitted mean, median, and variance.}
#' }

fit.Weibull = function(time,status,alpha=0.05){
  # Events
  n = length(time);
  nobs = sum(status);
  # Observed events
  tObs = time[status==1];
  # Profile score for shape
  pscore = function(a){
    # Components
    v = time^a;
    logTime = log(time);
    # Score
    Out = nobs/a-sum(v*logTime)/sum(v)*nobs+sum(status*logTime);
    return(Out);
  }
  # Initialize
  l0 = 1/quantile(x=tObs,probs=c(1-exp(-1)));
  a0 = digamma(1)/(log(l0)+mean(log(tObs)));
  # Search interval
  L = 0;
  sL = sign(pscore(0));
  # Note: Sign at L=0 is always positive
  U = a0;
  sU = sign(pscore(U));
  if(sL==sU){
    # Search for an upper bound
    while(sU==sL){
      U = 2*U;
      sU = sign(pscore(U));
    }
  }
  # Find MLE of alpha
  a1 = uniroot(f=pscore,lower=L,upper=U)$root;
  # MLE of lambda
  l1 = (sum(time^a1)/nobs)^(-1/a1);
  ## Information
  S1 = sum(time^(a1));
  S2 = sum(time^(a1)*log(time));
  S3 = sum(time^(a1)*(log(time))^2);
  # Information components
  Jaa = (nobs/(a1^2))+(l1^a1)*(log(l1))^2*S1+2*(l1^a1)*log(l1)*S2+(l1^a1)*S3;
  Jll = (nobs*a1)/(l1^2)+a1*(a1-1)*(l1^(a1-2))*S1;
  Jla = -(nobs/l1)+(l1^(a1-1))*S1+a1*(l1^(a1-1))*log(l1)*S1+a1*(l1^(a1-1))*S2;
  # Information matrix
  J = matrix(c(Jaa,Jla,Jla,Jll),nrow=2);
  colnames(J) = rownames(J) = c("a","l");
  Ji = fastInv(J);
  # Parameters
  P = data.frame(c("Shape","Rate"),c(a1,l1),sqrt(diag(Ji)));
  colnames(P) = c("Aspect","Estimate","SE");
  ## Outcome
  # Estimate mean
  mu = (1/l1)*gamma(1+1/a1);
  dg = -gamma(1+1/a1)*c(digamma(1+1/a1)/a1^2,1/l1^2);
  se.mu = sqrt(as.numeric(fastQF(dg,Ji)));
  # Estimate median
  me = (1/l1)*(log(2))^(1/a1);
  dg = -(1/l1)*(log(2))^(1/a1)*c(log(log(2))/a1^2,(1/l1));
  se.me = sqrt(as.numeric(fastQF(dg,Ji)));
  # Estimate variance
  v = (1/l1^2)*(gamma(1+2/a1)-gamma(1+1/a1)^2);
  dg = -(2/l1^2)*c((1/a1^2)*(gamma(1+2/a1)*digamma(1+2/a1)-(gamma(1+1/a1)^2)*digamma(1+1/a1)),
                   (1/l1)*(gamma(1+2/a1)-gamma(1+1/a1)^2));
  se.v = sqrt(as.numeric(fastQF(dg,Ji)));
  # Outcome characteristics
  Y = data.frame(c("Mean","Median","Variance"),c(mu,me,v),c(se.mu,se.me,se.v));
  colnames(Y) = c("Aspect","Estimate","SE");
  # CIs
  z = qnorm(1-alpha/2);
  P$L = P$Estimate-z*P$SE;
  P$U = P$Estimate+z*P$SE;
  Y$L = Y$Estimate-z*Y$SE;
  Y$U = Y$Estimate+z*Y$SE;
  ## Format Results
  Out = new(Class="fit",Distribution="Weibull",Parameters=P,Information=J,Outcome=Y);
  return(Out);
}
