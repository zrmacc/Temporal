# Purpose: Estimation of exponential distribution
# Updated: 20/03/07

########################
# Exponential Distribution
########################

#' Exponential Distribution Parameter Estimation
#'
#' Estimates parameters for exponential event times subject to non-informative
#' right censoring. The exponential distribution is parameterized in terms
#' of the rate \eqn{\lambda}: \deqn{f(t) = \lambda e^{-\lambda t}, t>0}
#'
#' @param time Numeric observation times.
#' @param status Status indicator, coded as 1 if an event was observed, 0 if censored.
#' @param sig Significance level, for CIs.
#' @param tau Optional truncation times for calculating RMSTs. 
#'
#' @return An object of class \code{fit} containing the following:
#' \describe{
#'  \item{Parameters}{The estimated model parameters.}
#'  \item{Information}{The observed information matrix.}
#'  \item{Outcome}{The fitted mean, median, and variance of the time to event distribution.}
#'  \item{RMST}{The estimated RMSTs, if tau was specified.}
#' }
#'
#' @importFrom methods new
#'
#' @seealso
#' \itemize{
#'   \item{Fitting function for parametric survival distributions \code{\link{fitParaSurv}}}
#' }
#'
#' @examples
#' # Generate exponential data with 20% censoring
#' data = genData(n=1e3,dist="exp",theta=c(2),p=0.2);
#' # Estimate
#' fit = fitParaSurv(time=data$time,status=data$status,dist="exp");

fit.Exp = function(time,status,sig=0.05,tau=NULL){
  # Events
  n = length(time);
  nobs = sum(status);
  # MLE of lambda
  l = nobs/sum(time);
  # Information
  J = nobs/(l^2);
  Ji = 1/J;

  # Parameter frame
  P = data.frame(c("Rate"),c(l),sqrt(Ji));
  colnames(P) = c("Aspect","Estimate","SE");

  ## Fitted Distribution
  # Estimate mean
  mu = 1/l;
  dg = -1/(l^2);
  se.mu = sqrt(dg*Ji*dg);

  # Estimate median
  me = (1/l)*log(2);
  dg = -1/(l^2)*log(2);
  se.me = sqrt(dg*Ji*dg);

  # Estimate variance
  v = 1/(l^2);
  dg = -2/(l^3);
  se.v = sqrt(dg*Ji*dg);
  
  # Outcome characteristics
  Y = data.frame(c("Mean","Median","Variance"),
                 c(mu,me,v),c(se.mu,se.me,se.v),
                 stringsAsFactors=FALSE);
  colnames(Y) = c("Aspect","Estimate","SE");

  # CIs
  z = qnorm(1-sig/2);
  P$L = P$Estimate-z*P$SE;
  P$U = P$Estimate+z*P$SE;
  Y$L = Y$Estimate-z*Y$SE;
  Y$U = Y$Estimate+z*Y$SE;

  # Fitted survival function
  S = function(t){exp(-l*t)};

  # Format Results
  J = matrix(J);
  rownames(J) = colnames(J) = c("l");
  Out = new(Class="fit",Distribution="exp",Parameters=P,Information=J,Outcome=Y,S=S);
  
  ## Add RMST if requested. 
  if(is.numeric(tau)){
    rmst = paraRMST(fit=Out,sig=sig,tau=tau);
    Out@RMST = rmst;
  }
  
  return(Out);
}
