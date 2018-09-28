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
#' @param sig Significance level, for CIs.
#' @importFrom methods new
#' @return An object of class \code{fit} containing the following:
#' \describe{
#'  \item{Parameters}{The estimated model parameters.}
#'  \item{Information}{The observed information matrix.}
#'  \item{Outcome}{The fitted mean, median, and variance of the time to event distribution.}
#' }

fit.Exp = function(time,status,sig=0.05){
  # Events
  n = length(time);
  nobs = sum(status);
  # MLE of lambda
  l1 = nobs/sum(time);
  # Information
  J = nobs/(l1^2);
  Ji = 1/J;

  # Parameter frame
  P = data.frame(c("Rate"),c(l1),sqrt(Ji));
  colnames(P) = c("Aspect","Estimate","SE");

  ## Fitted Distribution
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
  z = qnorm(1-sig/2);
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
