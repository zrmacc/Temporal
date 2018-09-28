# Purpose: Parameter estimation
# Updated: 180817

########################
# Master Fitting Function
########################

#' Fit Parametric Survival Distribution
#'
#' Estimates parametric survival distributions using event times subject to
#' non-informative right censoring. Available distributions include:
#' exponential, gamma, generalized gamma, log-logistic, log-normal, and Weibull.
#'
#' @param time Observation times.
#' @param status Status indicator, coded as 1 if an event was observed, 0 if
#'   censored.
#' @param dist Distribution to fit, selected from among: exp, gamma, gengamma,
#'   log-logistic, log-normal, and weibull.
#' @param sig Significance level, for CIs.
#' @param init List of initial parameter values. See individual distributions for
#'   naming convention.
#' @param eps Tolerance for Newton-Raphson iterations.
#' @param maxit Maximum number of NR iterations.
#' @param report Report fitting progress?
#'
#' @export
#'
#' @return An object of class \code{fit} containing the following:
#' \describe{
#'   \item{Parameters}{The estimated shape and rate parameters.}
#'   \item{Information}{The observed information matrix.}
#'   \item{Outcome}{The fitted mean, median, and variance.}
#' }
#'
#' @seealso
#' \itemize{
#'   \item{}{Between group comparison of survival experience \code{\link{compParaSurv}}}
#'   \item{}{Exponential distribution \code{\link{fit.Exp}}}
#'   \item{}{Gamma distribution \code{\link{fit.Gamma}}}
#'   \item{}{Generalized gamma distribution \code{\link{fit.GenGamma}}}
#'   \item{}{Log-logistic distribution \code{\link{fit.LogLogistic}}}
#'   \item{}{Log-normal distribution \code{\link{fit.LogNormal}}}
#'   \item{}{Weibull distribution \code{\link{fit.Weibull}}}
#' }
#'
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

fitParaSurv = function(time,status,dist="weibull",sig=0.05,init=NULL,eps=1e-6,maxit=10,report=F){
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
  choices = c("exp","gamma","gengamma","log-logistic","log-normal","weibull");
  if(!(dist %in% choices)){stop(c("Select distribution from among:\n",paste0(choices,collapse=" ")))};
  # Initialization
  if(!is.null(init)&!is.list(init)){stop("init is expected as a list.")};

  # Exponential
  if(dist=="exp"){
    M = fit.Exp(time=time,status=status,sig=sig);
  }
  # Gamma
  if(dist=="gamma"){
    M = fit.Gamma(time=time,status=status,sig=sig,init=init,eps=eps,maxit=maxit,report=report);
  }
  # Generalized gamma
  if(dist=="gengamma"){
    M = fit.GenGamma(time=time,status=status,sig=sig,init=init,eps=eps,maxit=maxit,report=report);
  }
  # Log logistic
  if(dist=="log-logistic"){
    M = fit.LogLogistic(time=time,status=status,sig=sig,init=init,eps=eps,maxit=maxit,report=report);
  }
  # Log normal
  if(dist=="log-normal"){
    M = fit.LogNormal(time=time,status=status,sig=sig,init=init,eps=eps,maxit=maxit,report=report);
  }
  # Weibull
  if(dist=="weibull"){
    M = fit.Weibull(time=time,status=status,sig=sig,init=init);
  }
  # Output
  return(M);
}
