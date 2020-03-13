# Purpose: Parameter estimation
# Updated: 20/03/07

########################
# Master Fitting Function
########################

#' Fit Parametric Survival Distribution
#'
#' Estimates parametric survival distributions using event times subject to
#' non-informative right censoring. Available distributions include:
#' exponential, gamma, generalized gamma, log-normal, and Weibull.
#'
#' @param time Numeric observation times.
#' @param status Status indicator, coded as 1 if observed, 0 if censored.
#' @param dist String, distribution to fit, selected from among: exp, gamma, gen-gamma
#'   log-normal, and weibull.
#' @param tau Optional truncation time for calculating RMSTs.
#' @param sig Significance level, for CIs.
#' @param init Numeric vector of initial parameters. See individual distributions for
#'   parameter order.
#' @param bL If dist="gen-gamma", lower limit on possible values for beta.
#' @param bU If dist="gen-gamma", upper limit on possible values for beta. 
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
#'  \item{RMST}{The estimated RMSTs, if tau was specified.}
#' }
#'
#' @seealso
#' \itemize{
#'   \item{Between group comparison of survival experience \code{\link{compParaSurv}}}
#'   \item{Exponential distribution \code{\link{fit.Exp}}}
#'   \item{Gamma distribution \code{\link{fit.Gamma}}}
#'   \item{Generalized gamma distribution \code{\link{fit.GenGamma}}}
#'   \item{Log-normal distribution \code{\link{fit.LogNormal}}}
#'   \item{Weibull distribution \code{\link{fit.Weibull}}}
#' }
#'
#' @examples
#' # Generate Gamma data with 20% censoring
#' D = genData(n=1e3,dist="gamma",theta=c(2,2),p=0.2);
#' # Fit gamma distribution
#' M = fitParaSurv(time=D$time,status=D$status,dist="gamma");
#'
#' # Generate Weibull data with 10% censoring
#' D = genData(n=1e3,dist="weibull",theta=c(2,2),p=0.1);
#' # Fit weibull distribution, calculate RMST at tau=0.5
#' M = fitParaSurv(time=D$time,status=D$status,dist="weibull",tau=0.5);

fitParaSurv = function(time,status=NULL,dist="weibull",tau=NULL,sig=0.05,
                       init=NULL,bL=0.1,bU=10,eps=1e-6,maxit=10,report=F){
  
  # Input checks
  n = length(time);
  
  ## Positivity
  if(min(time)<0){
    stop("Strictly positive observation times required.");
  }
 
  ## Status
  if(is.null(status)){
    status=rep(1,n);
    warning("Since status was not supplied, all events are assumed observed.");
  }
  pass = checkStatus(n=n,status=status);
  if(!pass){
    stop("Status check failed.");
  }
  
  ## Distribution
  pass = checkDist(dist=dist);
  if(!pass){
    stop("Distribution check failed.");
  }
  
  ## Initialization
  pass = checkInit(dist=dist,init=init);
  if(!pass){
    stop("Initialization check failed.");
  }

  # Model fitting
  if(dist=="exp"){
    fit = fit.Exp(time=time,status=status,sig=sig,tau=tau);
  } else if(dist=="gamma"){
    fit = fit.Gamma(time=time,status=status,sig=sig,tau=tau,init=init,eps=eps,maxit=maxit,report=report);
  } else if(dist=="gen-gamma"){
    fit = fit.GenGamma(time=time,status=status,sig=sig,tau=tau,init=init,eps=eps,maxit=maxit,report=report);
  } else if(dist=="log-normal"){
    fit = fit.LogNormal(time=time,status=status,sig=sig,tau=tau,init=init,eps=eps,maxit=maxit,report=report);
  } else if(dist=="weibull"){
    fit = fit.Weibull(time=time,status=status,sig=sig,tau=tau,init=init);
  }
  
  # Output
  return(fit);
}
