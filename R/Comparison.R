# Purpose: Two group comparison
# Updated: 180817

########################
# Auxiliary Functions
########################

#' Difference of Estimates
#'
#' Calculate CIs and p-value for the difference of estimated parameters
#'
#' @param t1 Treatment estimate
#' @param s1 Treatment standard error
#' @param t0 Reference estimate
#' @param s0 Reference standard error
#' @param sig Significance level
#'
#' @return Data.frame containing estimated difference, its standard error, lower and
#'   upper confidence bounds, and a p-value assessing the null hypothesis of no
#'   difference.
#'
#' @importFrom stats pnorm qnorm

estDiff = function(t1,s1,t0,s0,sig=0.05){
  # Critical value
  z = qnorm(1-sig/2);
  # Point
  d1 = t1-t0;
  # SE
  d1.se = sqrt(s1^2+s0^2);
  # CI
  d1.L = d1-z*d1.se;
  d1.U = d1+z*d1.se;
  # P
  p = 2*pnorm(abs(d1/d1.se),lower.tail=F);
  # Output
  Out = data.frame("Point"=d1,"SE"=d1.se,"L"=d1.L,"U"=d1.U,"p"=p);
  return(Out);
}

#' Ratio of Estimates
#'
#' Calculate CIs and p-value for the ratio of estimated parameters
#'
#' @param t1 Treatment estimate
#' @param s1 Treatment standard error
#' @param t0 Reference estimate
#' @param s0 Reference standard error
#' @param sig Significance level
#'
#' @return Data.frame containing estimated ratio, its standard error, lower and
#'   upper confidence bounds, and a p-value assessing the null hypothesis that
#'   the ratio is unity.
#'
#' @importFrom stats pnorm qnorm

estRatio = function(t1,s1,t0,s0,sig=0.05){
  # Critical value
  z = qnorm(1-sig/2);
  # Point
  r1 = log(t1)-log(t0);
  # SE, log scale
  r1.se = sqrt(s1^2/(t1^2)+s0^2/(t0^2));
  # SE, ratio scale
  r2.se = sqrt(s1^2/(t0^2)+(t1^2*s0^2)/(t0^4));
  # CI
  r1.L = r1-z*r1.se;
  r1.U = r1+z*r1.se;
  # P
  p = 2*pnorm(abs(r1/r1.se),lower.tail=F);
  # Output
  Out = data.frame("Point"=exp(r1),"SE"=r2.se,"L"=exp(r1.L),"U"=exp(r1.U),"p"=p);
  return(Out);
}

########################
# Distribution Comparison
########################

#' Compare Parametric Survival Distribution
#'
#' Compares the means and medians of parametric survival distributions fit to
#' each of two treatment arms. Available distributions include: exponential,
#' gamma, generalized gamma, log-logistic, log-normal, and Weibull.
#'
#' @param time Observation time.
#' @param status Status indicator, coded as 1 if an event was observed, 0 if
#'   censored.
#' @param arm Treatment indicator, coded as 1 for the target group, 0 for the
#'   reference group.
#' @param dist1 Distribution to fit for the target group. Selected from among:
#'   exp, gamma, gengamma, log-logistic, log-normal, and weibull.
#' @param dist0 Distribution to fit for the reference group. Same choices as for
#'   the target group. If omitted, defaults to the distribution specified for
#'   the target group.
#' @param sig Significance level, for CIs.
#' @param init1 List of initial parameter values for the target group. See
#'   individual distributions for naming convention.
#' @param init0 List of initial parameter values for the reference group. See
#'   individual distributions for naming convention.
#' @param eps Tolerance for Newton-Raphson iterations.
#' @param maxit Maximum number of NR iterations.
#' @param report Report fitting progress?
#'
#' @importFrom methods new
#'
#' @export
#' @return An object of class \code{contrast} containing the following:
#' \describe{
#'   \item{Model1}{The fitted model for the target group.}
#'   \item{Model0}{The fitted model for the reference group.}
#'   \item{Contrast}{Contrasts of means and medians.}
#' }
#'
#' @seealso
#' \itemize{
#'   \item{}{Fitting function for parametric survival distributions \code{\link{fitParaSurv}}}
#' }
#'
#' @examples
#' set.seed(100);
#' ## Weibull and Weibull, different means and medians
#' # Target group
#' D1 = rWeibull(n=1e3,a=1,l=1,p=0.2);
#' D1$arm = 1;
#' # Reference group
#' D0 = rWeibull(n=1e3,a=1,l=2,p=0.2);
#' D0$arm = 0;
#' # Overall data set
#' D = rbind(D1,D0);
#' # Comparison
#' E = compParaSurv(time=D$time,status=D$status,arm=D$arm,dist1="weibull");
#'
#' ## Gamma and Weibull, different means and medians
#' # Target group
#' D1 = rGamma(n=1e3,a=2,l=2,p=0.2);
#' D1$arm = 1;
#' # Reference group
#' D0 = rWeibull(n=1e3,a=2,l=2/sqrt(pi),p=0.2);
#' D0$arm = 0;
#' # Overall data set
#' D = rbind(D1,D0);
#' # Comparison
#' E = compParaSurv(time=D$time,status=D$status,arm=D$arm,dist1="gamma",dist0="weibull");
#'
#' ## Weibull and Log-normal, different means, same median
#' # Target group
#' D1 = rLogNormal(n=1e3,m=0,s=2,p=0.2);
#' D1$arm = 1;
#' # Reference group
#' D0 = rWeibull(n=1e3,a=2,l=sqrt(log(2)),p=0.2);
#' D0$arm = 0;
#' # Overall data set
#' D = rbind(D1,D0);
#' # Comparison
#' E = compParaSurv(time=D$time,status=D$status,arm=D$arm,dist1="log-normal",dist0="weibull");

compParaSurv = function(time,status,arm,dist1="weibull",dist0,
                        init1=NULL,init0=NULL,sig=0.05,eps=1e-6,maxit=10,report=F){
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
  # Distributions
  choices = c("exp","gamma","gengamma","log-logistic","log-normal","weibull");
  if(!(dist1 %in% choices)){stop(c("Select distribution from among:\n",paste0(choices,collapse=" ")))};
  if(missing(dist0)){dist0=dist1};
  if(!(dist0 %in% choices)){stop(c("Select distribution from among:\n",paste0(choices,collapse=" ")))};

  ## Partition data
  D0 = data.frame("time"=time[arm==0],"status"=status[arm==0]);
  D1 = data.frame("time"=time[arm==1],"status"=status[arm==1]);
  # Fit to reference arm
  M0 = fitParaSurv(time=D0$time,status=D0$status,dist=dist0,sig=sig,init=init0,eps=eps,maxit=maxit,report=report);
  m0 = M0@Outcome;
  # Fit to treatment arm
  M1 = fitParaSurv(time=D1$time,status=D1$status,dist=dist1,sig=sig,init=init1,eps=eps,maxit=maxit,report=report);
  m1 = M1@Outcome;

  # Difference of means
  Contrast.mu = estDiff(t1=m1[1,2],s1=m1[1,3],t0=m0[1,2],s0=m0[1,3],sig=sig);
  # Ratio of means
  Contrast.rmu = estRatio(t1=m1[1,2],s1=m1[1,3],t0=m0[1,2],s0=m0[1,3],sig=sig);
  # Difference of medians
  Contrast.me = estDiff(t1=m1[2,2],s1=m1[2,3],t0=m0[2,2],s0=m0[2,3],sig=sig);
  # Ratio of medians
  Contrast.rme = estRatio(t1=m1[2,2],s1=m1[2,3],t0=m0[2,2],s0=m0[2,3],sig=sig);

  # Contrast frame
  Contrast = rbind(Contrast.mu,Contrast.me,Contrast.rmu,Contrast.rme);
  Contrast = cbind("Contrast"=c("Mean1-Mean0","Med1-Med0","Mean1/Mean0","Med1/Med0"),Contrast);

  # Output
  Out = new(Class="contrast",Dist1=M1@Distribution,Dist0=M0@Distribution,Model1=M1,Model0=M0,Contrast=Contrast);
  return(Out);
}
