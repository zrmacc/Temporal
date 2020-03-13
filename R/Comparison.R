# Purpose: Two group comparison
# Updated: 20/03/11

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
  Out = data.frame("Point"=d1,"SE"=d1.se,"L"=d1.L,"U"=d1.U,"P"=p);
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
  Out = data.frame("Point"=exp(r1),"SE"=r2.se,"L"=exp(r1.L),"U"=exp(r1.U),"P"=p);
  return(Out);
}

########################
# Distribution Comparison
########################

#' Compare Parametric Survival Distribution
#'
#' Compares the means and medians of parametric survival distributions fit to
#' each of two treatment arms. Available distributions include: exponential,
#' gamma, generalized gamma, log-normal, and Weibull.
#' 
#' Status is encoded 0 for censored and 1 for observed. Arm is encoded 0 for 
#' reference, 1 for target. Tau is an optional numeric vector of truncation times 
#' for calculating restricted mean survival time, which is the area under the
#' survival curve up to the specified truncation point. 
#' 
#' 
#'
#' @param time Observation time.
#' @param status Status indicator, coded as 1 if an event was observed, 0 if
#'   censored.
#' @param arm Treatment indicator, coded as 1 for the target group, 0 for the
#'   reference group.
#' @param tau Optional truncation times for calculating RMST. 
#' @param dist1 Distribution to fit for the target group. Selected from among:
#'   exp, gamma, gengamma, log-normal, and weibull.
#' @param dist0 Distribution to fit for the reference group. Same choices as for
#'   the target group. If omitted, defaults to the distribution specified for
#'   the target group.
#' @param boot If provided, integer number of bootstraps for constructing CIs.
#' @param perm If provided, integer number of permutations for calculating p-values. 
#' @param sig Significance level, for CIs.
#' @param init1 Initial parameter values for the target group.
#' @param init0 Initial parameter values for the reference group. 
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
#'   \item{Fitting function for parametric survival distributions \code{\link{fitParaSurv}}}
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(100);
#' # Weibull and Weibull, different means and medians
#' n = 1e3;
#' ## Target group
#' d1 = genData(n=n,dist="weibull",theta=c(1,1),p=0.2);
#' d1$arm = 1;
#' ## Reference group
#' d0 = genData(n=n,dist="weibull",theta=c(1,2),p=0.2);
#' d0$arm = 0;
#' ## Overall data set
#' data = rbind(d1,d0);
#' ## Comparison
#' comp = compParaSurv(time=data$time,status=data$status,arm=data$arm,tau=0.5,dist1="weibull");
#'
#' # Gamma and Weibull, same mean, different medians, bootstrap CIs
#' ## Target group
#' d1 = genData(n=n,dist="gamma",theta=c(2,2),p=0.2);
#' d1$arm = 1;
#' ## Reference group
#' d0 = genData(n=n,dist="weibull",theta=c(2,sqrt(pi)/2),p=0.2);
#' d0$arm = 0;
#' ## Overall data set
#' data = rbind(d1,d0);
#' ## Comparison
#' comp = compParaSurv(time=data$time,status=data$status,arm=data$arm,
#'                     tau=0.5,dist1="gamma",dist0="weibull",boot=2e3);
#'
#' # Weibull and Log-normal, different means, same median, permutation p.values
#' ## Target group
#' d1 = genData(n=n,dist="log-normal",theta=c(0,2),p=0.2);
#' d1$arm = 1;
#' ## Reference group
#' d0 = genData(n=n,dist="weibull",theta=c(2,sqrt(log(2))),p=0.2);
#' d0$arm = 0;
#' ## Overall data set
#' data = rbind(d1,d0);
#' ## Comparison
#' comp = compParaSurv(time=data$time,status=data$status,arm=data$arm,
#'                     tau=0.5,dist1="log-normal",dist0="weibull",perm=2e3);
#' }

compParaSurv = function(time,status=NULL,arm,tau=NULL,
                        dist1="weibull",dist0=NULL,
                        init1=NULL,init0=NULL,
                        boot=NULL,perm=NULL,
                        sig=0.05,eps=1e-6,maxit=10,report=F){
  
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
  
  ## Arm
  pass = checkArm(arm=arm,n=n);
  
  if(!pass){
    stop("Arm check failed.")
  }
  
  ## Distribution
  if(is.null(dist0)){
    dist0 = dist1;
  }
  
  pass = checkDist(dist=dist1)&checkDist(dist=dist0);
  
  if(!pass){
    stop("Distribution check failed.");
  }
  
  ## Initialization
  pass = checkInit(dist=dist1,init=init1) & checkInit(dist=dist0,init=init0);
  
  if(!pass){
    stop("Initialization check failed.");
  }

  # Partition data
  d0 = data.frame("time"=time[arm==0],"status"=status[arm==0]);
  d1 = data.frame("time"=time[arm==1],"status"=status[arm==1]);
  ## Fit to reference arm
  fit0 = fitParaSurv(time=d0$time,status=d0$status,dist=dist0,tau=tau,
                     sig=sig,init=init0,eps=eps,maxit=maxit,report=report);
  y0 = fit0@Outcome;
  ## Fit to treatment arm
  fit1 = fitParaSurv(time=d1$time,status=d1$status,dist=dist1,tau=tau,
                     sig=sig,init=init1,eps=eps,maxit=maxit,report=report);
  y1 = fit1@Outcome;

  # Contrasts
  ## Difference of means
  Contrast.mu = estDiff(t1=y1[1,2],s1=y1[1,3],t0=y0[1,2],s0=y0[1,3],sig=sig);
  ## Ratio of means
  Contrast.rmu = estRatio(t1=y1[1,2],s1=y1[1,3],t0=y0[1,2],s0=y0[1,3],sig=sig);
  ## Difference of medians
  Contrast.me = estDiff(t1=y1[2,2],s1=y1[2,3],t0=y0[2,2],s0=y0[2,3],sig=sig);
  ## Ratio of medians
  Contrast.rme = estRatio(t1=y1[2,2],s1=y1[2,3],t0=y0[2,2],s0=y0[2,3],sig=sig);
  
  ## Difference of RMSTs
  if(is.numeric(tau)){
    # RMSTs
    m = length(tau);
    r0 = fit0@RMST;
    r1 = fit1@RMST;
    
    aux = function(i){
      contrast.rmst = estDiff(t1=r1[i,2],s1=r1[i,3],t0=r0[i,2],s0=r0[i,3],sig=sig);
      contrast.rrmst = estRatio(t1=r1[i,2],s1=r1[i,3],t0=r0[i,2],s0=r0[i,3],sig=sig);
      out = rbind(contrast.rmst,contrast.rrmst);
      return(out);
    }
    
    RMST = lapply(seq(1:m),aux);
    RMST = do.call(rbind,RMST);
    
    # Formatting
    RMST$Tau = rep(tau,each=2);
    RMST$Contrast = rep(c("RMST1-RMST0","RMST1/RMST0"),times=m);
    RMST = RMST[,c(6:7,1:5)];
  }
  
  # Location frame
  Location = rbind(Contrast.mu,Contrast.me,Contrast.rmu,Contrast.rme);
  Location = cbind("Contrast"=c("Mean1-Mean0","Med1-Med0","Mean1/Mean0","Med1/Med0"),Location);
  
  # Bootstrap
  if(is.numeric(boot)){
    Boot = bootCI(d1=d1,d0=d0,tau=tau,dist1=dist1,dist0=dist0,B=boot,sig=sig);
    
    # Location
    Location.boot = Boot$Location;
    Location = merge(x=Location,y=Location.boot,by=c("Contrast"));
    
    # RMST
    if(is.numeric(tau)){
      RMST.boot = Boot$RMST;
      RMST = merge(x=RMST,y=RMST.boot,by=c("Tau","Contrast"));
    }
  }

  # Permutation
  if(is.numeric(perm)){
    
    Perm = permP(d1=d1,d0=d0,tau=tau,dist1=dist1,dist0=dist0,B=perm-1);
    
    # Location
    Location.perm = Perm$Location;
    Location = merge(x=Location,y=Location.perm,by=c("Contrast"));
    
    # RMST
    if(is.numeric(tau)){
      RMST.perm = Perm$RMST;
      RMST = merge(x=RMST,y=RMST.perm,by=c("Tau","Contrast"));
    }
  }

  # Output
  Location = Location[order(Location$Contrast),];
  Out = new(Class="contrast",Dist1=fit1@Distribution,Dist0=fit0@Distribution,
            Model1=fit1,Model0=fit0,Location=Location);
  if(is.numeric(tau)){
    RMST = RMST[order(RMST$Tau,RMST$Contrast),];
    Out@RMST = RMST;
  }
  
  return(Out);
}

########################
# Bootstrap CI
########################

#' Resample Data
#' 
#' Generates a resampled data.frame with replacement.
#' 
#' @param data Input data.farme
#' 
#' @return Data.frame resampled from the input. 

resample = function(data){
  n = nrow(data);
  key = sample(x=n,size=n,replace=T);
  out = data[key,];
  return(out);
}

#' Bootstrap CIs
#'
#' Generates bootstrap confidence intervals for location and RMST estimates.
#'
#' @param d1 Target data.frame containing time and status. 
#' @param d0 Reference data.frame containing time and status. 
#' @param tau Optional truncation times for calculating restricted mean survival
#'   time.
#' @param dist1 String, target distribution.
#' @param dist0 String, reference distribution.
#' @param B Number of resamples.
#' @param sig Significance level, for CIs.
#'
#' @return List containing data.frames with the bootstrap CIs for the location
#'   and RMST estimates.

bootCI = function(d1,d0,tau=NULL,dist1,dist0,B=2000,sig){
  
  # bootstrap
  aux = function(b){
    
    # bootstrap data sets
    b1 = resample(d1);
    b0 = resample(d0);
    
    # fit weibull
    f1 = try(fitParaSurv(time=b1$time,status=b1$status,tau=tau,dist=dist1),silent=T);
    f0 = try(fitParaSurv(time=b0$time,status=b0$status,tau=tau,dist=dist0),silent=T);
    
    # success
    key1 = (class(f1)!="try-error");
    key0 = (class(f0)!="try-error");
    
    out = NA;
    
    if(key1&key0){
      
      loc.diff.b = f1@Outcome$Estimate[1:2]-f0@Outcome$Estimate[1:2];
      loc.ratio.b = f1@Outcome$Estimate[1:2]/f0@Outcome$Estimate[1:2];
      out = c(loc.diff.b,loc.ratio.b);
      
      if(is.numeric(tau)){
        rmst.diff.b = f1@RMST$Estimate-f0@RMST$Estimate;
        rmst.ratio.b = f1@RMST$Estimate/f0@RMST$Estimate;
        out = c(out,rmst.diff.b,rmst.ratio.b);
      }
      
    }
    
    return(out);
  }
  
  boot = lapply(seq(1:B),aux);
  boot = do.call(rbind,boot);
  
  # Confidence intervals
  lower = function(x){as.numeric(quantile(x,sig/2,na.rm=TRUE))};
  L = apply(boot,2,lower);
  upper = function(x){as.numeric(quantile(x,1-sig/2,na.rm=TRUE))};
  U = apply(boot,2,upper);
  
  # Location
  Out = list();
  Location = data.frame("Contrast"=c("Mean1-Mean0","Med1-Med0","Mean1/Mean0","Med1/Med0"));
  Location$L.boot = L[1:4];
  Location$U.boot = U[1:4];
  Out$Location = Location;
  
  # RMST
  if(is.numeric(tau)){
    
    L = L[5:length(L)];
    U = U[5:length(U)];
    m = length(tau);
    
    RMST = data.frame("Tau"=rep(tau,times=2),"Contrast"=rep(c("RMST1-RMST0","RMST1/RMST0"),each=m));
    RMST$L.boot = L;
    RMST$U.boot = U;
    Out$RMST = RMST;
  }

  return(Out);
}

########################
# Permutation
########################

#' Permutation P Value
#' 
#' Calculates permutation p-values for location and RMST estimates.
#' 
#' @param d1 Target data.frame containing time and status. 
#' @param d0 Reference data.frame containing time and status. 
#' @param tau Optional truncation times for calculating restricted mean survival
#'   time.
#' @param dist1 String, target distribution.
#' @param dist0 String, reference distribution.
#' @param B Number of resamples.
#' 
#' @return List containing data.frames with the bootstrap CIs for the location
#'   and RMST estimates.

permP = function(d1,d0,tau,dist1,dist0,B=1999){
  
  # Observed values
  f1 = try(fitParaSurv(time=d1$time,status=d1$status,tau=tau,dist=dist1),silent=T);
  f0 = try(fitParaSurv(time=d0$time,status=d0$status,tau=tau,dist=dist0),silent=T);
  
  loc.diff = f1@Outcome$Estimate[1:2]-f0@Outcome$Estimate[1:2];
  loc.ratio = log(f1@Outcome$Estimate[1:2]/f0@Outcome$Estimate[1:2]);
  obs = c(loc.diff,loc.ratio);
  
  if(is.numeric(tau)){
    rmst.diff = f1@RMST$Estimate-f0@RMST$Estimate;
    rmst.ratio = log(f1@RMST$Estimate/f0@RMST$Estimate);
    obs = c(obs,rmst.diff,rmst.ratio);
  }
  
  # Data.frame to permute
  d1$arm = 1;
  d0$arm = 0;
  data.perm = rbind(d1,d0);
  n = nrow(data.perm);
  
  # Permute
  aux = function(i){
    
    # permute treatment
    data.perm$arm = data.perm$arm[sample(n,n,replace=F)]; 
    
    b1 = data.perm[data.perm$arm==1,];
    b0 = data.perm[data.perm$arm==0,];
    
    # fit distribution
    f1 = try(fitParaSurv(time=b1$time,status=b1$status,tau=tau,dist=dist1),silent=T);
    f0 = try(fitParaSurv(time=b0$time,status=b0$status,tau=tau,dist=dist0),silent=T);
    
    # success
    key1 = (class(f1)!="try-error");
    key0 = (class(f0)!="try-error");
    
    out = NA;
    
    if(key1&key0){
      
      loc.diff.b = f1@Outcome$Estimate[1:2]-f0@Outcome$Estimate[1:2];
      loc.ratio.b = log(f1@Outcome$Estimate[1:2]/f0@Outcome$Estimate[1:2]);
      out = c(loc.diff.b,loc.ratio.b);
      
      if(is.numeric(tau)){
        rmst.diff.b = f1@RMST$Estimate-f0@RMST$Estimate;
        rmst.ratio.b = log(f1@RMST$Estimate/f0@RMST$Estimate);
        out = c(out,rmst.diff.b,rmst.ratio.b);
      }
      
    }
    
    # return
    return(out);
  }
  
  perm = lapply(seq(1:B),aux);
  perm = do.call(rbind,perm);
  
  # Confidence intervals
  aux = function(i){
    out = (1+sum(abs(perm[,i])>=abs(obs[i])))/(1+B);
    return(out);
  }
  
  pvalues = lapply(seq(1:ncol(perm)),aux);
  pvalues = do.call(c,pvalues);
  
  # Location
  Out = list();
  Location = data.frame("Contrast"=c("Mean1-Mean0","Med1-Med0","Mean1/Mean0","Med1/Med0"));
  Location$P.perm = pvalues[1:4];
  Out$Location = Location;
  
  # RMST
  if(is.numeric(tau)){
    
    pvalues = pvalues[c(5:length(pvalues))];
    m = length(tau);
    
    RMST = data.frame("Tau"=rep(tau,times=2),"Contrast"=rep(c("RMST1-RMST0","RMST1/RMST0"),each=m));
    RMST$P.perm = pvalues;
    Out$RMST = RMST;
  }
  
  # Output
  return(Out);
}

