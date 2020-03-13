# Purpose: Function to calculate RMST
# Updated: 20/03/11

#' Restricted Mean Survival Time
#' 
#' Calculates the tau-year RMST for a fitted parametric model. 
#' 
#' @param fit Fitted parametric survival distribution. 
#' @param sig Significance level, for CIs.
#' @param tau Numeric vector of truncation times. 
#' 
#' @return Data.table containing the estimated RMST at each truncation time.
#' 
#' @importFrom stats integrate
#' @importFrom numDeriv grad
#' @export 
#' 
#' @examples 
#' # Generate Weibull data with 20% censoring.
#' data = genData(n=1e3,dist="weibull",theta=c(2,0.5),p=0.2);
#' fit = fitParaSurv(time=data$time,status=data$status,dist="weibull");
#' rmst = paraRMST(fit=fit,tau=c(0.5,1.0,1.5,2.0));
#' 
#' # Generate Gamma data with 10% censoring.
#' data = genData(n=1e3,dist="gamma",theta=c(2,2),p=0.10);
#' fit = fitParaSurv(time=data$time,status=data$status,dist="gamma");
#' rmst = paraRMST(fit=fit,tau=c(0.5,1.0,1.5,2.0));

paraRMST = function(fit,sig=0.05,tau){
  
  # Input check
  if(class(fit)!="fit"){
    stop("Requires a fitted parametric survival distribution.");
  }
  if(!is.numeric(tau)|(min(tau)<0)){
    stop("Requires positive, numeric truncation times.");
  }
  
  ## Number of truncation times
  m = length(tau);
  ## Critical value
  z = qnorm(1-sig/2);
  ## Observed parameters
  theta.obs = fit@Parameters$Estimate;
  ## Inverse information
  Ji = matInv(fit@Information);
  ## RMST function
  rmst = function(theta,U){
    S = survFunc(dist=fit@Distribution,theta=theta);
    out = integrate(f=S,lower=0,upper=U)$value;
    return(out);
  };
  
  # Loop
  aux = function(i){
    ## Truncation time
    out = data.frame("Tau"=tau[i]);
    ## RMST as a function of theta.
    aux = function(theta){rmst(theta,U=tau[i])};
    ## Observed RMST
    out$Estimate = aux(theta.obs);
    ## SE
    dg = grad(func=aux,x=theta.obs);
    out$SE = as.numeric(sqrt(matQF(X=matrix(dg,ncol=1),A=Ji)));
    ## Output
    return(out);
  }
  
  Out = lapply(seq(1:m),aux);
  Out = do.call(rbind,Out);
  
  Out$L = Out$Estimate-z*Out$SE;
  Out$U = Out$Estimate+z*Out$SE;
  
  return(Out);
}

