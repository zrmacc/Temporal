# Purpose: Functions to select inputs. 

########################
# Parameter Defaults
########################

#' Set Default Parameters
#' 
#' Function to select default parameter values for each distribution. 
#' 
#' @param dist String, distribution name.
#' 
#' @return Numeric parameter vector

defaultParam = function(dist){
  theta = c();
  if(dist=="exp"){
    # Rate
    theta[1] = 1;
  } else if(dist=="gamma"){
    # Shape
    theta[1] = 1;
    # Rate
    theta[2] = 1;
  } else if(dist=="gen-gamma"){
    # Shape
    theta[1] = theta[2] = 1;
    # Rate
    theta[3] = 1;
  } else if(dist=="log-normal"){
    # Location
    theta[1] = 0;
    # Scale
    theta[2] = 1;
  } else if(dist=="weibull"){
    # Shape
    theta[1] = 1;
    # Rate
    theta[2] = 1;
  }
  return(theta);
}

########################
# Check Arm
########################

#' Check Arm
#' 
#' Check whether treatment arm is properly formatted. 
#'
#' @param arm 0/1, treatment arm. 
#' @param n Integer, sample size. 
#' 
#' @return Logical indication of whether arm was properly formatted. 

checkArm = function(arm,n){
  
  pass = TRUE;
  
  if(length(arm)!=n){
    pass = FALSE;
    warning("Arm should have the same length as time.");
  }
  
  arm.levels = sort(unique(arm));
  n.arm.levels = length(arm.levels);
  
  if(n.arm.levels==1){
    pass = FALSE;
    warning("compParaSurv is for comparing arms. See fitParaSurv for single arm estimation.")
  } else if(n.arm.levels==2){
    if(!all.equal(arm.levels,c(0,1))){
      pass = FALSE;
      warning("Numeric 0,1 encoding is expected for arm.")
    }
  } else {
    pass = FALSE;
    warning("Arm should have exactly two levels, encoded 0,1. ")
  }
  
  return(pass);
}

########################
# Check Distribution
########################

#' Check Distribution
#' 
#' Check whether the distribution selected is available.
#'
#' @param dist String, distribution name.
#' 
#' @return Logical indication of whether the distribution was available. 

checkDist = function(dist){
  
  pass = TRUE;
  choices = c("exp","gamma","gen-gamma","log-normal","weibull");
  valid = (dist %in% choices);
  
  if(!valid){
    pass = FALSE;
    warning("Select distribution from among the following choices:\n",
            paste(choices,collapse=", "));
    
  }
  
  return(pass);
}

########################
# Check Initialization
########################

#' Check Initialization
#' 
#' Check whether the initialization is valid. 
#'
#' @param dist String, distribution name.
#' @param init Numeric vector, initialization.
#' 
#' @return Logical indication of whether the initialization was valid. 

checkInit = function(dist,init){
  
  pass = TRUE;
  
  # Only check initialization if not null
  if(!is.null(init)){
    
    # Check for NA
    if(any(is.na(init))){
      pass = FALSE;
      warning("All parameters must be initialized.");
    } else {
      len = length(init);
      # Exponential 
      if (dist=="exp"){
        warning("Exponential does not require initialization.");
      } else if(dist=="log-normal"){
        # Log normal
        if(len!=2){
          pass = FALSE;
          warning("Log-normal requires 2 parameters.");
        } else if(!init[2]>0){
          pass = FALSE;
          warning("Scale parameter must be positive.");
        }
      } else {
        # Remaining distributions
        ## Ensure positivity
        if(!all(init>0)){
          pass = FALSE;
          warning("Initialize all parameters at positive values");
        }
        ## Check length
        if(dist=="gamma"&len!=2){
          pass = FALSE;
          warning("Gamma requires 2 parameters.");
        } else if(dist=="gen-gamma"&len!=3){
          pass = FALSE;
          warning("Generalized gamma requires 3 parameters.");
        } else if(dist=="weibull"&len!=1){
          warning("Weibull initialization only requires the shape parameter.");
        }
      }
    }
  }
  
  return(pass);
}

########################
# Check Status
########################

#' Status Check
#'
#' Function to ensure the status indicator is properly formatted
#'
#' @param n Integer, sample size.
#' @param status 0/1 status indicator.
#'
#' @return Logical indicator of whether the status indicator was properly
#'   formatted.

checkStatus = function(n,status){
  
  pass = TRUE;
  
  if(length(status)!=n){
    pass = FALSE;
    warning("Status should have the same length as time.");
  }
  
  status.levels = sort(unique(status));
  n.status.levels = length(status.levels);
  
  if(n.status.levels==1){
    if(status.levels!=1){
      pass = FALSE;
      warning("If only one status level is present, all observations should have status 1.")
    }
  } else if(n.status.levels==2){
    if(!all.equal(status.levels,c(0,1))){
      pass = FALSE;
      warning("Numeric 0,1 encoding is expected for status.")
    }
  } else {
    pass = FALSE;
    warning("Status should have at most two levels, encoded 0,1. ")
  }
  
  return(pass);
}


########################
# Check Theta
########################

#' Check Theta
#'
#' Function to check the appropriate number of parameters are supplied for the
#' selected distribution.
#'
#' @param dist String, distribution.
#' @param theta Numeric, parameter vector.
#'
#' @return Logical indication of whether the appropriate number of parameters
#'   was provided.

checkTheta = function(dist,theta){
  
  pass = TRUE;
  len = length(theta);
  
  # Exponential
  if(dist=="exp" & len!=1){
    pass = FALSE;
    warning("Exponential requires a single rate parameter.");
  } else if(dist=="gamma" & len!=2){
    pass = FALSE;
    warning("Gamma requires a shape and a rate parameter.");
  } else if(dist=="gen-gamma" & len!=3){
    pass = FALSE; 
    warning("Generalized gamma requires two shape and one rate parameters.");
  } else if(dist=="log-normal" & len!=2){
    pass = FALSE;
    warning("Log-normal requires a location and a scale parameter.");
  } else if(dist=="weibull" & len!=2){
    pass = FALSE;
    warning("Weibull requires a shape and a rate parameter.");
  }
  
  # Output
  return(pass);
}


