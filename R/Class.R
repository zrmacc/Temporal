#' Fitted Survival Distribution
#'
#' Defines the object class returned by fitting functions.
#'
#' @slot Distribution Fitted distribution, string.
#' @slot Parameters Parameters, data.frame. 
#' @slot Information Information components, matrix. 
#' @slot Outcome Properties of the fitted distribution, data.frame. 
#' @slot RMST Estimated restricted mean survival times, data.frame
#' @slot S Fitted survival function, function.
#' @name fit-class
#' @rdname fit-class
#' @exportClass fit

setClass(Class="fit",representation=representation(Distribution="character",
                                                   Parameters="data.frame",
                                                   Information="matrix",
                                                   Outcome="data.frame",
                                                   RMST="data.frame",
                                                   S="function"));

########################
# Print Method
########################

#' Print Method for Fitted Survival Distributions
#'
#' Print method for objects of class \code{fit}.
#'
#' @param x An object of class \code{fit}.
#' @param ... Unused.
#' @export

print.fit = function(x,...){
  # Function to round data.frames
  aux = function(v){
    if(is.numeric(v)){return(signif(v,digits=3))}
    else{return(v)};
    };
  # Distribution
  dist = x@Distribution;
  # Parameters
  P = x@Parameters;
  P[] = lapply(P,aux);
  # Outcome characteristics
  Y = x@Outcome;
  Y[] = lapply(Y,aux);
  # RMST
  if(length(x@RMST)>0){
    R = x@RMST;
    R[] = lapply(R,aux);
  }
  
  # Display
  cat(paste0("Fitted ",dist," Distribution."),"\n");
  cat("Estimated Parameters:\n");
  print(P);
  cat("\n");
  cat("Distribution Properties:\n");
  print(Y);
  cat("\n");
  if(length(x@RMST)>0){
    cat("Restricted Mean Survival Times:\n");
    print(R);
    cat("\n");
  }
}

########################
# Show Method
########################

#' Show Method for Fitted Survival Distributions
#'
#' @param object An object of class \code{fit}.
#' @rdname fit-method
#' @importFrom methods show

setMethod(f="show",
          signature=c(object="fit"),
          definition=function(object){print.fit(x=object)});

########################
# Model Contrast
########################

#' Contrast of Survival Distributions.
#'
#' Defines the object class returned by the comparison function.
#'
#' @slot Dist1 Distribution fit to the target group, string.
#' @slot Dist0 Distribution fit to the reference group, string.
#' @slot Model1 Fitted model for the target group, fit.
#' @slot Model0 Fitted model for the reference group, fit.
#' @slot Location Contrasts of means and medians, data.frame.
#' @slot RMST Contrasts of RMSTs, data.frame. 
#' @name contrast-class
#' @rdname contrast-class
#' @exportClass contrast

setClass(Class="contrast",
         representation=representation(Dist1="character",
                                       Dist0="character",
                                       Model1="fit",
                                       Model0="fit",
                                       Location="data.frame",
                                       RMST="data.frame"));

########################
# Print Method
########################

#' Print Method for a Contrast of Survival Distributions.
#'
#' Print method for an object of class \code{contrast}.
#'
#' @param x A \code{contrast} object.
#' @param ... Unused.
#' @export

print.contrast = function(x,...){
  # Function to round data.frames
  aux = function(v){
    if(is.numeric(v)){return(signif(v,digits=3))}
    else{return(v)};
  };
  # Distribution
  dist1 = x@Dist1;
  dist0 = x@Dist0;
  flag = (dist1==dist0);
  # Group 1
  G1 = x@Model1@Outcome;
  G1[] = lapply(G1,aux);
  # Group 0
  G0 = x@Model0@Outcome;
  G0[] = lapply(G0,aux);
  # Location
  Location = x@Location;
  Location[] = lapply(Location,aux);
  # RMST
  if(length(x@RMST>0)){
    RMST = x@RMST;
    RMST[] = lapply(RMST,aux);
  }
  
  # Display
  if(flag){
    cat(paste0("Contrast of Fitted ",dist1," Distributions."),"\n\n");
  } else {
    cat(paste0("Contrast of Fitted ",dist1," v. ",dist0),"\n\n");
  }
  cat("Fitted Characteristics for Group 1:\n");
  print(G1);
  cat("\n");
  cat("Fitted Characteristics for Group 0:\n");
  print(G0);
  cat("\n");
  cat("Location:\n");
  print(Location);
  cat("\n");
  if(length(x@RMST>0)){
    cat("RMST:\n");
    print(RMST);
    cat("\n");
  }
}

########################
# Show Method
########################

#' Show Method for a Contrast of Survival Distributions.
#'
#' @param object An object of class \code{contrast}.
#' @rdname contrast-method
#' @importFrom methods show

setMethod(f="show",
          signature=c(object="contrast"),
          definition=function(object){print.contrast(x=object)});
