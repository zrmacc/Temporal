#' Model Fit
#'
#' Defines the object class returned by fitting functions.
#'
#' @slot Distribution Fitted distribution.
#' @slot Parameters Parameters.
#' @slot Information Information components.
#' @slot Outcome Properties of fitted distribution.
#' @name fit-class
#' @rdname fit-class
#' @exportClass fit

setClass(Class="fit",representation=representation(Distribution="character",Parameters="data.frame",Information="matrix",Outcome="data.frame"));

########################
# Print Method
########################

#' Print for Fitted Survival Models
#'
#' Print method for objects of class \code{fit}.
#'
#' @param x A \code{fit} object.
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
  # Display
  cat(paste0("Fitted ",dist," Distribution."),"\n");
  cat("Estimated Parameters:\n");
  print(P);
  cat("\n");
  cat("Distribution Properties:\n");
  print(Y);
  cat("\n");
}

########################
# Show Method
########################

#' Show for Fitted Survival Models
#' @param object A \code{fit} object.
#' @rdname fit-method
#' @importFrom methods show

setMethod(f="show",signature=c(object="fit"),definition=function(object){print.fit(x=object)});

########################
# Model Contrast
########################

#' Model Contrast
#'
#' Defines the object class returned by fitting functions.
#'
#' @slot Dist1 Distribution fit to the target group.
#' @slot Dist0 Distribution fit to the reference group.
#' @slot Model1 Fitted model for the target group.
#' @slot Model0 Fitted model for the reference group. 
#' @slot Contrast Model contrasts. 
#' @name contrast-class
#' @rdname contrast-class
#' @exportClass contrast

setClass(Class="contrast",representation=representation(Dist1="character",Dist0="character",Model1="fit",Model0="fit",Contrast="data.frame"));

########################
# Print Method
########################

#' Print for Contrasts
#'
#' Print method for objects of class \code{contrast}.
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
  # Contrasts
  Contrast = x@Contrast;
  Contrast[] = lapply(Contrast,aux);
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
  cat("Contrasts:\n");
  print(Contrast);
  cat("\n");
}

########################
# Show Method
########################

#' Show for Contrasts
#' @param object A \code{contrast} object.
#' @rdname contrast-method
#' @importFrom methods show

setMethod(f="show",signature=c(object="contrast"),definition=function(object){print.contrast(x=object)});
