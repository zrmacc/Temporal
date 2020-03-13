#' Survival Functions
#' 
#' Constructs the survival function for a parameter distribution.
#' 
#' The parameter vector theta should contain the following elements, in order,
#' according to the distribution:
#' \describe{
#'  \item{Exponential}{Rate \eqn{\lambda}.}
#'  \item{Gamma}{Shape \eqn{\alpha}, rate \eqn{\lambda}.}
#'  \item{Generalized Gamma}{Shape 1 \eqn{\alpha}, shape 2 \eqn{\beta}, rate \eqn{\lambda}.}
#'  \item{Log-Normal}{Locaion \eqn{\mu}, scale \eqn{\sigma}.}
#'  \item{Weibull}{Shape \eqn{\alpha}, rate \eqn{\lambda}.}
#' }
#' 
#' @param dist String, distribution name.
#' @param theta Numeric parameter vector.
#' 
#' @return Survival function. 
#' 
#' @importFrom expint gammainc
#' @export 

survFunc = function(dist,theta){
  
  # Input check
  ## Distribution
  pass = checkDist(dist);
  if(!pass){
    stop("Distribution check failed.");
  }
  
  ## Parameter
  pass = checkTheta(dist=dist,theta=theta);
  if(!pass){
    stop("Parameter check failed.");
  }
  
  S = NULL;
  
  # Define survival function
  if(dist=="exp"){
    l = theta[1];
    S = function(t){
      exp(-l*t);
    }
  } else if(dist=="gamma"){
    a = theta[1];
    l = theta[2];
    S = function(t){
      expint::gammainc(a,l*t)/gamma(a);
    }
  } else if(dist=="gen-gamma"){
    a = theta[1];
    b = theta[2];
    l = theta[3];
    S = function(t){
      expint::gammainc(a,(l*t)^b)/gamma(a);
    }
  } else if(dist=="log-normal"){
    m = theta[1];
    s = theta[2];
    S = function(t){
      pnorm(q=log(t),mean=m,sd=s,lower.tail=F);
    }
  } else if(dist=="weibull"){
    a = theta[1];
    l = theta[2];
    S = function(t){
      exp(-(l*t)^a);
    }
  }
  
  # Output
  return(S);
}
