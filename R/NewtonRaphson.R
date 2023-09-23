# Purpose: General purpose Newton-Raphson.
# Updated: 2021-07-12

#' Quadratic Form
#' 
#' @param x Numeric vector.
#' @param A Numeric matrix.
#' @return Numeric scalar.
QF <- function(x, A) {
  x <- matrix(x, ncol = 1)
  out <- t(x) %*% A %*% x
  return(as.numeric(out))
}


# -----------------------------------------------------------------------------

#' Newton Raphson Update Iteration
#' 
#' @param obj Objective function.
#' @param state List containing the parameter vector `theta`.
#' @return List containing the updated parameter vector `theta` and the
#'   objective increment `delta`.
NRUpdate <- function(obj, state) {
  theta <- state$theta
  
  # Initial objective.
  obj0 <- obj(theta)
  
  # Score.
  score <- matrix(numDeriv::grad(obj, theta), ncol = 1)
  
  # Hessian.
  hess <- numDeriv::hessian(obj, theta)
  
  # Proposal.
  prop <- theta - as.numeric(solve(hess, score))
  
  # Final objective.
  obj1 <- obj(prop)
  
  # Output state.
  out <- list(
    theta = prop,
    delta = obj1 - obj0
  )
  return(out)
}


#' Newton Raphson Estimation
#' 
#' @param init Initial value.
#' @param obj Objective function.
#' @param eps Tolerance for Newton-Raphson iterations.
#' @param maxit Maximum number of NR iterations.
#' @param report Report fitting progress?
#' @return Numeric parameter estimate.
#' @export
NewtonRaphson <- function(
  init,
  obj,
  eps = 1e-6, 
  maxit = 10, 
  report = FALSE
) {
  
  # Maximization.
  state <- list(theta = init)
  for (i in 1:maxit) {
    new_state <- NRUpdate(obj, state)
    
    # Accept update if increment is positive.
    if (new_state$delta > 0) {
      state <- new_state
      if (report) {
        cat("Objective increment: ", signif(new_state$d, digits = 3), "\n")
      }
    }
    
    # Terminate if increment is below tolerance.
    if (new_state$delta < eps) {break}
  }
  
  ## Fitting report
  if (report) {
    if (i < maxit) {
      cat(paste0(i - 1, " update(s) performed before tolerance limit.\n\n"))
    } else {
      cat(paste0(i, " update(s) performed without reaching tolerance limit.\n\n"))
    }
  }

  return(state$theta)
}
