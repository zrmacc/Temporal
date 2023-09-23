# Purpose: Comparison of two fitted, parametric survival distributions.
# Updated: 2021-07-18

# -----------------------------------------------------------------------------
# Analytical contrast of two arms.
# -----------------------------------------------------------------------------

#' Difference of Estimates
#'
#' Calculate CIs and p-value for the difference of estimates.
#'
#' @param est1 Arm 1 estimate.
#' @param se1 Arm 1 standard error.
#' @param est0 Arm 0 estimate.
#' @param se0 Arm 0 standard error.
#' @param sig Significance level.
#' @return Data.frame containing estimated difference, its standard error, lower and
#'   upper confidence bounds, and a p-value assessing the null hypothesis of no
#'   difference.

EstDiff <- function(est1, se1, est0, se0, sig = 0.05) {
  z <- stats::qnorm(1 - sig / 2)
  
  # Difference and standard error.
  delta <- est1 - est0
  se <- sqrt(se1^2 + se0^2)
  
  # CI and p-value.
  lower <- delta - z * se
  upper <- delta + z * se
  p <- 2 * stats::pnorm(abs(delta / se), lower.tail = FALSE)
  
  # Output.
  out <- data.frame(
    Point = delta,
    SE = se,
    L = lower,
    U = upper,
    P = p
  )
  return(out)
}


#' Ratio of Estimates
#'
#' Calculate CIs and p-value for the ratio of estimates.
#'
#' @param est1 Arm 1 estimate.
#' @param se1 Arm 1 standard error.
#' @param est0 Arm 0 estimate.
#' @param se0 Arm 0 standard error.
#' @param sig Significance level.
#' @return Data.frame containing estimated ratio, its standard error, lower and
#'   upper confidence bounds, and a p-value assessing the null hypothesis that
#'   the ratio is unity.

EstRatio <- function(est1, se1, est0, se0, sig = 0.05) {
  z <- stats::qnorm(1 - sig / 2)
  
  # Ratio and standard error.
  log_rho <- log(est1) - log(est0)
  log_rho_se <- sqrt(se1^2 / (est1^2) + se0^2 / (est0^2))
  rho_se <- sqrt(se1^2 / (est0^2) + (est1^2 * se0^2) / (est0^4))
  
  # CI and p-value.
  lower <- exp(log_rho - z * log_rho_se)
  upper <- exp(log_rho + z * log_rho_se)
  p <- 2 * stats::pnorm(abs(log_rho / log_rho_se), lower.tail = FALSE)
  
  # Output.
  out <- data.frame(
    Point = exp(log_rho), 
    SE = rho_se, 
    L = lower, 
    U = upper, 
    P = p
  )
  return(out)
}


# -----------------------------------------------------------------------------
# Compare RMSTs
# -----------------------------------------------------------------------------

#' Contrast Locations
#' 
#' Compare the means and medians of the fitted distributions
#' for two treatment arms. 
#' @param fit1 Fitted parametric survival distribution for arm 1.
#' @param fit0 Fitted parametric survival distribution for arm 0.
#' @param sig Significance level.
#' @return Data.frame contrasting the difference and ratio of the
#'   mean and median at each time point.

ContrastLocs <- function(fit1, fit0, sig = 0.05) {
  
  # Contrast means.
  mu1 <- fit1@Outcome$Estimate[1]
  mu1_se <- fit1@Outcome$SE[1]
  mu0 <- fit0@Outcome$Estimate[1]
  mu0_se <- fit0@Outcome$SE[1]
  
  diff_means <- EstDiff(mu1, mu1_se, mu0, mu0_se, sig = sig)
  ratio_means <- EstRatio(mu1, mu1_se, mu0, mu0_se, sig = sig)
  
  # Contrast medians.
  me1 <- fit1@Outcome$Estimate[2]
  me1_se <- fit1@Outcome$SE[2]
  me0 <- fit0@Outcome$Estimate[2]
  me0_se <- fit0@Outcome$SE[2]
  
  diff_meds <- EstDiff(me1, me1_se, me0, me0_se, sig = sig)
  ratio_meds <- EstRatio(me1, me1_se, me0, me0_se, sig = sig)
  
  # Output.
  out <- data.frame(
    Contrast = c("Mean1-Mean0", "Mean1/Mean0", "Med1-Med0", "Med1/Med0")
  )
  out <- cbind(
    out, 
    rbind(diff_means, ratio_means, diff_meds, ratio_meds)
  )
  return(out)
}


#' Contrast RMSTs
#' 
#' Compare the restricted mean survival times of the fitted distributions
#' for two treatment arms.
#' 
#' @param fit1 Fitted parametric survival distribution for arm 1.
#' @param fit0 Fitted parametric survival distribution for arm 0.
#' @param sig Significance level, for 
#' @return Data.frame contrasting the difference and ratio of RMSTs
#'   at each time point.
#' @importFrom dplyr "%>%"

ContrastRMSTs <- function(fit1, fit0, sig = 0.05) {
  
  arm <- NULL
  Estimate <- NULL
  SE <- NULL
  Tau <- NULL
  
  rmst1 <- fit1@RMST %>% 
    dplyr::select(Tau, Estimate, SE) %>%
    dplyr::mutate(arm = 1)
  
  rmst0 <- fit0@RMST %>% 
    dplyr::select(Tau, Estimate, SE) %>%
    dplyr::mutate(arm = 0)
  
  rmst <- rbind(rmst1, rmst0) %>%
    tidyr::pivot_wider(
      id_cols = Tau, 
      names_from = arm,
      names_sep = "",
      values_from = c("Estimate", "SE")
    )
  
  split_rmst <- split(x = rmst, f = seq(1:nrow(rmst)))
  contrast_rmst <- lapply(split_rmst, function(x) {
    diff <- EstDiff(x$Estimate1, x$SE1, x$Estimate0, x$SE0, sig = sig)
    ratio <- EstRatio(x$Estimate1, x$SE1, x$Estimate0, x$SE0, sig = sig)
    out <- data.frame(
      Tau = x$Tau,
      Contrast = c("RMST1-RMST0", "RMST1/RMST0")
    )
    out <- cbind(out, rbind(diff, ratio))
    return(out)
  })
  out <- do.call(rbind, contrast_rmst)
  return(out)  
}

# -----------------------------------------------------------------------------
# Permutation
# -----------------------------------------------------------------------------

#' Extract Observed Estimates
#' 
#' Helper function for permutation inference.
#' 
#' @param fit1 Fitted parametric survival distribution for arm 1.
#' @param fit0 Fitted parametric survival distribution for arm 0.
#' @return Numeric vector.

ExtractObsEst <- function(fit1, fit0) {
  
  # Means.
  mu1 <- fit1@Outcome$Estimate[1]
  mu0 <- fit0@Outcome$Estimate[1]
  mu_diff <- mu1 - mu0
  mu_ratio <- mu1 / mu0
  
  # Medians.
  me0 <- fit0@Outcome$Estimate[2]
  me1 <- fit1@Outcome$Estimate[2]
  me_diff <- me1 - me0
  me_ratio <- me1 / me0
  
  # Observed values.
  obs <- c(
    mu_diff = mu_diff,
    mu_ratio_log = log(mu_ratio),
    me_diff = me_diff,
    me_ratio_log = log(me_ratio)
  )
  
  # RMST.
  if (nrow(fit1@RMST) > 0) {
    rmst0 <- fit0@RMST$Estimate
    rmst1 <- fit1@RMST$Estimate
    
    rmst_diff <- rmst1 - rmst0
    names(rmst_diff) <- paste0("tau_", fit0@RMST$Tau, "_diff")
    
    rmst_ratio <- rmst1 / rmst0
    names(rmst_ratio) <- paste0("tau_", fit0@RMST$Tau, "_ratio_log")
    
    interleave <- c(rbind(rmst_diff, log(rmst_ratio))) 
    obs <- c(obs, interleave)
  }
  
  return(obs)
}


#' Permutation P Value
#'
#' Calculates permutation p-values for location and RMST estimates.
#'
#' @param df1 Target data.frame containing time and status.
#' @param df0 Reference data.frame containing time and status.
#' @param fit1 Fitted parametric survival distribution for arm 1.
#' @param fit0 Fitted parametric survival distribution for arm 0.
#' @param eps Tolerance for Newton-Raphson iterations.
#' @param init1 Initial parameter values for the target group.
#' @param init0 Initial parameter values for the reference group.
#' @param maxit Maximum number of Newton-Raphson iterations.
#' @param reps Number of permutation replicates.
#' @param tau Optional truncation times for calculating RMST.
#' @return Numeric vector of permutation p-values.

PermP <- function(df1, df0, fit0, fit1, eps, init1, init0, maxit, reps, tau) {
  
  # Extract distributions.
  dist0 <- fit0@Distribution
  dist1 <- fit1@Distribution
  obs <- ExtractObsEst(fit1, fit0)
  
  # Data.frame to permute.
  df1$arm <- 1
  df0$arm <- 0
  df_orig  <- rbind(df1, df0)
  n <- nrow(df_orig)
  
  # Permutations.
  perm <- lapply(seq_len(reps), function(i) {
    
    # Permute treatment arm.
    df_perm <- df_orig
    df_perm$arm <- df_orig$arm[sample(n, n, replace = FALSE)]
    
    # Fit distributions.
    arm <- NULL 
    
    fit1_perm <- try(
      FitParaSurv(
        df_perm %>% dplyr::filter(arm == 1),
        dist = dist1,
        eps = eps,
        init = init1,
        maxit = maxit,
        tau = tau), 
      silent = TRUE
    )
    fit1_failed <- methods::is(fit1_perm, "try-error")
   
    fit0_perm <- try(
       FitParaSurv(
        df_perm %>% dplyr::filter(arm == 0),
        dist = dist0,
        eps = eps,
        init = init0,
        maxit = maxit,
        tau = tau), 
      silent = TRUE
    )
    fit0_failed <- methods::is(fit0_perm, "try-error")
    
    if (!fit1_failed & !fit0_failed) {
      out <- ExtractObsEst(fit1_perm, fit0_perm)
    } else {
      out <- NULL
    }
    return(out)
  })
  
  perm <- do.call(rbind, perm)
  if (is.null(perm)){
    warning("Calculation of permutation p-value failed.")
    return(NULL)
  }
  
  # P-values.
  pvalues <- sapply(seq_len(ncol(perm)), function(i) {
    return((1 + sum(abs(perm[, i]) >= abs(obs[i]))) / (1 + nrow(perm)))
  })
  return(pvalues)
}


# -----------------------------------------------------------------------------
# Distribution Comparison
# -----------------------------------------------------------------------------

#' Compare Parametric Survival Distribution
#'
#' Compares the means and medians of parametric survival distributions fit to
#' two treatment arms. Available distributions include: exponential, gamma,
#' generalized gamma, log-normal, and Weibull.
#'
#' Status should be coded as 0 for censored and 1 for observed. Arm is coded as
#' 0 for reference, 1 for target. Tau is an optional numeric vector of
#' truncation times for calculating restricted mean survival time, which is the
#' area under the survival curve up to the specified truncation point.
#'
#' @param data Data.frame.
#' @param arm_name Name of the column containing the treatment group, coded as 1 for
#'   treatment, 0 for reference.
#' @param dist1 Distribution to fit for the target group. Selected from among:
#'   exp, gamma, gengamma, log-normal, and weibull.
#' @param dist0 Distribution to fit for the reference group. Same choices as for
#'   the target group. If omitted, defaults to the distribution specified for
#'   the target group.
#' @param eps Tolerance for Newton-Raphson iterations.
#' @param init1 Initial parameter values for the target group.
#' @param init0 Initial parameter values for the reference group.
#' @param maxit Maximum number of Newton-Raphson iterations.
#' @param report Report fitting progress?
#' @param reps Number of permutation replicates, if requesting permutation
#'   p-values.
#' @param sig Significance level, for constructing confidence intervals.
#' @param status_name Name of the status indicator, 1 if observed, 0 if
#'   censored.
#' @param tau Optional truncation times for calculating RMST.
#' @param time_name Name of column containing the time to event.
#' @return An object of class \code{contrast} containing the following:
#' \describe{
#'   \item{Model1}{The fitted model for the target group.}
#'   \item{Model0}{The fitted model for the reference group.}
#'   \item{Contrast}{Contrasts of means and medians.}
#'   \item{RMST}{Contrasts of the RMSTs, if `tau` was specified.}
#' }
#' @export
#' @examples
#' \donttest{
#' set.seed(100)
#' # Weibull and Weibull, different means and medians.
#' n <- 1e3
#' 
#' # Generate data.
#' df1 <- GenData(n = n, dist = "weibull", theta = c(1, 1), p = 0.2)
#' df1$arm <- 1
#' df0 <- GenData(n = n, dist = "weibull", theta = c(1, 2), p = 0.2)
#' df0$arm <- 0
#' data <- rbind(df1, df0)
#' 
#' # Comparison.
#' comp <- CompParaSurv(data, dist1 = "weibull")
#' 
#' # Add RMST at time 1.
#' comp <- CompParaSurv(data, dist1 = "weibull", tau = 1)
#' 
#' # Calculate permutation p-values (slow).
#' comp <-  CompParaSurv(data, dist1 = "weibull", tau = 1, reps = 100)
#' }
CompParaSurv <- function(
  data,
  arm_name = "arm",
  dist1 = "weibull", 
  dist0 = NULL,
  eps = 1e-6,
  init1 = NULL, 
  init0 = NULL,
  maxit = 10,
  report = FALSE,
  reps = NULL,
  sig = 0.05, 
  status_name = "status",
  tau = NULL,
  time_name = "time"
) {
  
  # Formatting.
  arm <- NULL
  status <- NULL
  time <- NULL
  df <- data %>%
    dplyr::rename(
      arm = {{ arm_name }},
      status = {{ status_name }},
      time = {{ time_name }}
    )

  # Positivity.
  if (min(df$time) < 0) {
    stop("Strictly positive observation times required.")
  }

  # Input Checks.
  CheckArm(df$arm)
  CheckStatus(df$status)
  if (is.null(dist0)) {
    dist0 <- dist1
  }
  CheckDist(dist1) 
  CheckDist(dist0)
  CheckInit(dist1, init1)
  CheckInit(dist0, init0)
  
  # Partition data.
  df0 <- df %>% dplyr::filter(arm == 0)
  df1 <- df %>% dplyr::filter(arm == 1)
  
  # Fit to reference arm.
  fit0 <- FitParaSurv(
    data = df0, 
    dist = dist0, 
    eps = eps,
    init = init0,
    maxit = maxit,
    report = report,
    sig = sig,
    tau = tau
  )
  
  # Fit to treatment arm.
  fit1 <- FitParaSurv(
    data = df1, 
    dist = dist1, 
    eps = eps,
    init = init1,
    maxit = maxit,
    report = report,
    sig = sig,
    tau = tau
  )

  # Contrasts.
  loc_contrast <- ContrastLocs(fit1, fit0, sig = sig)
  
  # RMSTs.
  if (!is.null(tau)) {
    rmst_contrast <- ContrastRMSTs(fit1, fit0, sig = sig)
  }
  
  # Permutation.
  if (!is.null(reps)) {
    perm_p <- PermP(df0, df1, fit0, fit1, eps, init1, init0, maxit, reps, tau)
    loc_contrast$Perm_P <- perm_p[1:4]
    if (!is.null(tau)) {
      rmst_contrast$Perm_P <- perm_p[5:length(perm_p)]
    }
  }

  # Output.
  out <- methods::new(
    Class = "contrast", 
    Dist1 = fit1@Distribution, 
    Dist0 = fit0@Distribution,
    Model1 = fit1, 
    Model0 = fit0,
    Location = loc_contrast
  )
  if (!is.null(tau)) {
    out@RMST <- rmst_contrast
  }
  return(out)
}