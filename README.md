---
title: "README"
author: "Zachary McCaw"
date: "2020-03-13"
output: 
  html_document: 
    keep_md: TRUE
--- 

# Package Vignette




# Contents

* [Introduction](#introduction)
* [Simulation](#simulation)
* [Estimation](#estimation)
* [Group Contrasts](#group-contrasts)

# Introduction

## Overview

This package performs estimation and inference on parametric survival curves. Supported distributions include the exponential, gamma, generalized gamma, log-normal, and Weibull. Data are expected in time-to-event format, with a status indicator to accommodate non-informative right censoring. The function `fitParaSurv` provides maximum likelihood estimates (MLEs) and confidence intervals (CIs) for the distribution of interest. Estimates are presented for model parameters, and for characteristics of the event time distribution (the mean, median, variance, and restricted mean). The function `compParaSurv` compares the fitted survival distributions of two treatment arms with respect to the difference and ratio of means, medians, and restricted means. For each contrast, point estimates, CIs, and p-values assessing the null hypothesis of no between-group difference are presented. While by default the between-group difference is evaluated using asymptotic CIs and p-values, options are provided for estimating the CIs via bootstrap, and the p-values via permutation. 

## Setting

Suppose the data consist of pairs $(U_{i},\delta_{i})$, where the observation time $U_{i} = \min(T_{i},C_{i})$ is the minimum of an event time $T_{i}$ and a censoring time $C_{i}$; and $\delta_{i}=I(T_{i}\leq C_{i})$ is an indicator that an event was observed. The event times $T_{i}$ are assumed to follow a survival distribution $S_{T}$ parameterized by $\theta$. The distribution of censoring times $C_{i}$ is left unspecified. Estimation of $\theta$ proceeds by maximizing the right-censored log likelihood $\ell(\theta) = \sum_{i=1}^{n}\delta_{i}\ln h_{i} + \ln S_{i}$, where $S_{i}$ and $h_{i}$ denote the survival and hazard contributions of the $i$th subject. The asymptotic covariance of the MLE $\hat{\theta}$ is estimated using the (inverse) observed information, and standard errors for functions of the MLE (e.g. means and medians) are obtained from the $\Delta$-method.

# Estimation

## Overview {#fit}

The function `fitParaSurv` requires the observation times `time` and status indicators `status` as input. `status` is coded as 1 if an event was observed, and as 0 if censored. The distribution of interest is specified using `dist`, which defaults to Weibull. The output is an object of class `fit`. The slot `fit@Parameters` contains the parameter estimates and CIs. The slot `fit@Information` contains the observed sample information matrix. The slot `fit@Outcome` contains characteristics of the fitted event time distribution. If truncation times `tau` are provided, the corresponding RMSTs are stored in `fit@RMST`. Functions are presented below for simulating from and estimating the exponential, gamma, generalized gamma, log-normal, and Weibull distributions. If an expected censoring proportion `p` is provided, then the event time are subject to non-informative random right censoring. 

## Exponential

The exponential distribution is parameterized in terms of the rate $\lambda$. The density is:

$$
f(t) = \lambda e^{-\lambda t},\ t>0
$$

In the following, $n=10^{3}$ exponential event times are simulated, with rate $\lambda=2$ and expected censoring proportion 20%. Generative parameters are recovered using `fitParaSurv` with `dist="exp"`.


```r
# Generate exponential time to event data
data = genData(n=1e3,dist="exp",theta=c(2),p=0.2)
# Estimate parameters
fit= fitParaSurv(time=data$time,status=data$status,dist="exp");
show(fit);
```

```
## Fitted exp Distribution. 
## Estimated Parameters:
##   Aspect Estimate     SE    L    U
## 1   Rate     2.12 0.0741 1.98 2.27
## 
## Distribution Properties:
##     Aspect Estimate     SE     L     U
## 1     Mean    0.471 0.0164 0.439 0.503
## 2   Median    0.327 0.0114 0.304 0.349
## 3 Variance    0.222 0.0155 0.192 0.252
```

## Gamma 

The gamma distribution is parameterized in terms of the shape $\alpha$ and rate $\lambda$. The density is:

$$
f(t) = \frac{\lambda}{\Gamma(\alpha)}(\lambda t)^{\alpha-1}e^{-\lambda t},\ t>0
$$

In the following, $n=10^{3}$ gamma event times are simulated, with shape $\alpha=2$, rate $\lambda=2$, and expected censoring proportion 25%. Generative parameters are recovered using `fitParaSurv` with `dist="gamma"`. RMST at $\tau=0.5$ is requested.


```r
# Generate gamma time to event data
data = genData(n=1e3,dist="gamma",theta=c(2,2),p=0.25);
# Estimate parameters
fit = fitParaSurv(time=data$time,status=data$status,dist="gamma",tau=0.5);
show(fit);
```

```
## Fitted gamma Distribution. 
## Estimated Parameters:
##   Aspect Estimate     SE    L    U
## a  Shape     1.90 0.0854 1.74 2.07
## l   Rate     1.85 0.1020 1.67 2.07
## 
## Distribution Properties:
##     Aspect Estimate     SE     L     U
## 1     Mean    1.020 0.0267 0.972 1.080
## 2   Median    0.851 0.0224 0.807 0.895
## 3 Variance    0.552 0.0404 0.473 0.631
## 
## Restricted Mean Survival Times:
##   Tau Estimate      SE    L     U
## 1 0.5    0.447 0.00335 0.44 0.453
```

## Generalized Gamma

The generalized gamma distribution is parameterized in terms of shapes $\alpha$ and $\beta$, and rate $\lambda$. The density is:

$$
f(t) = \frac{\beta\lambda}{\Gamma(\alpha)}(\lambda t)^{\alpha\beta-1}e^{-(\lambda t)^{\beta}},\ t>0
$$

The standard gamma and Weibull distributions are nested within the generalized gamma. Setting $\beta=1$ recovers the standard gamma, while setting $\alpha=1$ recovers the Weibull. In the following, $n=10^{4}$ generalized gamma event times are simulated, with shapes $\alpha=2$ and $\beta=2$, rate $\lambda=2$, and expected censoring proportion 10%. Generative parameters are recovered using `fitParaSurv` with `dist="gen-gamma"`. Accurate estimation of the generalized gamma generally requires larger sample sizes. Fitting progress is monitored by setting `report=T`. 


```r
set.seed(102);
# Generate generalized gamma time to event data
data = genData(n=1e4,dist="gen-gamma",theta=c(2,2,2),p=0.1);
# Estimate parameters
fit = fitParaSurv(time=data$time,status=data$status,dist="gen-gamma",report=T);
show(fit);
```

```
## Objective increment:  29.5 
## Objective increment:  0.103 
## Objective increment:  6.59e-06 
## Objective increment:  1.18e-11 
## 3 update(s) performed before tolerance limit.
## 
## Fitted gen-gamma Distribution. 
## Estimated Parameters:
##   Aspect Estimate     SE    L    U
## a ShapeA     2.02 0.1460 1.75 2.33
## b ShapeB     2.00 0.0833 1.84 2.17
## l   Rate     2.02 0.1120 1.81 2.25
## 
## Distribution Properties:
##     Aspect Estimate       SE      L     U
## a     Mean   0.6610 0.002470 0.6560 0.666
## b   Median   0.6450 0.002710 0.6390 0.650
## l Variance   0.0572 0.000892 0.0555 0.059
```

The final parameter estimates for the generalized gamma distribution may be sensitive to the initial values. If the Newton-Raphson iteration is initialized too far from the optimum of the log-likelihood, the search may not reach the maximum. If the search halts prematurely, the fitting procedure will indicate that the information matrix was not positive definite. Supplying initial parameter values, as in the following, may help to improve the final estimates.


```r
# Initialization
fitParaSurv(time=data$time,status=data$status,dist="gen-gamma",init=c(2,2,2));
```

## Log-Normal

The log-normal distribution is parameterized in terms of the location $\mu$ and scale $\sigma$. The density is:

$$
f(t) = \frac{1}{t\sigma\sqrt{2\pi}}e^{-\frac{(\ln t-\mu)^2}{2\sigma^2}},\ t>0
$$

In the following, $n=10^{3}$ log-normal event times are simulated, with location $\mu=1$, scale $\sigma=2$, and expected censoring proportion 15%. Generative parameters are recovered using `fitParaSurv` with `dist="log-normal"`. RMSTs at $\tau=(5,10,25)$ are requested.


```r
# Generate log-normal time to event data
data = genData(n=1e3,dist="log-normal",theta=c(1,2),p=0.15);
# Estimate parameters
fit = fitParaSurv(time=data$time,status=data$status,dist="log-normal",tau=c(5,10,25));
show(fit);
```

```
## Fitted log-normal Distribution. 
## Estimated Parameters:
##     Aspect Estimate     SE     L    U
## 1 Location     1.04 0.0676 0.907 1.17
## 2    Scale     2.08 0.0514 1.980 2.18
## 
## Distribution Properties:
##     Aspect Estimate       SE       L       U
## 1     Mean    24.60 3.24e+00   18.30    31.0
## 2   Median     2.83 1.91e-01    2.45     3.2
## 3 Variance 45400.00 2.10e+04 4190.00 86700.0
## 
## Restricted Mean Survival Times:
##   Tau Estimate     SE    L    U
## 1   5     2.83 0.0587 2.72 2.95
## 2  10     4.45 0.1190 4.22 4.68
## 3  25     7.40 0.2790 6.85 7.94
```

## Weibull

The Weibull distribution is parameterized in terms of the shape $\alpha$ and rate $\lambda$. The density is:

$$
f(t) = \alpha\lambda(\lambda t)^{\alpha-1}e^{-(\lambda t)^{\alpha}},\ t>0
$$

In the following, $n=10^{3}$ Weibull event times are simulated, with shape $\alpha=2$, rate $\lambda=2$, and expected censoring proportion 30%. Generative parameters are recovered using `fitParaSurv` with `dist="weibull"`. 


```r
# Generate Weibull time to event data
data = genData(n=1e3,dist="weibull",theta=c(2,2),p=0.3);
# Estimate parameters
fit = fitParaSurv(time=data$time,status=data$status,dist="weibull");
show(fit);
```

```
## Fitted weibull Distribution. 
## Estimated Parameters:
##   Aspect Estimate     SE    L    U
## 1  Shape     1.99 0.0596 1.88 2.11
## 2   Rate     1.96 0.0374 1.89 2.03
## 
## Distribution Properties:
##     Aspect Estimate      SE      L      U
## 1     Mean   0.4520 0.00860 0.4350 0.4690
## 2   Median   0.4240 0.00848 0.4080 0.4410
## 3 Variance   0.0561 0.00370 0.0489 0.0633
```

# Group Contrasts

## Overview

The function `compParaSurv` takes the observation times `time`, status indicators `status`, and the treatment arms `arm` as inputs. `arm` is coded as 1 for the target group, and 0 for the reference group. The distributions for the target and reference groups are selected using `dist1` and `dist0`, respectively. If unspecified, the distribution for the reference group defaults to the distribution for the target group, and the distribution for the target group defaults to Weibull. The output is an object of class `contrast`. The slots `contrast@Model1` and `contrast@Model2` contain the fitted models for the target and reference groups. See [fit](#fit) for a description of class `fit` objects. The slot `contrast@Location` contains the location parameter contrasts, CIs, and p-values. If truncation times `tau` are specified, the slot `contrast@RMST` will contain these contrasts. 

By default, all CIs are p-values are asymptotic. If a number of bootstraps `boot` is specified, then bootstrap confidence intervals are constructed using that many resamples. Likewise, if a number of permutations `perm` is specified, then permutation p-values are calculated using that many reshuffles of the treatment assignments. 

## Examples

### Comparison of Distinct Gammas

The target group consists of $10^{3}$ observations from a gamma distribution with shape $\alpha=2$, rate $\lambda=1$, and 25% censoring. The reference group consists of $10^{3}$ observations from a gamma distribution with the same shape ($\alpha=2$) but different rate ($\lambda=2$) and censoring frequency (15%). The true difference of mean survival times is $1.0$, and the true ratio is $2.0$. The true difference of median survival times is $0.84$, and the true ratio is $2$. 


```r
set.seed(101);
# Target group
d1 = genData(n=1e3,dist="gamma",theta=c(2,1),p=0.25);
d1$arm = 1;
# Reference group 
d0 = genData(n=1e3,dist="gamma",theta=c(2,2),p=0.15);
d0$arm = 0;
# Overall data 
data = rbind(d1,d0);
# Comparison
comp = compParaSurv(time=data$time,status=data$status,arm=data$arm,dist1="gamma",dist0="gamma");
cat("\n");
show(comp);
```

```
## 
## Contrast of Fitted gamma Distributions. 
## 
## Fitted Characteristics for Group 1:
##     Aspect Estimate     SE    L    U
## 1     Mean     2.05 0.0553 1.94 2.15
## 2   Median     1.68 0.0458 1.60 1.77
## 3 Variance     2.31 0.1750 1.96 2.65
## 
## Fitted Characteristics for Group 0:
##     Aspect Estimate     SE     L     U
## 1     Mean    1.020 0.0240 0.970 1.060
## 2   Median    0.859 0.0209 0.818 0.899
## 3 Variance    0.502 0.0335 0.436 0.567
## 
## Location:
##      Contrast Point     SE     L     U        P
## 1 Mean1-Mean0 1.030 0.0603 0.910 1.150 3.02e-65
## 3 Mean1/Mean0 2.010 0.0722 1.870 2.160 2.25e-84
## 2   Med1-Med0 0.826 0.0503 0.728 0.925 1.19e-60
## 4   Med1/Med0 1.960 0.0715 1.830 2.110 2.03e-76
```

### Comparison of Equivalent Weibulls

Both the target and the reference groups consist of $10^{3}$ observations from the Weibull distribution with shape $\alpha=2$ and rate $\lambda=2$. However, the target group is subject to 50% censoring, while the reference group is uncensored. The true difference in means and medians is $0.0$. The true ratio of means and medians is $1.0$. Both asymptotic and bootstrap confidence intervals are constructed. 


```r
# Target group
d1 =genData(n=1e3,dist="weibull",theta=c(2,2),p=0.5);
d1$arm = 1;
# Reference group
d0 = genData(n=1e3,dist="weibull",theta=c(2,2),p=0.0);
d0$arm = 0;
# Overall data
data = rbind(d1,d0);
# Comparison
comp = compParaSurv(time=data$time,status=data$status,arm=data$arm,dist1="weibull",dist0="weibull",boot=1e3);
cat("\n");
show(comp);
```

```
## 
## Contrast of Fitted weibull Distributions. 
## 
## Fitted Characteristics for Group 1:
##     Aspect Estimate      SE      L      U
## 1     Mean   0.4280 0.00957 0.4100 0.4470
## 2   Median   0.4030 0.00887 0.3850 0.4200
## 3 Variance   0.0495 0.00411 0.0415 0.0576
## 
## Fitted Characteristics for Group 0:
##     Aspect Estimate      SE      L      U
## 1     Mean   0.4340 0.00712 0.4200 0.4480
## 2   Median   0.4080 0.00753 0.3930 0.4230
## 3 Variance   0.0508 0.00238 0.0461 0.0555
## 
## Location:
##      Contrast    Point     SE       L      U     P  L.boot U.boot
## 1 Mean1-Mean0 -0.00569 0.0119 -0.0291 0.0177 0.633 -0.0286 0.0191
## 2 Mean1/Mean0  0.98700 0.0274  0.9350 1.0400 0.634  0.9360 1.0400
## 3   Med1-Med0 -0.00540 0.0116 -0.0282 0.0174 0.643 -0.0294 0.0189
## 4   Med1/Med0  0.98700 0.0283  0.9330 1.0400 0.643  0.9300 1.0500
```

### Log-Normal v. Weibull, Different Means, Same Medians

The target group consists of $10^{3}$ observations from a log-normal distribution with location $\mu=0$ and scale $\sigma=\sqrt{2\ln2}$. The reference group consists of $10^{3}$ observations from a Weibull distribution with shape $\alpha=2$ and rate $\lambda = \sqrt{\ln(2)}$. In each case the expected censoring proportion is 10%. The true difference of mean survival times is $0.94$, and the true ratio is $1.88$. The true difference of median survival times is $0.0$, and the true ratio is $1.0$. Both asymptotic and permutation p-values are estimated. 


```r
set.seed(105);
# Target group
d1 = genData(n=1e3,dist="log-normal",theta=c(0,sqrt(2*log(2))),p=0.1);
d1$arm = 1;
# Reference group
d0 = genData(n=1e3,dist="weibull",theta=c(2,sqrt(log(2))),p=0.1);
d0$arm = 0;
# Overall data
data = rbind(d1,d0);
# Comparison
comp = compParaSurv(time=data$time,status=data$status,arm=data$arm,dist1="log-normal",dist0="weibull",perm=1e3);
cat("\n");
show(comp);
```

```
## 
## Contrast of Fitted log-normal v. weibull 
## 
## Fitted Characteristics for Group 1:
##     Aspect Estimate     SE    L     U
## 1     Mean     1.96 0.0972 1.77  2.15
## 2   Median     1.02 0.0376 0.95  1.10
## 3 Variance    10.30 1.7400 6.89 13.70
## 
## Fitted Characteristics for Group 0:
##     Aspect Estimate     SE     L     U
## 1     Mean    1.060 0.0176 1.030 1.100
## 2   Median    1.010 0.0184 0.971 1.040
## 3 Variance    0.284 0.0147 0.256 0.313
## 
## Location:
##      Contrast  Point     SE       L      U        P P.perm
## 1 Mean1-Mean0 0.9000 0.0988  0.7060 1.0900 8.40e-20  0.001
## 2 Mean1/Mean0 1.8500 0.0965  1.6700 2.0500 7.73e-32  0.001
## 3   Med1-Med0 0.0167 0.0418 -0.0653 0.0987 6.90e-01  0.996
## 4   Med1/Med0 1.0200 0.0417  0.9380 1.1000 6.89e-01  0.996
```

### Gamma v. Exponential, Same Means, Different Medians

The target group consists of $10^{3}$ observations from a gamma distribution with shape $\alpha=4$ and scale $\lambda=4$. The reference group consists of $10^{3}$ observations from an exponential distribution with rate $\lambda = 1$. In each case the expected censoring proportion is 20%. The true difference of mean survival times is $0.0$, and the true ratio is $1.0$. The true difference of median survival times is $0.22$, and the true ratio is $1.32$. Both asymptotic and permutation p-values are estimated. RMSTs at $\tau=(0.5,1.0,1.5)$ are requested. 


```r
set.seed(106);
# Target group
d1 = genData(n=1e3,dist="gamma",theta=c(4,4),p=0.2);
d1$arm = 1;
# Reference group
d0 = genData(n=1e3,dist="exp",theta=c(1),p=0.2);
d0$arm = 0;
# Overall data
data = rbind(d1,d0);
# Comparison
comp = compParaSurv(time=data$time,status=data$status,arm=data$arm,dist1="gamma",dist0="exp",tau=c(0.5,1.0,1.5),perm=1e3);
cat("\n");
show(comp);
```

```
## 
## Contrast of Fitted gamma v. exp 
## 
## Fitted Characteristics for Group 1:
##     Aspect Estimate     SE     L     U
## 1     Mean    0.993 0.0170 0.960 1.030
## 2   Median    0.913 0.0158 0.882 0.944
## 3 Variance    0.242 0.0144 0.214 0.271
## 
## Fitted Characteristics for Group 0:
##     Aspect Estimate     SE     L    U
## 1     Mean    0.984 0.0350 0.915 1.05
## 2   Median    0.682 0.0243 0.634 0.73
## 3 Variance    0.968 0.0690 0.833 1.10
## 
## Location:
##      Contrast   Point     SE       L      U        P P.perm
## 1 Mean1-Mean0 0.00928 0.0390 -0.0671 0.0856 8.12e-01  0.892
## 2 Mean1/Mean0 1.01000 0.0399  0.9340 1.0900 8.12e-01  0.890
## 3   Med1-Med0 0.23100 0.0290  0.1740 0.2880 1.50e-15  0.001
## 4   Med1/Med0 1.34000 0.0530  1.2400 1.4500 1.71e-13  0.001
## 
## RMST:
##   Tau    Contrast  Point      SE      L      U         P P.perm
## 1 0.5 RMST1-RMST0 0.0894 0.00366 0.0822 0.0966 1.90e-131  0.001
## 2 0.5 RMST1/RMST0 1.2300 0.01110 1.2100 1.2500 3.35e-115  0.001
## 3 1.0 RMST1-RMST0 0.1760 0.01200 0.1520 0.1990  3.07e-48  0.001
## 4 1.0 RMST1/RMST0 1.2800 0.02270 1.2400 1.3300  3.35e-44  0.001
## 5 1.5 RMST1-RMST0 0.1680 0.02010 0.1290 0.2080  6.13e-17  0.001
## 6 1.5 RMST1/RMST0 1.2200 0.02980 1.1600 1.2800  6.01e-16  0.001
```

