---
title: "README"
author: "Zachary McCaw"
date: "2018-09-28"
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

This package performs estimation and inference on parametric survival curves. Supported distributions include the exponential, gamma, generalized gamma, log-logistic, log-normal, and Weibull. Data are expected in time-to-event format, with a status indicator to accommodate non-informative right censoring. The function `fitParaSurv` provides maximum likelihood point estimates (MLEs) and confidence intervals (CIs) for the distribution of interest. Estimates are presented for model parameters, and for characteristics of the event time distribution (the mean, median, and variance). The function `compParaSurv` compares the fitted survival distributions of two treatment groups. The groups are contrasted using the difference and ratio of means and medians. For each contrast, point estimates and CIs are presented. *p*-values are calculated assessing the null hypothesis of no difference between the treatment groups. 

## Setting

Suppose the data consist of pairs $(U_{i},\delta_{i})$, where the observation time $U_{i} = \min(T_{i},C_{i})$ is the minimum of the event time $T_{i}$ and a non-informative censoring time $C_{i}$. The status $\delta_{i}=I[T_{i}\leq C_{i}]$ is an indicator that an event was observed. The event times $T_{i}$ are assumed to follow a survival distribution $S_{T}$ parameterized by $\theta$. The distribution of censoring times $C_{i}$ is left unspecified. Estimation of $\theta$ proceeds by maximizing the right-censored log likelihood $\ell(\theta) = \sum_{i=1}^{n}\delta_{i}\ln h_{i} + \ln S_{i}$. Here $S_{i}$ and $h_{i}$ denote the survival and hazard contributions of the $i$th subject. The asymptotic covariance of the MLE $\hat{\theta}$ is estimated using the (inverse) observed information. Standard errors for functions of the MLE are derived from the $\Delta$-method.

# Estimation

## Overview {#fit}

The function `fitParaSurv` requires observation times `time` and the status indicator `status` as input. `status` is coded as 1 if an event was observed, and as 0 if censored. The distribution of interest is specified using `dist`, which defaults to Weibull. The output is an object of class `fit`. The slot `fit@Parameters` contains the parameter estimates and CIs. The slot `fit@Information` contains the observed sample information matrix. The slot `fit@Outcome` contains characteristics of the fitted event time distribution. Functions are presented below for simulating event times from the gamma, log-logistic, log-normal, and Weibull distributions. If an expected censoring proportion `p` is provided, the event time are subject to non-informative random right censoring. 

## Exponential

The exponential distribution is parameterized in terms of the rate $\lambda$. The density is:

$$
f(t) = \lambda e^{-\lambda t},\ t>0
$$

Exponential event times may be simulated using the `rWeibull` function with the shape parameter $\alpha=1$. In the following, $n=10^{3}$ exponential event times are simulated, with rate $\lambda=2$ and expected censoring proportion $20$%. Generative parameters are recovered using `fitParaSurv` with `dist="exp"`. 


```r
# Generate exponential time to event data
D = rWeibull(n=1e3,a=1,l=2,p=0.2);
# Estimate parameters
M = fitParaSurv(time=D$time,status=D$status,dist="exp");
show(M);
```

```
## Fitted Exponential Distribution. 
## Estimated Parameters:
##   Aspect Estimate     SE    L    U
## 1   Rate     1.92 0.0687 1.78 2.05
## 
## Distribution Properties:
##     Aspect Estimate     SE     L     U
## 1     Mean    0.522 0.0187 0.485 0.558
## 2   Median    0.361 0.0129 0.336 0.387
## 3 Variance    0.272 0.0195 0.234 0.310
```

## Gamma 

The gamma distribution is parameterized in terms of the shape $\alpha$ and rate $\lambda$. The density is:

$$
f(t) = \frac{\lambda}{\Gamma(\alpha)}(\lambda t)^{\alpha-1}e^{-\lambda t},\ t>0
$$

In the following, $n=10^{3}$ gamma event times are simulated, with shape $\alpha=2$, rate $\lambda=2$, and expected censoring proportion $20$%. Generative parameters are recovered using `fitParaSurv` with `dist="gamma"`. 


```r
# Generate gamma time to event data
D = rGamma(n=1e3,a=2,l=2,p=0.2);
# Estimate parameters
M = fitParaSurv(time=D$time,status=D$status,dist="gamma");
show(M);
```

```
## Fitted Gamma Distribution. 
## Estimated Parameters:
##   Aspect Estimate     SE    L    U
## 1  Shape     1.94 0.0869 1.77 2.11
## 2   Rate     1.93 0.1050 1.74 2.15
## 
## Distribution Properties:
##     Aspect Estimate     SE     L     U
## 1     Mean    1.000 0.0254 0.952 1.050
## 2   Median    0.835 0.0215 0.793 0.878
## 3 Variance    0.518 0.0372 0.445 0.591
```

## Generalized Gamma

The generalized gamma distributionis parameterized in terms of shapes $\alpha$ and $\beta$, and rate $\lambda$. The density is:

$$
f(t) = \frac{\beta\lambda}{\Gamma(\alpha)}(\lambda t)^{\alpha\beta-1}e^{-(\lambda t)^{\beta}},\ t>0
$$

The standard gamma and Weibull distributions are nested within the generalized gamma. Setting $\beta=1$ recovers the standard gamma, while setting $\alpha=1$ recovers the Weibull. In the following, $n=10^{3}$ generalized gamma event times are simulated, with shapes $\alpha=2$ and $\beta=2$, rate $\lambda=2$, and expected censoring proportion $20$%. Generative parameters are recovered using `fitParaSurv` with `dist="gengamma"`.


```r
set.seed(100);
# Generate generalized gamma time to event data
D = rGenGamma(n=1e4,a=2,b=2,l=2,p=0.2);
# Estimate parameters
M = fitParaSurv(time=D$time,status=D$status,dist="gengamma",report=T);
show(M);
```

```
## Objective increment:  90.2 
## Objective increment:  0.827 
## Objective increment:  0.103 
## Objective increment:  0.00407 
## Objective increment:  9.08e-06 
## Objective increment:  4e-11 
## 5 update(s) performed before tolerance limit.
## 
## Fitted GenGamma Distribution. 
## Estimated Parameters:
##   Aspect Estimate    SE    L    U
## 1 ShapeA     2.03 0.156 1.75 2.36
## 2 ShapeB     1.99 0.089 1.82 2.17
## 3   Rate     2.04 0.120 1.81 2.29
## 
## Distribution Properties:
##     Aspect Estimate      SE     L      U
## 1     Mean   0.6590 0.00256 0.654 0.6640
## 2   Median   0.6430 0.00279 0.637 0.6480
## 3 Variance   0.0569 0.00096 0.055 0.0588
```

The final parameter estimates for the generalized gamma distribution are sensitive to the initial values. If the Newton-Raphson iteration is initialized too far from the optimum of the log-likelihood, the search may not reach the maximum. If the search halts prematurely, the fitting procedure will indicate that the information matrix was not positive definite, and will return robust standard errors. Supplying initial parameter values may improve the final estimates.


```r
set.seed(103);
# Generate generalized gamma time to event data
D = rGenGamma(n=1e3,a=2,b=2,l=2,p=0.2);
# Estimate parameters
fitParaSurv(time=D$time,status=D$status,dist="gengamma",report=T);
```

```
## Objective increment:  9.56 
## 1 update(s) performed before tolerance limit.
```

```
## Warning in fit.GenGamma(time = time, status = status, sig = sig, init =
## init, : Observed information was not positive definite. Consider another
## parameter initialization.
```

```
## Eigenvalues of log-scale information:
## 11046.41 1956.48 1.76 
## Eigenvalues of original-scale information:
## 2290.04 470.49 -1.36
```

```
## Fitted GenGamma Distribution. 
## Estimated Parameters:
##   Aspect Estimate    SE     L     U
## 1 ShapeA     2.43 0.447 0.507 11.60
## 2 ShapeB     1.82 0.212 0.749  4.40
## 3   Rate     2.33 0.392 0.582  9.34
## 
## Distribution Properties:
##     Aspect Estimate       SE      L      U
## 1     Mean   0.6640 0.000725 0.6630 0.6660
## 2   Median   0.6460 0.011200 0.6240 0.6670
## 3 Variance   0.0572 0.002070 0.0531 0.0612
```

```r
# Initialization
init0 = list("la"=log(2),"lb"=log(2),"ll"=log(2));
fitParaSurv(time=D$time,status=D$status,dist="gengamma",init=init0,report=T);
```

```
## Objective increment:  0.0245 
## Objective increment:  4.45e-05 
## Objective increment:  6.96e-09 
## 2 update(s) performed before tolerance limit.
```

```
## Fitted GenGamma Distribution. 
## Estimated Parameters:
##   Aspect Estimate    SE    L    U
## 1 ShapeA     1.98 0.476 1.23 3.17
## 2 ShapeB     2.02 0.286 1.54 2.67
## 3   Rate     1.98 0.359 1.39 2.82
## 
## Distribution Properties:
##     Aspect Estimate      SE      L      U
## 1     Mean   0.6650 0.00814 0.6490 0.6810
## 2   Median   0.6490 0.00888 0.6310 0.6660
## 3 Variance   0.0576 0.00307 0.0516 0.0637
```

## Log-Logistic

The log-logistic distribution is parameterized in terms of the shape $\alpha$ and rate $\lambda$. The density is:

$$
f(t) = \frac{\alpha\lambda(\lambda t)^{\alpha-1}}{[1+(\lambda t)^{\alpha}]^{2}},\ t>0
$$

In the following, $n=10^{3}$ log-logistic event times are simulated, with shape $\alpha=4$, rate $\lambda=2$, and expected censoring proportion $20$%. Generative parameters are recovered using `fitParaSurv` with `dist="log-logistic"`. 


```r
# Generate log-logistic time to event data
D = rLogLogistic(n=1e3,a=4,l=2,p=0.2);
# Estimate parameters
M = fitParaSurv(time=D$time,status=D$status,dist="log-logistic");
show(M);
```

```
## Fitted Log-Logistic Distribution. 
## Estimated Parameters:
##   Aspect Estimate     SE    L    U
## 1  Shape     4.02 0.1160 3.79 4.25
## 2   Rate     1.99 0.0281 1.93 2.04
## 
## Distribution Properties:
##     Aspect Estimate      SE      L     U
## 1     Mean   0.5580 0.00878 0.5410 0.576
## 2   Median   0.5030 0.00710 0.4890 0.517
## 3 Variance   0.0843 0.00800 0.0686 0.100
```

## Log-Normal

The log-normal distribution is parameterized in terms of the location $\mu$ and scale $\sigma$. The density is:

$$
f(t) = \frac{1}{t\sigma\sqrt{2\pi}}e^{-\frac{(\ln t-\mu)^2}{2\sigma^2}},\ t>0
$$

In the following, $n=10^{3}$ log-normal event times are simulated, with location $\mu=1$, scale $\sigma=2$, and expected censoring proportion $20$%. Generative parameters are recovered using `fitParaSurv` with `dist="log-normal"`. 


```r
# Generate log-normal time to event data
D = rLogNormal(n=1e3,m=1,s=2,p=0.2)
# Estimate parameters
M = fitParaSurv(time=D$time,status=D$status,dist="log-normal");
show(M);
```

```
## Fitted Log-Normal Distribution. 
## Estimated Parameters:
##     Aspect Estimate     SE    L    U
## 1 Location    0.999 0.0656 0.87 1.13
## 2    Scale    1.990 0.0504 1.89 2.09
## 
## Distribution Properties:
##     Aspect Estimate       SE       L        U
## 1     Mean    19.60    2.470   14.80    24.50
## 2   Median     2.72    0.178    2.37     3.06
## 3 Variance 19700.00 8660.000 2720.00 36700.00
```

## Weibull

The Weibull distribution is parameterized in terms of the shape $\alpha$ and rate $\lambda$. The density is:

$$
f(t) = \alpha\lambda(\lambda t)^{\alpha-1}e^{-(\lambda t)^{\alpha}},\ t>0
$$

In the following, $n=10^{3}$ Weibull event times are simulated, with shape $\alpha=2$, rate $\lambda=2$, and expected censoring proportion $20$%. Generative parameters are recovered using `fitParaSurv` with `dist="weibull"`. 


```r
# Generate Weibull time to event data
D = rWeibull(n=1e3,a=2,l=2,p=0.2);
# Estimate parameters
M = fitParaSurv(time=D$time,status=D$status,dist="weibull");
show(M);
```

```
## Fitted Weibull Distribution. 
## Estimated Parameters:
##   Aspect Estimate     SE    L    U
## 1  Shape     2.12 0.0584 2.00 2.23
## 2   Rate     1.94 0.0331 1.87 2.00
## 
## Distribution Properties:
##     Aspect Estimate      SE      L      U
## 1     Mean   0.4570 0.00781 0.4420 0.4730
## 2   Median   0.4340 0.00801 0.4190 0.4500
## 3 Variance   0.0516 0.00289 0.0459 0.0573
```

# Group Contrasts

## Overview

The function `compParaSurv` takes observation time `time`, the status indicator `status`, and the treatment group `arm` as inputs. `arm` is coded as 1 for the target group, and 0 for the reference group. The distributions for the target and reference groups are selected using `dist1` and `dist0`, respectively. By default, the Weibull is adopted for both groups. The output is an object of class `contrast`. The slots `contrast@Model1` and `contrast@Model2` contain the fitted models for the target and reference groups. See [fit](#fit) for a description of class `fit` objects. The slot `contrast@Contrast` contains the estimated contrasts, CIs, and *p*-values. 

## Examples

### Comparison of Distinct Gammas

In the first example, the target group consists of $10^{2}$ observations from the gamma distribution, with shape $\alpha=2$ and rate $\lambda=1$. The reference group consists of $10^{2}$ observations from the gamma distribution, with the same shape ($\alpha=2$) but rate $\lambda=2$. In each case the expected censoring proportion is $20$%. The true difference of mean survival times is $0.5$, and the true ratio is $2.0$. The true difference of median survival times is $0.42$, and the true ratio is $2$. 


```r
# Target group
D1 = rGamma(n=1e2,a=2,l=2,p=0.2);
D1$arm = 1;
# Reference group 
D0 = rGamma(n=1e2,a=2,l=4,p=0.2);
D0$arm = 0;
# Overall data 
D = rbind(D1,D0);
# Comparison
E = compParaSurv(time=D$time,status=D$status,arm=D$arm,dist1="gamma",dist0="gamma");
cat("\n");
show(E);
```

```
## 
## Contrast of Fitted Gamma Distributions. 
## 
## Fitted Characteristics for Group 1:
##     Aspect Estimate     SE     L     U
## 1     Mean    0.988 0.0842 0.823 1.150
## 2   Median    0.804 0.0698 0.668 0.941
## 3 Variance    0.568 0.1320 0.308 0.827
## 
## Fitted Characteristics for Group 0:
##     Aspect Estimate     SE      L     U
## 1     Mean    0.479 0.0383 0.4040 0.554
## 2   Median    0.395 0.0326 0.3310 0.459
## 3 Variance    0.126 0.0276 0.0724 0.181
## 
## Contrasts:
##      Contrast Point     SE     L     U        p
## 1 Mean1-Mean0 0.509 0.0925 0.328 0.690 3.74e-08
## 2   Med1-Med0 0.410 0.0770 0.259 0.561 1.03e-07
## 3 Mean1/Mean0 2.060 0.2410 1.640 2.590 5.87e-10
## 4   Med1/Med0 2.040 0.2440 1.610 2.580 2.78e-09
```

### Comparison of Equivalent Weibulls

In this example, both the target and reference groups each consist of $10^{2}$ observations from the Weibull distribution with shape $\alpha=2$ and rate $\lambda=2$. However, the target group is subject to $50$% censoring, while the reference group is uncensored. The true difference in means and medians is $0.0$. The true ratio of means and medians is $1.0$. 


```r
# Target group
D1 = rWeibull(n=1e2,a=2,l=2,p=0.5);
D1$arm = 1;
# Reference group
D0 = rWeibull(n=1e2,a=2,l=2,p=0.0);
D0$arm = 0;
# Overall data
D = rbind(D1,D0);
# Comparison
E = compParaSurv(time=D$time,status=D$status,arm=D$arm,dist1="weibull",dist0="weibull");
cat("\n");
show(E);
```

```
## 
## Contrast of Fitted Weibull Distributions. 
## 
## Fitted Characteristics for Group 1:
##     Aspect Estimate     SE      L      U
## 1     Mean   0.4500 0.0238 0.4030 0.4970
## 2   Median   0.4420 0.0237 0.3960 0.4890
## 3 Variance   0.0317 0.0069 0.0182 0.0452
## 
## Fitted Characteristics for Group 0:
##     Aspect Estimate      SE      L      U
## 1     Mean   0.4900 0.02360 0.4440 0.5370
## 2   Median   0.4680 0.02520 0.4190 0.5180
## 3 Variance   0.0559 0.00813 0.0399 0.0718
## 
## Contrasts:
##      Contrast   Point     SE       L      U     p
## 1 Mean1-Mean0 -0.0404 0.0335 -0.1060 0.0253 0.229
## 2   Med1-Med0 -0.0261 0.0346 -0.0938 0.0417 0.451
## 3 Mean1/Mean0  0.9180 0.0656  0.7980 1.0600 0.230
## 4   Med1/Med0  0.9440 0.0717  0.8140 1.1000 0.450
```

### Log-Logistic v. Gamma, Same Means, Different Medians

In this example, the target group consists of $10^{2}$ observations from the log-logistic distribution with shape $\alpha=4$ and rate $\lambda=\pi/(2\sqrt{2})$. The reference group consists of $10^{2}$ observations from the gamma distribution with shape $\alpha=1/2$ and rate $\lambda=1/2$. In each case the expected censoring proportion is $20$%. The true difference of mean survival times is $0.0$, and the true ratio is $1.0$. The true difference of median survival times is $0.45$, and the true ratio is $1.98$. 


```r
# Target group
D1 = rLogLogistic(n=1e2,a=4,l=pi/(2*sqrt(2)),p=0.2);
D1$arm = 1;
# Reference group
D0 = rGamma(n=1e2,a=(1/2),l=(1/2),p=0.2);
D0$arm = 0;
# Overall data
D = rbind(D1,D0);
# Comparison
E = compParaSurv(time=D$time,status=D$status,arm=D$arm,dist1="log-logistic",dist0="gamma");
cat("\n");
show(E);
```

```
## 
## Contrast of Fitted Log-Logistic v. Gamma 
## 
## Fitted Characteristics for Group 1:
##     Aspect Estimate     SE      L     U
## 1     Mean    0.978 0.0431 0.8940 1.060
## 2   Median    0.902 0.0361 0.8320 0.973
## 3 Variance    0.187 0.0521 0.0849 0.289
## 
## Fitted Characteristics for Group 0:
##     Aspect Estimate    SE     L    U
## 1     Mean    1.280 0.229 0.828 1.73
## 2   Median    0.618 0.113 0.396 0.84
## 3 Variance    3.030 1.250 0.584 5.48
## 
## Contrasts:
##      Contrast  Point    SE       L     U      p
## 1 Mean1-Mean0 -0.299 0.233 -0.7560 0.158 0.2000
## 2   Med1-Med0  0.284 0.119  0.0514 0.517 0.0168
## 3 Mean1/Mean0  0.766 0.142  0.5330 1.100 0.1490
## 4   Med1/Med0  1.460 0.274  1.0100 2.110 0.0436
```

### Log-Normal v. Exponential, Different Means, Same Medians

In this example, the target group consists of $10^{2}$ observations from the log-normal distribution with location $\mu=0$ and scale $\sigma=\sqrt{2\ln2}$. The reference group consists of $10^{2}$ observations from the exponential distribution with rate $\lambda=1$. In each case the expected censoring proportion is $20$%. The true difference of mean survival times is $1.0$, and the true ratio is $2.0$. The true difference of median survival times is $0.0$, and the true ratio is $1.0$. 


```r
set.seed(100);
# Target group
D1 = rLogNormal(n=1e2,m=0,s=sqrt(2*log(2)),p=0.2);
D1$arm = 1;
# Reference group
D0 = rWeibull(n=1e2,a=1,l=1,p=0.2);
D0$arm = 0;
# Overall data
D = rbind(D1,D0);
# Comparison
E = compParaSurv(time=D$time,status=D$status,arm=D$arm,dist1="log-normal",dist0="exp");
cat("\n");
show(E);
```

```
## 
## Contrast of Fitted Log-Normal v. Exponential 
## 
## Fitted Characteristics for Group 1:
##     Aspect Estimate    SE      L     U
## 1     Mean    1.840 0.302  1.250  2.43
## 2   Median    0.956 0.113  0.735  1.18
## 3 Variance    9.190 5.180 -0.951 19.30
## 
## Fitted Characteristics for Group 0:
##     Aspect Estimate     SE     L     U
## 1     Mean    1.030 0.1200 0.798 1.270
## 2   Median    0.716 0.0833 0.553 0.879
## 3 Variance    1.070 0.2480 0.581 1.550
## 
## Contrasts:
##      Contrast Point    SE       L     U       p
## 1 Mean1-Mean0 0.808 0.325  0.1710 1.440 0.01290
## 2   Med1-Med0 0.240 0.140 -0.0351 0.514 0.08740
## 3 Mean1/Mean0 1.780 0.358  1.2000 2.640 0.00405
## 4   Med1/Med0 1.330 0.221  0.9650 1.850 0.08150
```
