---
title: "README"
author: "Zachary McCaw"
date: "2018-08-18"
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

This package performs estimation and inference on parametric survival curves. Supported distributions include the exponential, gamma, log-logistic, log-normal, and Weibull. Data are expected in time-to-event format, with a status indicator to accommodate non-informative right censoring. The function `fitParaSurv` provides maximum likelihood point estimates (MLEs) and confidence intervals (CIs) for the distribution of interest. Estimates are presented for model parameters, and for characteristics of the event time distribution (the mean, median, and variance). The function `compParaSurv` compares the fitted survival distributions of two treatment groups. The groups are contrasted using the difference and ratio of means and medians. For each contrast, point estimates and CIs are presented. *p*-values are calculated assessing the null hypothesis of no difference between the treatment groups. 

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

Exponential event times may be simulated using the `rWeibull` function with the shape parameter $\alpha=1$. In the following, $n=10^{3}$ exponential event times are simulated, with rate $\lambda=2$ and expected censoring proportion $20\%$. Generative parameters are recovered using `fitParaSurv` with `dist="exp"`. 


```r
# Generate exponential time to event data
D = rWeibull(n=1e3,a=1,l=2,p=0.2);
# Estimate rate
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

In the following, $n=10^{3}$ gamma event times are simulated, with shape $\alpha=2$, rate $\lambda=2$, and expected censoring proportion $20\%$. Generative parameters are recovered using `fitParaSurv` with `dist="gamma"`. 


```r
# Generate gamma time to event data
D = rGamma(n=1e3,a=2,l=2,p=0.2);
# Estimate rate
M = fitParaSurv(time=D$time,status=D$status,dist="gamma");
show(M);
```

```
## 4 update(s) performed before tolerance limit. 
## Fitted Gamma Distribution. 
## Estimated Parameters:
##   Aspect Estimate     SE    L    U
## 1  Shape     1.94 0.0869 1.77 2.11
## 2   Rate     1.93 0.1050 1.73 2.14
## 
## Distribution Properties:
##     Aspect Estimate     SE     L     U
## 1     Mean    1.000 0.0254 0.952 1.050
## 2   Median    0.835 0.0215 0.793 0.878
## 3 Variance    0.518 0.0372 0.445 0.591
```

## Log-Logistic

The log-logistic distribution is parameterized in terms of the shape $\alpha$ and rate $\lambda$. The density is:

$$
f(t) = \frac{\alpha\lambda(\lambda t)^{\alpha-1}}{[1+(\lambda t)^{\alpha}]^{2}},\ t>0
$$

In the following, $n=10^{3}$ log-logistic event times are simulated, with shape $\alpha=4$, rate $\lambda=2$, and expected censoring proportion $20\%$. Generative parameters are recovered using `fitParaSurv` with `dist="log-logistic"`. 


```r
# Generate log-logistic time to event data
D = rLogLogistic(n=1e3,a=4,l=2,p=0.2);
# Estimate rate
M = fitParaSurv(time=D$time,status=D$status,dist="log-logistic");
show(M);
```

```
## 3 update(s) performed before tolerance limit. 
## Fitted Log-Logistic Distribution. 
## Estimated Parameters:
##   Aspect Estimate     SE    L    U
## 1  Shape     4.24 0.1240 4.00 4.48
## 2   Rate     2.00 0.0266 1.95 2.05
## 
## Distribution Properties:
##     Aspect Estimate      SE      L      U
## 1     Mean   0.5490 0.00808 0.5330 0.5650
## 2   Median   0.5000 0.00665 0.4870 0.5130
## 3 Variance   0.0706 0.00648 0.0579 0.0833
```

## Log-Normal

The log-normal distribution is parameterized in terms of the location $\mu$ and scale $\sigma$. The density is:

$$
f(t) = \frac{1}{t\sigma\sqrt{2\pi}}e^{-\frac{(\ln t-\mu)^2}{2\sigma^2}},\ t>0
$$

In the following, $n=10^{3}$ log-normal event times are simulated, with location $\mu=1$, scale $\sigma=2$, and expected censoring proportion $20\%$. Generative parameters are recovered using `fitParaSurv` with `dist="log-normal"`. 


```r
# Generate log-normal time to event data
D = rLogNormal(n=1e3,m=1,s=2,p=0.2)
# Estimate rate
M = fitParaSurv(time=D$time,status=D$status,dist="log-normal");
show(M);
```

```
## 4 update(s) performed before tolerance limit. 
## Fitted Log-Normal Distribution. 
## Estimated Parameters:
##     Aspect Estimate     SE     L    U
## 1 Location    0.954 0.0663 0.824 1.08
## 2    Scale    2.010 0.0509 1.910 2.11
## 
## Distribution Properties:
##     Aspect Estimate      SE        L        U
## 1     Mean     19.7    2.54 1.48e+01    24.70
## 2   Median      2.6    1.31 2.73e-02     5.16
## 3 Variance  22200.0 9970.00 2.65e+03 41700.00
```

## Weibull

The Weibull distribution is parameterized in terms of the shape $\alpha$ and rate $\lambda$. The density is:

$$
f(t) = \alpha\lambda(\lambda t)^{\alpha-1}e^{-(\lambda t)^{\alpha}},\ t>0
$$

In the following, $n=10^{3}$ Weibull event times are simulated, with shape $\alpha=2$, rate $\lambda=2$, and expected censoring proportion $20\%$. Generative parameters are recovered using `fitParaSurv` with `dist=""`. 


```r
# Generate Weibull time to event data
D = rWeibull(n=1e3,a=2,l=2,p=0.2);
# Estimate rate
M = fitParaSurv(time=D$time,status=D$status,dist="weibull");
show(M);
```

```
## Fitted Weibull Distribution. 
## Estimated Parameters:
##   Aspect Estimate     SE    L    U
## 1  Shape     1.91 0.0531 1.81 2.02
## 2   Rate     2.01 0.0374 1.93 2.08
## 
## Distribution Properties:
##     Aspect Estimate      SE      L      U
## 1     Mean   0.4420 0.00815 0.4260 0.4580
## 2   Median   0.4110 0.00828 0.3950 0.4280
## 3 Variance   0.0578 0.00341 0.0511 0.0645
```

# Group Contrasts

## Overview

The function `compParaSurv` takes observation time `time`, the status indicator `status`, and the treatment group `arm` as inputs. `arm` is coded as 1 for the target group, and 0 for the reference group. The distributions for the target and reference groups are selected using `dist1` and `dist0`, respectively. By default, the Weibull is adopted for both groups. The output is an object of class `contrast`. The slots `contrast@Model1` and `contrast@Model2` contain the fitted models for the target and reference groups. See [fit](#fit) for a description of class `fit` objects. The slot `contrast@Contrast` contains the estimated contrasts, CIs, and *p*-values. 

## Examples

### Comparison of Distinct Gammas

In the first example, the target group consists of $10^{2}$ observations from the gamma distribution, with shape $\alpha=2$ and rate $\lambda=1$. The reference group consists of $10^{2}$ observations from the gamma distribution, with the same shape ($\alpha=2$) but rate $\lambda=2$. In each case the expected censoring proportion is $20\%$. The true difference of mean survival times is $0.5$, and the true ratio is $2.0$. The true difference of median survival times is $0.42$, and the true ratio is $2$. 


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
## 4 update(s) performed before tolerance limit. 
## 5 update(s) performed before tolerance limit. 
## 
## Contrast of Fitted Gamma Distributions. 
## 
## Fitted Characteristics for Group 1:
##     Aspect Estimate     SE     L     U
## 1     Mean    0.926 0.0808 0.768 1.080
## 2   Median    0.763 0.0666 0.633 0.894
## 3 Variance    0.471 0.1160 0.244 0.698
## 
## Fitted Characteristics for Group 0:
##     Aspect Estimate     SE     L     U
## 1     Mean    0.535 0.0390 0.459 0.612
## 2   Median    0.466 0.0342 0.399 0.533
## 3 Variance    0.114 0.0252 0.065 0.164
## 
## Contrasts:
##      Contrast Point     SE     L     U        p
## 1 Mean1-Mean0 0.391 0.0897 0.215 0.567 1.32e-05
## 2   Med1-Med0 0.297 0.0748 0.150 0.444 7.23e-05
## 3 Mean1/Mean0 1.730 0.1970 1.380 2.160 1.41e-06
## 4   Med1/Med0 1.640 0.1870 1.310 2.050 1.53e-05
```

### Comparison of Equivalent Weibulls

In this example, both the target and reference groups each consist of $10^{2}$ observations from the Weibull distribution with shape $\alpha=2$ and rate $\lambda=2$. However, the target group is subject to $50\%$ censoring, while the reference group is uncensored. The true difference in means and medians is $0.0$. The true ratio of means and medians is $1.0$. 


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
## 1     Mean   0.4310 0.0328 0.3670 0.4950
## 2   Median   0.3990 0.0295 0.3410 0.4570
## 3 Variance   0.0571 0.0153 0.0271 0.0872
## 
## Fitted Characteristics for Group 0:
##     Aspect Estimate      SE      L      U
## 1     Mean   0.4000 0.02210 0.3570 0.4430
## 2   Median   0.3690 0.02340 0.3230 0.4150
## 3 Variance   0.0505 0.00784 0.0351 0.0659
## 
## Contrasts:
##      Contrast Point     SE       L     U     p
## 1 Mean1-Mean0 0.031 0.0395 -0.0464 0.108 0.432
## 2   Med1-Med0 0.030 0.0376 -0.0437 0.104 0.425
## 3 Mean1/Mean0 1.080 0.1010  0.8960 1.300 0.426
## 4   Med1/Med0 1.080 0.1050  0.8940 1.310 0.422
```

### Log-Logistic v. Gamma, Same Means, Different Medians

In this example, the target group consists of $10^{2}$ observations from the log-logistic distribution with shape $\alpha=4$ and rate $\lambda=\pi/(2\sqrt{2})$. The reference group consists of $10^{2}$ observations from the gamma distribution with shape $\alpha=1/2$ and rate $\lambda=1/2$. In each case the expected censoring proportion is $20\%$. The true difference of mean survival times is $0.0$, and the true ratio is $1.0$. The true difference of median survival times is $0.45$, and the true ratio is $1.98$. 


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
## 5 update(s) performed before tolerance limit. 
## 3 update(s) performed before tolerance limit. 
## 
## Contrast of Fitted Log-Logistic v. Gamma 
## 
## Fitted Characteristics for Group 1:
##     Aspect Estimate     SE     L     U
## 1     Mean    1.030 0.0543 0.928 1.140
## 2   Median    0.930 0.0428 0.846 1.010
## 3 Variance    0.297 0.0956 0.110 0.485
## 
## Fitted Characteristics for Group 0:
##     Aspect Estimate     SE     L     U
## 1     Mean    0.950 0.1480 0.660 1.240
## 2   Median    0.459 0.0796 0.303 0.615
## 3 Variance    1.680 0.5910 0.527 2.840
## 
## Contrasts:
##      Contrast  Point     SE      L     U        p
## 1 Mean1-Mean0 0.0842 0.1580 -0.225 0.393 5.93e-01
## 2   Med1-Med0 0.4710 0.0904  0.294 0.648 1.85e-07
## 3 Mean1/Mean0 1.0900 0.1790  0.789 1.500 6.06e-01
## 4   Med1/Med0 2.0300 0.3640  1.430 2.880 8.27e-05
```

### Log-Normal v. Exponential, Different Means, Same Medians

In this example, the target group consists of $10^{2}$ observations from the log-normal distribution with location $\mu=0$ and scale $\sigma=\sqrt{2\ln2}$. The reference group consists of $10^{2}$ observations from the exponential distribution with rate $\lambda=1$. In each case the expected censoring proportion is $20\%$. The true difference of mean survival times is $1.0$, and the true ratio is $2.0$. The true difference of median survival times is $0.0$, and the true ratio is $1.0$. 


```r
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
## 3 update(s) performed before tolerance limit. 
## 
## Contrast of Fitted Log-Normal v. Exponential 
## 
## Fitted Characteristics for Group 1:
##     Aspect Estimate    SE      L     U
## 1     Mean    1.580 0.239  1.110  2.04
## 2   Median    0.899 0.175  0.557  1.24
## 3 Variance    5.160 2.700 -0.136 10.50
## 
## Fitted Characteristics for Group 0:
##     Aspect Estimate     SE     L     U
## 1     Mean    1.010 0.1080 0.801 1.220
## 2   Median    0.702 0.0748 0.555 0.848
## 3 Variance    1.020 0.2180 0.596 1.450
## 
## Contrasts:
##      Contrast Point    SE       L    U      p
## 1 Mean1-Mean0 0.565 0.262  0.0516 1.08 0.0310
## 2   Med1-Med0 0.198 0.190 -0.1750 0.57 0.2980
## 3 Mean1/Mean0 1.560 0.288  1.0800 2.24 0.0166
## 4   Med1/Med0 1.280 0.284  0.8300 1.98 0.2630
```
