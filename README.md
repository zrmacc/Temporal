# Temporal: Parametric Comparison of Time to Event Distributions

Zachary McCaw <br>
Updated: 2021-07-19

### Description

This package performs maximum likelihood based estimation and inference on time to event data subject to non-informative right censoring. `FitParaSurv` provides maximum likelihood estimation of model parameters and distributional characteristics, including the mean, median, variance, and restricted mean. `CompParaSurv` compares the mean, median, and restricted mean survival experiences of two treatment groups. Available distributions include the exponential, gamma, generalized gamma, log-normal, and Weibull. 

### Installation


```r
devtools::install_github(repo = 'zrmacc/Temporal')
```

### Vignette

A detailed vignette including the setting, distribution parameterizations, and usage examples is available [here](https://github.com/zrmacc/Temporal/blob/master/vignettes/Temporal.pdf).

### Compact Example

The follow example compares data from two Weibull distributions with the same shape but different rate parameters. 


```r
library(Temporal)
set.seed(100)

n <- 1000
# Weibull data with shape = 1, rate = 1, and 20% censoring.
df1 <- GenData(n = n, dist = "weibull", theta = c(1, 1), p = 0.20)
df1$arm <- 1
fit <- FitParaSurv(df1, dist = "weibull")
show(fit)
```

```
## Fitted Weibull Distribution. 
## Estimated Parameters:
##   Aspect Estimate    SE     L     U
## 1  Shape    0.968 0.026 0.917 1.019
## 2   Rate    1.086 0.040 1.008 1.164
## 
## Distribution Properties:
##     Aspect Estimate    SE     L     U
## 1     Mean    0.934 0.034 0.867 1.000
## 2   Median    0.630 0.025 0.581 0.679
## 3 Variance    0.931 0.090 0.754 1.107
```

```r
# Weibull data with shape = 1, rate = w, and 25% censoring.
df0 <- GenData(n = n, dist = "weibull", theta = c(1, 2), p = 0.25)
df0$arm <- 0
fit <- FitParaSurv(df0, dist = "weibull")
show(fit)
```

```
## Fitted Weibull Distribution. 
## Estimated Parameters:
##   Aspect Estimate    SE     L     U
## 1  Shape    1.033 0.030 0.975 1.091
## 2   Rate    1.995 0.071 1.856 2.133
## 
## Distribution Properties:
##     Aspect Estimate    SE     L     U
## 1     Mean    0.495 0.018 0.460 0.530
## 2   Median    0.352 0.013 0.325 0.378
## 3 Variance    0.230 0.023 0.184 0.275
```

```r
# Comparison of Weibulls.
data <- rbind(df1, df0)
comp <- CompParaSurv(data, dist1 = "weibull")
show(comp)
```

```
## Contrast of Fitted Weibull Distributions. 
## 
## Fitted Characteristics for Group 1:
##     Aspect Estimate    SE     L     U
## 1     Mean    0.934 0.034 0.867 1.000
## 2   Median    0.630 0.025 0.581 0.679
## 3 Variance    0.931 0.090 0.754 1.107
## 
## Fitted Characteristics for Group 0:
##     Aspect Estimate    SE     L     U
## 1     Mean    0.495 0.018 0.460 0.530
## 2   Median    0.352 0.013 0.325 0.378
## 3 Variance    0.230 0.023 0.184 0.275
## 
## Location:
##      Contrast Point    SE     L     U P
## 1 Mean1-Mean0 0.439 0.038 0.364 0.514 0
## 2 Mean1/Mean0 1.887 0.097 1.707 2.086 0
## 3   Med1-Med0 0.279 0.028 0.223 0.334 0
## 4   Med1/Med0 1.793 0.099 1.610 1.997 0
```
