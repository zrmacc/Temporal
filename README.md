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
df1 <- GenData(n = n, dist = "weibull", theta = c(1, 2), p = 0.20)
df1$arm <- 1
fit <- FitParaSurv(df1, dist = "weibull")
show(fit)
```

```
## Fitted Weibull Distribution. 
## Estimated Parameters:
##   Aspect Estimate    SE     L     U
## 1  Shape    0.968 0.026 0.917 1.019
## 2   Rate    2.173 0.080 2.017 2.329
## 
## Distribution Properties:
##     Aspect Estimate    SE     L     U
## 1     Mean    0.467 0.017 0.434 0.500
## 2   Median    0.315 0.013 0.291 0.340
## 3 Variance    0.233 0.023 0.189 0.277
```

```r
# Weibull data with shape = 1, rate = 2, and 25% censoring.
df0 <- GenData(n = n, dist = "weibull", theta = c(2, 2), p = 0.25)
df0$arm <- 0
fit <- FitParaSurv(df0, dist = "weibull")
show(fit)
```

```
## Fitted Weibull Distribution. 
## Estimated Parameters:
##   Aspect Estimate    SE     L     U
## 1  Shape    2.066 0.059 1.950 2.182
## 2   Rate    1.997 0.035 1.928 2.067
## 
## Distribution Properties:
##     Aspect Estimate    SE     L     U
## 1     Mean    0.444 0.008 0.428 0.459
## 2   Median    0.419 0.008 0.404 0.435
## 3 Variance    0.051 0.003 0.045 0.057
```

```r
data <- rbind(df1, df0)
comp <- CompParaSurv(data, dist1 = "weibull")
show(comp)
```

```
## Contrast of Fitted Weibull Distributions. 
## 
## Fitted Characteristics for Group 1:
##     Aspect Estimate    SE     L     U
## 1     Mean    0.467 0.017 0.434 0.500
## 2   Median    0.315 0.013 0.291 0.340
## 3 Variance    0.233 0.023 0.189 0.277
## 
## Fitted Characteristics for Group 0:
##     Aspect Estimate    SE     L     U
## 1     Mean    0.444 0.008 0.428 0.459
## 2   Median    0.419 0.008 0.404 0.435
## 3 Variance    0.051 0.003 0.045 0.057
## 
## Location:
##      Contrast  Point    SE      L      U     P
## 1 Mean1-Mean0  0.023 0.019 -0.013  0.060 0.212
## 2 Mean1/Mean0  1.053 0.043  0.972  1.140 0.204
## 3   Med1-Med0 -0.104 0.015 -0.133 -0.075 0.000
## 4   Med1/Med0  0.752 0.033  0.689  0.819 0.000
```
