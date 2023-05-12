
<!-- README.md is generated from README.Rmd. Please edit that file -->

# chemdeg

<!-- badges: start -->
<!-- badges: end -->

The goal of chemdeg is to assist chemical degradation kinetic data
analysis. It provides a rigorous statistical method for the
determination of the reaction order for n<sup>th</sup>-order kinetic
models, which is often prone to subjectivity. The first-order
multi-target model (FOMT) is introduced as a genuine alternative to
Weibull distribution. In addition, some goodness-of-fit statistics have
been added for nonlinear models.

## Installation

You can install the development version of chemdeg like so:

``` r
install.packages("chemdeg")
```

## Usage

``` r
library(chemdeg)
```

#### Determination of the reaction order

The determination of the reaction order of a generic degradation kinetic
model is often prone to subjectivity and can be misleading. By analyzing
the phase space of the dynamical system it is possible to rigorously
determine the reaction order thought statistics. In the package
`chemdeg` this analysis can be performed with the function
`det_order()`:

``` r
res <- det_order(ord1)
#> Reaction order estimated: 1
```

where `ord1` is built in the package data-frame from a
1<sup>st</sup>-order kinetic model.

The summary of the results can be accessed with the function results:

``` r
results(res)
#> 
#> Linear regression in the phase space: 
#> log(dx/dt)= 0.97 log(x) + ( -0.46 )
#> 
#> Estimate of n:
#> 
#>     Estimate   Std. Error      t value     Pr(>|t|) 
#> 9.701806e-01 4.235714e-03 2.290477e+02 1.835114e-07 
#> 
#> Confidence interval of n: 
#>     2.5 %    97.5 % 
#> 0.9567007 0.9836605 
#> 
#> Statistical analysis indicates that an order 1 degradation kineitc model is likely to describe the data.
#> The null hypothesis H0:
#> "The process is described by an order 1kinetic model"
#>  cannot be rejected.
#> 
#> Non-linear least squares regression was performed with an order  1  kinetic model:
#>  
#>  Estimate of k: 
#>    Estimate  Std. Error  t value     Pr(>|t|)
#> k 0.6878765 0.002626997 261.8489 1.541636e-11
#> In attesa che venga eseguita la profilazione...
#> Confidence interval of k: 
#>      2.5%     97.5% 
#> 0.6812276 0.6947116 
#> 
#> Goodness-of-fit:
#>                  Value
#> AIC:        -11.933893
#> AICc:       -10.933893
#> BIC:        -12.350374
#> RMSE:         2.002812
#> Chi-sq_red:   5.379602
#> -----------------------------------------------------
```

and the plot of both the phase space and the kinetic degradation data
with relative regressions can be accessed with `plot_ord()`:

``` r
plot_ord(res)
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" /><img src="man/figures/README-unnamed-chunk-4-2.png" width="100%" />
\#### First-order multi-target model fitting To fit the first-order
multi-target model, use the function `FOMT()`.

``` r
fit_fomt<-FOMT(fomtdata)
summary(fit_fomt)
#> 
#> Formula: y ~ 1 - (1 - exp(-k * t))^n
#> 
#> Parameters:
#>   Estimate Std. Error t value Pr(>|t|)    
#> k 0.056836   0.008206   6.926 0.000449 ***
#> n 9.478174   3.926148   2.414 0.052280 .  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error: 0.04243 on 6 degrees of freedom
#> 
#> Number of iterations to convergence: 10 
#> Achieved convergence tolerance: 3.955e-06
```

where `fomtdata` is an example data-frame provided with the package.

Alternatively, if the function does not converge the function `nls` from
package stats with the function `FOMTm()` as RHS of the formula and
initial values inserted arbitrarily:

``` r
fit_fomt1<-nls(y ~ FOMTm(t, k, n),
               data = list(y = fomtdata$tCQA_AA, t = fomtdata$time_h),
               start = list(n = 10, k = 0.05))
summary(fit_fomt1)
#> 
#> Formula: y ~ FOMTm(t, k, n)
#> 
#> Parameters:
#>   Estimate Std. Error t value Pr(>|t|)    
#> n 9.478164   3.926142   2.414 0.052280 .  
#> k 0.056836   0.008206   6.926 0.000449 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error: 0.04243 on 6 degrees of freedom
#> 
#> Number of iterations to convergence: 9 
#> Achieved convergence tolerance: 5.33e-06
```

#### Goodness-of-fit statistics

`chemdeg` implements the chi-squared reduced ($\chi^2_{red}=\chi^2/df$
where $df$ are the degrees of freedom) and Akaike Information Criterion
with correction for small sample size (AICc) as goodness of fit
measures. Those measures can be accessed with the functions
`chiquad_red()`and `AICC()`,respectively.

To access a table with Bayesian Information Criterion (BIC), Akaike
Information Criterion (AIC), AICc, Root Means Squared Error (RMSE) and
$\chi^2_{red}$ (from both package `stats` and `chemdeg`) call the
function `goodness_of_fit()`:

``` r
goodness_of_fit(fit_fomt)
#>                    Value
#> AIC:        -24.15681095
#> AICc:       -21.75681095
#> BIC:        -23.91848633
#> RMSE:         0.04243006
#> Chi-sq_red:           NA
```

For more details see the vignette [chemdeg
basics](articles/chemdeg_basics.html).
