
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tidykosmic R Package

<!-- badges: start -->
<!-- badges: end -->

## Overview

The tidykosmic R package estimates reference intervals from routinely
collected laboratory data. It does this by fitting a normal distribution
to the central part of the observed data. The [original C++
library](https://gitlab.miracum.org/kosmic/kosmic) was written by Jakob
Zierk and others.

The documentation can be found on the [tidykosmic R package
website](https://www.divinenephron.co.uk/tidykosmic/).

## Installation

Tidykosmic is still under development. You can install the development
version from GitHub.

``` r
# install.packages("devtools")
devtools::install_github("divinenephron/tidykosmic")
```

## Usage

Here a reference interval is estimated for the simulated `haemoglobin`
dataset which contains randomly generated physiological results and
contamination from randomly generated pathological result. The true
reference interval for this data is 12.0-16.0 g/dL, and the estimates
from Kosmic are close.

``` r
library(tidykosmic)
k <- kosmic(haemoglobin$result, decimals = 1)
plot(k)
```

<img src="man/figures/README-example-1.png" width="100%" />

``` r
summary(k)
#> An estimated distribution of physiological results
#> with the following quantiles:
#>     2.5%    50.0%    97.5% 
#> 12.14464 13.94471 15.95217
```
