
<!-- README.md is generated from README.Rmd. Please edit that file -->

# InterOpt

<!-- badges: start -->

<!-- badges: end -->

The goal of InterOpt is to provide optimal weights for a weighted mean
of multiple internal controls in qPCR experiments.

## Installation

You can install the development version from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("asalimih/InterOpt")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(InterOpt)
## basic example code
x = matrix(1:6, 2, 3)
w = calcWeight(x, ctVal = TRUE, weight_method = 'geom_sd')
```
