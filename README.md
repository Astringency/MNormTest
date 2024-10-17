
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MNormTest

<!-- badges: start -->

[![R-CMD-check](https://github.com/Astringency/MNormTest/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Astringency/MNormTest/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/Astringency/MNormTest/graph/badge.svg)](https://app.codecov.io/gh/Astringency/MNormTest)
<!-- badges: end -->

MNormTest包的目标是提供一套用于多元正态分布参数假设检验的函数，包括检验单个均值向量、两个均值向量、多个均值向量、单个协方差矩阵、多个协方差矩阵、同时检验均值和协方差矩阵，以及检验多元正态随机向量的独立性.

The objective of MNormTest is to provide a set of functions for
hypothesis testing of the parameters of multivariate normal
distributions, including the testing of a single mean vector, two mean
vectors, multiple mean vectors, a single covariance matrix, multiple
covariance matrices, a mean and a covariance matrix simultaneously, and
the testing of independence of multivariate normal random vectors.

## 下载 (Installation)

- CRAN:

首先，您可以从[CRAN](https://CRAN.R-project.org)安装MNormTest的发布版本:

Firstly, you can install the released version of MNormTest from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("MNormTest")
```

- GitHub:

其次，您可以从[GitHub](https://github.com/)安装MNormTest的开发版本:

What’s more, You can install the development version of MNormTest from
[GitHub](https://github.com/) with:

``` r
install.packages("pak")
pak::pak("Astringency/MNormTest")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
# library(MNormTest)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
