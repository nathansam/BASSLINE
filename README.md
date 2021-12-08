# BASSLINE  <img src="man/figures/sticker.svg" align="right" width="150" />

<!-- badges: start -->
[![R build status](https://github.com/nathansam/BASSLINE/workflows/R-CMD-check/badge.svg)](https://github.com/nathansam/BASSLINE/actions)
[![codecov](https://codecov.io/gh/nathansam/BASSLINE/branch/master/graph/badge.svg)](https://codecov.io/gh/nathansam/BASSLINE)
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/BASSLINE)](https://cran.r-project.org/package=BASSLINE)
![Top language](https://img.shields.io/github/languages/top/nathansam/BASSLINE)
[![License: GPL-3](https://img.shields.io/badge/License-GPL3-green.svg)](https://opensource.org/licenses/GPL-3.0)
<!-- badges: end -->

## Overview

Mixtures of life distributions provide a convenient framework for survival
analysis: particularly when standard models such as the Weibull or the
log-normal are unable to capture some features from the data. These mixtures
can also account for unobserved heterogeneity or outlying observations.  

BASSLINE (**BA**ye**S**ian **S**urvival ana**L**ys**I**s usi**N**g shap**E**
mixtures of log-normal distributions) uses shape mixtures of log-normal 
distributions to fit data with fat tails and has been adapted from code written
by Vallejos & Steel[1]. Some of the functions have been rewritten in C++ for
increased performance.

5 distributions from the log-normal family are supported by BASSLINE:

* The log-normal distribution
* The log student's T distribution
* The log-logistic distribution
* The log-Laplace distribution
* The log-exponential power distribution

As well as MCMC (Markov chain Monte Carlo) algorithms for the 5
distributions, additional functions which allow log-marginal likelihood
estimators and deviance information  criteria to be calculated are provided.
Case deletion analysis and outlier detection are also supported.


## Installation

BASSLINE is currently not available on CRAN but can be installed via the
devtools package

```R
devtools::install_github("nathansam/BASSLINE")
```

## References 
- [1] <a href="http://dx.doi.org/10.1080/01621459.2014.923316">Vallejos & Steel (2015). Journal of the American Statistical Association. </a>
