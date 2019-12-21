# BASSLINE  <img src="man/figures/sticker.svg" align="right" width="150" />

 <!-- badges: start -->
  [![Travis build status](https://travis-ci.org/nathansam/BASSLINE.svg?branch=master)](https://travis-ci.org/nathansam/BASSLINE)
   [![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/nathansam/SMLN?branch=master&svg=true)](https://ci.appveyor.com/project/nathansam/SMLN)
  [![codecov](https://codecov.io/gh/nathansam/BASSLINE/branch/master/graph/badge.svg)](https://codecov.io/gh/nathansam/SMLN)
  <!-- badges: end -->

BAyeSian Survival anaLysIs usiNg shapE mixtures of log-normal distributions.

A project to convert code written by Vallejos & Steel[1] into an R library using C++ to speed-up functions.  

## Overview

Mixtures of life distributions provide a convienient framework for survival analysis: particularly when standard models such as the Weibull or the log-normal are unable to capture some features from the data. These mixtures can also account for unobserved heterogeneity or outlying observations.  

BASSLINE uses shape mixtures of log-normal distributions to fit data with various tail behaviour. Accelerated failure time (AFT) regressions are used over proportional hazard models which provides a clearer interpretation of the regression parameters.  


## Installation

The recommended way to install BASSLINE is via pak which will handle installing all dependencies: 
```R
if (!requireNamespace("pak", quietly = TRUE))
    install.packages("pak")
pak::pkg_install("nathansam/BASSLINE")
```

Alternatively, BASSLINE can be installed via the devtools package

```R
devtools::install_github("nathansam/BASSLINE")
```

## References 
- [1] <a href="http://dx.doi.org/10.1080/01621459.2014.923316">Vallejos & Steel (2015). Journal of the American Statistical Association. </a>
