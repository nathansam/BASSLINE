# BASSLINE
BAyeSian Survival anaLysIs usiNg shapE mixtures of log-normal distributions

A project to convert code written by Vallejos & Steel[1] into a library using C++ to speed-up functions.  

### Build states


 <!-- badges: start -->
  [![Travis build status](https://travis-ci.org/nathansam/SMLN.svg?branch=rcpp)](https://travis-ci.org/nathansam/SMLN)
  [![codecov](https://codecov.io/gh/nathansam/BASSLINE/branch/rcpp/graph/badge.svg)](https://codecov.io/gh/nathansam/SMLN)
  <!-- badges: end -->

### Installation

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

### References 
- [1] <a href="http://dx.doi.org/10.1080/01621459.2014.923316">Vallejos & Steel (2015). Journal of the American Statistical Association. </a>
